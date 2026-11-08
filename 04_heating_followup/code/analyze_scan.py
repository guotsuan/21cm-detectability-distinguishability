#!/usr/bin/env python3
"""Strict baseline regression before interpreting an exploratory heating scan."""
import argparse,sys
from common import *
from interpolation import interp
from likelihood import expected_2logLR_diagonal
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
from signal_core import fft_power,spherical
from run_scan import jobs
EDGES=np.array([.05,.075,.1,.125,.15,.175,.2,.225,.25,.275,.3,.325,.35])

def main(a):
    runs=a.runs/f'seed{a.seed}';out=output(a.out);inputs=[];specs=[];hist=[];signals={};reg=[];density={};density_rows=[];runtime=None
    for tag,expected in jobs(a.frequencies,4,False,a.seed):
        nu=expected['frequency_MHz'];folder=runs/'production'/f'{nu}MHz'/tag/'real'
        require((folder/'COMPLETE.json').exists(),f'Missing completed job {folder}; finish run_scan first')
        cfg=json.loads((folder/'requested.json').read_text());left=json.loads(json.dumps(expected));right=json.loads(json.dumps(cfg));left['user_params'].pop('N_THREADS');right['user_params'].pop('N_THREADS');require(left==right,'Run configuration differs: '+str(folder))
        for name in ['requested.json','COMPLETE.json','resolved.json','environment.json','accepted_parameters.json']:inputs.append(folder/name)
        path=folder/'cube.npz';done=json.loads((folder/'COMPLETE.json').read_text());require(sha(path)==done['cube_sha256'],'Cube checksum mismatch');inputs.append(path)
        env=json.loads((folder/'environment.json').read_text());require(any(str(k).lower()=='21cmfast' and str(v)=='3.4.0' for k,v in env['packages'].items()),'Expected 21cmFAST 3.4.0 runtime')
        require('loaded_modules' in env,'Missing fresh runtime provenance')
        if runtime is None:runtime=env['loaded_modules']
        else:require(runtime==env['loaded_modules'],'Loaded simulation modules differ across jobs; use one unchanged environment')
        require(not (folder/'reuse.json').exists(),'Old reused controls are not permitted in this followup')
        accepted=json.loads((folder/'accepted_parameters.json').read_text())
        for key in ['L_X','M_TURN']:require(accepted[key]['log10_input']==cfg['astro_params'][key],'Accepted input differs')
        resolved=json.loads((folder/'resolved.json').read_text())
        for key in ['L_X','M_TURN']:require(accepted[key]['resolved']==resolved['astro_params'][key],'Accepted/resolved values differ')
        require(not resolved['flag_options']['APPLY_RSDS'] and not resolved['flag_options']['SUBCELL_RSD'],'RSD must remain disabled')
        with np.load(path,allow_pickle=False) as z:
            require(z['brightness_temp'].shape==(128,128,128) and float(z['box_len_cMpc'])==150 and int(z['random_seed'])==a.seed,'Wrong production geometry/seed')
            require(abs(float(z['redshift'])-cfg['redshift'])<1e-10 and str(z['coordinate_space'])=='real','Wrong redshift/space')
            dh=hashlib.sha256(np.ascontiguousarray(z['density']).tobytes()).hexdigest()
            if nu in density:
                base=density[nu];current=np.asarray(z['density'],dtype=float)
                delta=current-base['array'];rms=float(np.sqrt(np.mean(delta**2)));std=float(np.std(base['array']))
                density_rows.append(dict(seed=a.seed,frequency_MHz=nu,tag=tag,reference='F_control',exact_match=dh==base['hash'],rms_difference=rms,relative_rms=rms/std if std else None,max_abs_difference=float(np.max(np.abs(delta))),density_sha256=dh))
            else:density[nu]=dict(hash=dh,array=np.asarray(z['density'],dtype=float))
            k,p,mean,vr,vf=fft_power(z['brightness_temp'],150);km,pm,de,count,sums=spherical(k,p,EDGES)
            require(np.all(count>0) and np.all(de>0),'Unsupported/zero spectrum')
            hist.append(dict(seed=a.seed,frequency_MHz=nu,tag=tag,L_X=cfg['astro_params']['L_X'],M_TURN=cfg['astro_params']['M_TURN'],mean_Tb_mK=mean,xHI_volume=float(z['xH_box'].mean()),mean_Tk_K=float(z['Tk_box'].mean()),mean_Ts_K=float(z['Ts_box'].mean()),density_sha256=dh,elapsed_seconds=done['elapsed_seconds']))
        t=pd.DataFrame(dict(k_Mpc_inv=km,Delta2_mK2=de,fft_cells=count))
        for row in t.to_dict('records'):specs.append(dict(seed=a.seed,frequency_MHz=nu,tag=tag,**row))
        d,npz=native(a.forecast_root,nu);inputs.append(npz);sig=interp(t,d['k_Mpc_inv']);signals[nu,tag]=sig
        if a.seed==12345 and tag in ['F_control','Mplus_LX40.00']:
            model='F' if tag=='F_control' else 'Mplus';use=(d['k_Mpc_inv']>=.1)&(d['k_Mpc_inv']<=.3)
            require(np.isfinite(sig[use]).all(),'Missing native analysis support')
            error=float(np.max(abs(sig[use]/d[model+'_signal_mK2'][use]-1)))
            reg.append(dict(seed=a.seed,frequency_MHz=nu,model=model,max_relative_native_signal_error=error,tolerance=a.regression_rtol,passed=error<=a.regression_rtol))
    pd.DataFrame(specs).to_csv(out/'spectra.csv',index=False);pd.DataFrame(hist).to_csv(out/'thermal_history.csv',index=False);pd.DataFrame(reg).to_csv(out/'control_regression.csv',index=False)
    pd.DataFrame(density_rows).to_csv(out/'density_diagnostics.csv',index=False)
    if not all(r['exact_match'] for r in density_rows):
        finish(out,inputs,['DENSITY_MISMATCH: diagnostics saved; no science metrics accepted.'])
        raise ValueError('Fresh density fields differ. Send the review archive; do not bypass this gate.')
    if not all(r['passed'] for r in reg):
        finish(out,inputs,['CONTROL_REGRESSION_FAILED: do not interpret scan. Check version/configuration/evolution anchor and support.'])
        raise ValueError('Baseline regression failed; diagnostics saved. Do not relax tolerance without finding the cause.')
    repeat=[]
    oldpath=ROOT/'evidence/reviewed140_spectra.csv'
    if a.seed==12345 and 140 in a.frequencies and oldpath.is_file():
        old=pd.read_csv(oldpath);inputs.append(oldpath)
        d,_=native(a.forecast_root,140);use=mask(d,'no_wedge')
        for tag in ['F_control','Mplus_LX40.00','Mplus_LX40.25']:
            previous=old[(old.frequency_MHz==140)&(old.tag==tag)]
            require(len(previous)>1,'Missing reviewed reference spectrum '+tag)
            oldsig=interp(previous,d['k_Mpc_inv']);fresh=signals[140,tag]
            require(np.isfinite(oldsig[use]).all(),'Reference spectrum lacks common support')
            err=float(np.max(np.abs(fresh[use]/oldsig[use]-1)))
            repeat.append(dict(frequency_MHz=140,tag=tag,max_relative_native_signal_error=err,diagnostic_tolerance=a.regression_rtol,within_diagnostic_tolerance=err<=a.regression_rtol))
        pd.DataFrame(repeat).to_csv(out/'reviewed140_repeatability.csv',index=False)
    rows=[]
    for nu in a.frequencies:
        d,_=native(a.forecast_root,nu);f=d['F_signal_mK2'];s=d['F_total_sigma_mK2']
        for lx in ([40.,40.20,40.25,40.30] if a.refine else [40.,40.25]):
            tag=f'Mplus_LX{lx:.2f}';m=signals[nu,tag]
            for case in CASES:
                use=mask(d,case);require(np.isfinite(m[use]).all(),'No extrapolation allowed')
                Af,Df,Dsf=fit(signals[nu,'F_control'][use],m[use],s[use])
                A,D,Ds=fit(f[use],m[use],s[use]);sc=d['thermal_sigma_mK2'][use]+d['sample_coefficient'][use]*m[use]
                rows.append(dict(seed=a.seed,frequency_MHz=nu,L_X=lx,case=case,D_fixed_F=D,A_profiled=A,D_shape=Ds,D_fresh_F=Df,A_fresh_F=Af,D_shape_fresh_F=Dsf,expected_2logLR_gaussian_Ftruth=expected_2logLR_diagonal(f[use],m[use],s[use]**2,sc**2),entries=int(use.sum())))
    t=pd.DataFrame(rows)
    t['D_archived_F_secondary']=t.D_fixed_F
    t['D_shape_archived_F_secondary']=t.D_shape
    t['expected_2logLR_archived_Ftruth_secondary']=t.expected_2logLR_gaussian_Ftruth
    t.to_csv(out/'scan_metrics.csv',index=False);best=[]
    for (nu,case),g in t.groupby(['frequency_MHz','case']):
        r=g.loc[g.D_fresh_F.idxmin()];baseline=g[g.L_X==40].iloc[0]
        best.append(dict(seed=a.seed,frequency_MHz=nu,case=case,selection_metric='D_fresh_F',target='same_seed_fresh_F',covariance='archived_seed12345_F',best_sampled_L_X=r.L_X,D_best_sampled_fresh=r.D_fresh_F,D_baseline_fresh=baseline.D_fresh_F,distance_ratio_fresh=r.D_fresh_F/baseline.D_fresh_F if baseline.D_fresh_F>0 else np.nan,best_on_boundary=r.L_X in [g.L_X.min(),g.L_X.max()],A_fresh_at_best=r.A_fresh_F,D_shape_fresh_at_best=r.D_shape_fresh_F,D_archived_F_secondary_at_fresh_best=r.D_fixed_F))
    transfer=t[t.L_X==40.25].copy()
    baseline=t[t.L_X==40.00][['frequency_MHz','case','D_fixed_F','D_fresh_F']].rename(columns={'D_fixed_F':'baseline_D','D_fresh_F':'baseline_D_fresh'})
    transfer=transfer.merge(baseline,on=['frequency_MHz','case'],validate='one_to_one')
    transfer['distance_ratio_archived_secondary']=transfer.D_fixed_F/transfer.baseline_D
    transfer['distance_ratio_fresh']=transfer.D_fresh_F/transfer.baseline_D_fresh
    transfer['primary_metric']='D_fresh_F'
    transfer['target']='same_seed_fresh_F'
    transfer['covariance']='archived_seed12345_F'
    transfer.to_csv(out/'fixed_candidate_transfer.csv',index=False)
    pd.DataFrame(best).to_csv(out/'best_sampled.csv',index=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
    for nu in a.frequencies:
        q=t[(t.frequency_MHz==nu)&(t.case=='horizon_buffer')].sort_values('L_X')
        ax[0].plot(q.L_X,q.D_fresh_F,'o-',label=f'{nu} MHz');ax[1].plot(q.L_X,q.D_shape_fresh_F,'o-',label=f'{nu} MHz')
    for x,title in zip(ax,['Same-seed fresh F; archived F covariance','Same target/covariance; amplitude profiled']):x.set(xlabel='21cmFAST L_X input (log10)',ylabel='Conditional distance',title=title);x.grid(alpha=.2);x.legend()
    fig.savefig(out/'heating_scan.png',dpi=180);plt.close(fig)
    finish(out,inputs,['Exploratory discrete followup, not an observationally allowed prior or posterior.','Primary outputs use same-seed fresh F target with original archived seed12345 F covariance; Mplus M_TURN=9.5 fixed, only L_X changes.','D_fixed_F is retained for compatibility as SECONDARY archived seed12345 F mean diagnostic; explicit alias D_archived_F_secondary. All plots and main summaries use D_fresh_F.','D_shape_fresh_F profiles nonnegative amplitude for the primary matched-seed target; D_shape is a secondary archived-target diagnostic.','Gaussian expected_2logLR includes trace, determinant and mean terms but uses unvalidated diagonal covariance; not exclusion significance.','No interpolation between parameter samples; boundary best is unresolved.','No joint-band inference or physically calibrated gain prior.'])
    status=dict(status='PASS_FRESH_CHECKS',seed=a.seed,primary_metric='D_fresh_F',primary_ratio='distance_ratio_fresh',target='same_seed_fresh_F',covariance='archived_seed12345_F',secondary_archived_target_columns=['D_fixed_F','D_archived_F_secondary','D_shape','D_shape_archived_F_secondary','distance_ratio_archived_secondary','expected_2logLR_gaussian_Ftruth'],archived_control_regression=('passed' if a.seed==12345 else 'not_applicable_different_realization'),reviewed140_repeatability=('not_tested' if not repeat else 'within_diagnostic_tolerance' if all(r['within_diagnostic_tolerance'] for r in repeat) else 'REVIEW_REQUIRED'),frequencies=a.frequencies,refine=a.refine,interpretation='Conditional coeval distances, not calibrated exclusions or joint-band likelihoods')
    (out/'status.json').write_text(json.dumps(status,indent=2))
    if repeat and not all(r['within_diagnostic_tolerance'] for r in repeat):
        print('REVIEW REQUIRED: old/new spectra differ beyond the diagnostic threshold. Fresh metrics are saved, but do not claim repeatability or update the paper before review.')
    print('PRIMARY: same-seed fresh F target, fixed archived seed12345 F covariance; no physical parameter re-optimization.')
    print(pd.DataFrame(best).drop(columns=['D_archived_F_secondary_at_fresh_best']).to_string(index=False))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--forecast-root',type=Path,default=ROOT/'evidence/frequency_analysis');p.add_argument('--frequencies',type=int,nargs='+',default=[140],choices=[130,140]);p.add_argument('--regression-rtol',type=float,default=3e-5);p.add_argument('--seed',type=int,choices=[12345,23456,34567],default=12345);p.set_defaults(refine=False);a=p.parse_args();require(a.regression_rtol>0,'Tolerance must be positive');main(a)
