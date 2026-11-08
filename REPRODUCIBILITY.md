# Reproduction coverage and exclusions

## Verified here
- Four module synthetic/dry-run self-tests (no new physical simulations).
- Five selected final figures regenerated from packaged small arrays/tables.
- Same-operator output distances recomputed from the separately distributed response archive, with archived power/covariance reconstruction checks.
- Python syntax parsing and a limited credential-signature scan. These are not guarantees of complete software security or physics correctness.

See `VALIDATION.json` for actual execution results. A failed or unrun operation must not be inferred successful from this list.

## Not reproduced here
- Full 434-snapshot simulation grid and 29-band instrumental recalculation.
- Physical 21cmFAST runs, full OSKAR beam generation, and raw cube regeneration.
- All paper figures; missing plotting dependency closures are explicitly listed in `FIGURE_COVERAGE.csv`.
- Independent absolute sensitivity, complete cosmological/foreground/calibration covariance, or joint parameter inference.

## Data tiers
1. Git-tracked small inputs and CSV/NPZ derived products support the included lightweight replays.
2. The release asset `matched_response_review.zip` contains the frozen two-frequency visibility results. Its embedded historical absolute paths are provenance, not live paths to be followed.
3. Large raw coeval/lightcone cubes, beam caches and production containers are outside this release. Their access arrangements must be coordinated with the authors; this repository does not imply that these data are publicly downloadable.
4. `SKA_paper_code_20261105.zip` contained historical manuscript drafts, reviewer correspondence/instructions and thousands of duplicated nested records. These are intentionally excluded from the public code repository. They are not prerequisites for the documented lightweight replay.

Only path handling in shell entrypoints and the selected matched-figure input lookup was adapted during curation. Scientific formulae, physical parameters and reference arrays were not changed. The packaging/plot replay utility is new. Preserve the original reference hashes: they identify original inputs and should not be replaced with hashes of modified calculations.
