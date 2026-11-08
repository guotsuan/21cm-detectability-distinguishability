# Third-party dependencies and reference products

The MIT license in this repository applies to author-contributed code. It does not relicense third-party dependencies, station-layout products, sensitivity lookup tables, or bundled derived/reference scientific inputs.

- 21cmFAST: https://github.com/21cmfast/21cmFAST — obtain version 3.4.0 and its upstream license.
- 21cmSense: https://github.com/21cmfast/21cmSense — obtain version 3.0.0 and its upstream license.
- OSKAR: https://github.com/OxfordSKA/OSKAR — version 2.13.0; binaries/container not included.
- SKAO array configuration: https://gitlab.com/ska-telescope/ost/ska-ost-array-config — reference configuration version 5.0.0.
- SKAO sensitivity lookup: https://developer.skao.int/projects/ska-ost-senscalc/en/11.1.1/_modules/ska_ost_senscalc/low/sefd_lookup.html — original HDF5 table not included, and its historical upstream release is unknown.

The `vendor/` directories retain numerical helper modules from the authors' earlier calculation packages; they are not installations of the external packages above. The extracted 21cmFAST AstroParams source text in the original collection is not redistributed here. Synthetic catalogue entries are test sources, not an observed catalogue. Small array/station reference inputs retain provenance identifying their source.
