# SpecCosmo

A [Cactus Code](https://cactuscode.org/) arrangement for general relativistic
magnetohydrodynamics (GRMHD) simulations of the early universe, developed by
David Garrison (University of Houston–Clear Lake).

## Thorns

| Thorn | Purpose |
|---|---|
| `SpecGRMHD` | Core GRMHD evolution thorn for the SpecCosmo arrangement. |
| `MHD_INIT` | Initial data setup for the MHD evolution. |
| `MHD_Analysis` | Analysis routines, including spectral (FFT) diagnostics, scalar/vector spectra, and pointwise/coordinate analysis. |
| `MoL` | Method of Lines generic time integrators. Originally by Ian Hawke and the Cactus team; modified by David Garrison for use in SpecCosmo. |

## Config Files

`Config Files/` contains machine-specific Cactus configuration options
(`*.cfg`), submission/environment scripts, and the arrangement's `ThornList`
for the clusters this code has been built on (Lonestar, Ranger, Maxwell,
Atlantis, Athena, Xanadu, Sabine, macOS, generic GNU/Linux, Singularity).

## Building

These thorns are intended to be checked out into the `arrangements/` directory
of a Cactus source tree. Use `Config Files/ThornList` as the thorn list and the
appropriate `.cfg` file as the options file:

```
./simfactory/bin/sim build speccosmo --thornlist "Config Files/ThornList" \
    --optionlist "Config Files/<machine>.cfg"
```

## License

`MoL` is distributed under the GNU GPL v2 (see `MoL/COPYRIGHT`). The remaining
thorns are copyright their listed authors.

## Contact

David Garrison — <garrison@uhcl.edu>
