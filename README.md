# Zundel cation PIMD and PIGSTA workflows

This repository contains input files and auxiliary scripts used to run
path-integral molecular dynamics (PIMD) simulations of the Zundel cation
(H5O2+) and to perform post-processing and analysis using the PIGSTA
(path-integral generalized spatial tetrahedral analysis) framework.

The scripts and workflows provided here reflect the setup used in our
study and are not intended to function as a polished, portable software
package. In particular, the scripts may rely on system-specific paths,
environment settings, and local workflow conventions, and may therefore
require adaptation to run on different computing environments.

Nevertheless, we believe that making these files publicly available
can be useful for researchers interested in reproducing or extending
PIMD-based simulations of protonated water clusters, or in implementing
related PIGSTA-style analyses.

## Associated data

The raw and processed simulation data corresponding to this repository
are archived on Zenodo:

**DOI:** https://doi.org/10.5281/zenodo.18302429

## Scope and limitations

- This repository is provided for transparency and reproducibility.
- The scripts are not guaranteed to run out of the box on other systems.
- No claim is made that the code represents a general or optimized
  implementation of PIMD or PIGSTA methods.

Users are encouraged to treat the contents as a reference implementation
and adapt them as needed for their own computational setups.
