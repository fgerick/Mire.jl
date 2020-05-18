# Mire.jl

*Modes In Rotating Ellipsoids*

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fgerick.github.io/Mire.jl/dev/)
[![Build Status](https://travis-ci.com/fgerick/Mire.jl.svg?token=NJNkFC9qALxxCxMBhjwi&branch=master)](https://travis-ci.com/fgerick/Mire.jl)
[![codecov](https://codecov.io/gh/fgerick/Mire.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/fgerick/Mire.jl)



Solve for (magneto)hydrodynamic modes in rapidly rotating ellipsoids using a Galerkin method with basis vectors constructed from Cartesian monomials. The code supports fully 3-D, fully quasi-geostrophic (QG) and a hybrid model with QG velocities and 3-D magnetic field perturbations.

Check out the [documentation](https://fgerick.github.io/Mire.jl/dev/) for help on installation and examples (the documentation will improve with time).

## Questions

If you have trouble using this library or find errors in the examples, please don't hesitate to file an issue on GitHub or get in touch with me directly (via [mail](mailto:felix.gerick@univ-grenoble-alpes.fr) or [researchgate](https://www.researchgate.net/profile/Felix_Gerick)).

## Citation

If you use this software, please cite

F Gerick, D Jault, J Noir, and J Vidal, (2020). Pressure torque of torsional Alfvén modes acting on an ellipsoidal mantle, Geophysical Journal International, [10.1093/gji/ggaa166](https://doi.org/10.1093/gji/ggaa166)

```
@article{10.1093/gji/ggaa166,
    author = {Gerick, F and Jault, D and Noir, J and Vidal, J},
    title = "{Pressure torque of torsional Alfvén modes acting on an ellipsoidal mantle}",
    journal = {Geophysical Journal International},
    volume = {222},
    number = {1},
    pages = {338-351},
    year = {2020},
    month = {04},
    issn = {0956-540X},
    doi = {10.1093/gji/ggaa166},
    url = {https://doi.org/10.1093/gji/ggaa166},
    eprint = {https://academic.oup.com/gji/article-pdf/222/1/338/33181751/ggaa166.pdf},
}
```
