# Mire.jl

*Modes In Rotating Ellipsoids*

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fgerick.github.io/Mire.jl/dev/)
[![](https://github.com/fgerick/Mire.jl/workflows/CI/badge.svg)](https://github.com/fgerick/Mire.jl/actions)
[![codecov](https://codecov.io/gh/fgerick/Mire.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/fgerick/Mire.jl)



Solve for inviscid (magneto)hydrodynamic modes in rapidly rotating spheres and ellipsoids using a Galerkin method with basis vectors constructed from Cartesian monomials. The code supports fully 3-D, fully quasi-geostrophic (QG) and a hybrid model with QG velocities and 3-D magnetic field perturbations. The magnetic field satisfies the perfectly conducting boundary condition in the ellipsoid. In the sphere the insulating boundary condition is also available. The velocity field satisfies the non-penetrating boundary condition, valid for inviscid fluids.

Check out the [documentation](https://fgerick.github.io/Mire.jl/dev/) for help on installation and examples (the documentation will improve with time).

## Questions

If you have trouble using this library or find errors in the examples, please don't hesitate to file an issue on GitHub or get in touch with me directly (via felix[dot]gerick[at]observatory[dot]be or [researchgate](https://www.researchgate.net/profile/Felix_Gerick)).

## Citation

If you use this software in your research, please consider citing one of the relevant articles:

Gerick F., Jault D., Noir J., and Vidal J., (2020). Pressure torque of torsional Alfvén modes acting on an ellipsoidal mantle, Geophysical Journal International, [10.1093/gji/ggaa166](https://doi.org/10.1093/gji/ggaa166)

Gerick F., Jault D., and Noir J., (2021). Fast Quasi-Geostrophic Magneto-Coriolis Modes in the Earth's Core, Geophysical Research Letters, [10.1029/2020GL090803](https://doi.org/10.1029/2020GL090803)

```
@article{gerick_pressure_2020,
    author = {Gerick, F. and Jault, D. and Noir, J. and Vidal, J.},
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

@article{gerick_fast_2021,
  title = {Fast Quasi-Geostrophic Magneto-Coriolis Modes in the Earth's Core},
  author = {Gerick, F. and Jault, D. and Noir, J.},
  year = {2021},
  journal = {Geophysical Research Letters},
  volume = {48},
  number = {4},
  pages = {e2020GL090803},
  issn = {1944-8007},
  doi = {10.1029/2020GL090803},
  langid = {english},
  annotation = {\_eprint: https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2020GL090803}
}
```
