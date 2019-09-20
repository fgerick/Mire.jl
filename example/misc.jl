
"""
    angularmom(u,a,b,c)

Calculates z-component of angular momentum.
"""
angularmom(u,a,b,c) = int_polynomial_ellipsoid((u×r)[3],a,b,c)
angularmom(u,cmat) = int_polynomial_ellipsoid((u×r)[3],cmat)
