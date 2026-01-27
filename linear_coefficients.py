import numpy as np

def linear_coefficients(r, ext_coef, sca_coef):
    """
    Compute linear interpolation coefficients y = a*r + b
    for extinction and scattering in each spherical layer.

    Parameters
    ----------
    r : ndarray (nhkm)
        Radii at layer boundaries [km]
    ext_coef : ndarray (nhkm)
        Extinction coefficient at boundaries
    sca_coef : ndarray (nhkm)
        Scattering coefficient at boundaries

    Returns
    -------
    aext, bext : ndarray (nhkm-1)
        Linear coefficients for extinction
    asca, bsca : ndarray (nhkm-1)
        Linear coefficients for scattering
    """

    nhkm = len(r)
    nlr = nhkm - 1

    aext = np.zeros(nlr)
    bext = np.zeros(nlr)
    asca = np.zeros(nlr)
    bsca = np.zeros(nlr)

    for ilr in range(nlr):
        dr = r[ilr + 1] - r[ilr]

        aext[ilr] = (ext_coef[ilr + 1] - ext_coef[ilr]) / dr
        bext[ilr] =  ext_coef[ilr] - aext[ilr] * r[ilr]

        asca[ilr] = (sca_coef[ilr + 1] - sca_coef[ilr]) / dr
        bsca[ilr] =  sca_coef[ilr] - asca[ilr] * r[ilr]

    return aext, bext, asca, bsca