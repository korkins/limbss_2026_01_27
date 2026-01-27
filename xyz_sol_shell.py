import numpy as np

def xyz_sol_shell(mu0, xs, ys, zs, r):
    """
    Compute the intersection of a solar beam with a sphere.

    Parameters
    ----------
    mu0 : float
        Cosine of solar zenith angle (cos(SZA)), must be > 0.
    xs, ys, zs : float
        Coordinates of a point on the line of sight irradiated by the Sun.
    Re : float
        Earth radius
    r : float
        Radius of the sphere centered in the origin

    Returns
    -------
    xp, yp, zp : float
        Coordinates of the intersection point with zp > zs.
        Returns np.nan if no intersection exists.

    Notes
    -----
    1. In the given system of coordinates, lo = {sqrt(1 - mu0**2), 0, -mu0}
       - note the minus: Z-axis points "up", the soalr beam goes "down".

    2. We consider only zp >= zs
    
    """
    
    tiny_distance_km = 0.010 # 10 meters
    
    rs = np.sqrt(xs**2 + ys**2 + zs**2)
    
    if rs > r:
        print("P = {xs, ys, zs} is outside the h-sphere:")
        print(f"Distance to the point: {rs:.3f}")
        print(f"Radius of the sphere, Re+h: {r:.3f}")
        return np.nan, np.nan, np.nan
    else:
        smu0 = np.sqrt(1.0 - mu0**2)   # sin(SZA)
    
        # Quadratic coefficients (A = 1 dropped)
        B = 2.0 * (smu0 * xs - mu0 * zs)
        C = xs**2 + ys**2 + zs**2 - r**2
        
        D = B**2 - 4.0 * C
        
        if D < 0:
            # No intersection - this cant happen in the 'else': the point belongs to the sphere
            return np.nan, np.nan, np.nan
        elif D < tiny_distance_km:
            # D = 0: one point (or 2 extremely close points) - tangent case
            t = -B / 2.0
            xp = xs + smu0*t
            yp = ys
            zp = zs - mu0*t
        else:
            # 2 points, general case
            sqrtD = np.sqrt(D)
            t_hi_z = (-B - sqrtD) / 2.0 # this gives higher zp
            #t_lo_z = (-B + sqrtD) / 2.0
        
            z_hi =  zs - mu0*t_hi_z
            #z_lo =  zs - mu0*t_lo_z
            
            xp = xs + smu0*t_hi_z
            yp = ys
            zp = z_hi
        # if ... else
    
    return xp, yp, zp
#-EOF-