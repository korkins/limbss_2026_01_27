import time
import numpy as np
from antideriv import antideriv

def tau_segment(xp, yp, zp, xq, yq, zq, exta, extb):
    """
    Compute optical depth (tau) along a straight segment from P to Q
    within a single spherical layer where extinction varies linearly
    with radius: ext(r) = exta*r + extb.

    Uses the exact antiderivative of r(s) ds via 'antideriv'.

    Parameters
    ----------
    xp, yp, zp : float
        Coordinates of the starting point P
    xq, yq, zq : float
        Coordinates of the end point Q
    exta, extb : float
        Linear extinction coefficients for the layer

    Returns
    -------
    tau : float
        Optical depth along the segment P->Q
    """

    # Vector from P to Q
    dx = xq - xp
    dy = yq - yp
    dz = zq - zp
    L = np.sqrt(dx*dx + dy*dy + dz*dz)

    if L == 0.0:
        return 0.0

    # Unit direction vector
    lx = dx / L
    ly = dy / L
    lz = dz / L

    # Radius at P
    rp = np.sqrt(xp*xp + yp*yp + zp*zp)

    # Cosine between radius vector and LOS
    mu = (xp*lx + yp*ly + zp*lz) / rp

    # Exact integral
    tau = (
        exta * (antideriv(L, rp, mu) - antideriv(0.0, rp, mu))
        + extb * L
    )

    return tau

if __name__ == "__main__":

    t0 = time.time()

    xp =   0.0000
    yp =   0.0000
    zp = 106.0000
    
    xq =   0.0000
    yq =   0.0000
    zq = 105.0000
    
    aext = 0.0006
    bext = 0.0600
    
    tau_sol = tau_segment(xp, yp, zp, xq, yq, zq, aext, bext)
    print(f"tau_sol = {tau_sol: 12.4e}")

    # --- Elapsed Time ---
    elapsed_time_sec = time.time() - t0
    print("total elapsed time:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time_sec)))
#-EOF-