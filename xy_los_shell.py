import time
import numpy as np
from getidx import getix

def xy_los_shell(r, z_tp, caz, saz):
    
    #diagnostic_flag = False
    
    nr = len(r)
    
    #ir_above_tp = np.searchsorted(r, z_tp, side="right") # first > value
    code, ir_below_tp, ir_above_tp = getix(r, z_tp)
     
    # number of shells above TP
    count = np.sum(r > z_tp)
    
    # number of LOS-shell intersections, both sides - always even
    num_los_shell = 2 * count
    
    x_los_r = np.zeros(num_los_shell)
    y_los_r = np.zeros_like(x_los_r)
    
    idx_right = count
    idx_left = idx_right - 1
    
    for ir in range(ir_above_tp, nr):
        d = np.sqrt(r[ir]**2 - z_tp**2)
        x = d * caz
        y = d * saz
        x_los_r[idx_right] =  x
        y_los_r[idx_right] =  y
        x_los_r[idx_left]  = -x
        y_los_r[idx_left]  = -y
        idx_right += 1
        idx_left  -= 1
        
    return x_los_r, y_los_r
        
#------------------------------------------------------------------------------
if __name__ == "__main__":

    t0 = time.time()

    Re   = 100.0
    hkm  = np.arange(7)
    h_tp = 2.5
    raz_deg = 45.0
    
    cos_raz =  np.cos(raz_deg * np.pi / 180.0)
    sin_raz =  np.sin(raz_deg * np.pi / 180.0)

    x, y = xy_los_shell(Re + hkm, Re + h_tp, cos_raz, sin_raz)
    
    npnts = len(x)

    # --- Print results ---
    print("Heights [km]:")
    for i, h in enumerate(hkm):
        print(f"{i:3d}  {h:6.2f}")

    print("\nTangent height [km]:", h_tp)

    print("\n(x, y) for LOS-shell intersection:")
    for ipnt in range(npnts):
        print(f"point {ipnt:2d} : x = {x[ipnt]:12.6e}, y = {y[ipnt]:12.6e}")

    # --- Elapsed Time ---
    elapsed_time_sec = time.time() - t0
    print("total elapsed time:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time_sec)))
#-EOF-

    