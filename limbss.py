import time
import numpy as np
from xyz_sol_shell import xyz_sol_shell
from tau_segment import tau_segment
from refine_rgrid import refine_rgrid
from scipy.integrate import simpson
from xy_los_shell import xy_los_shell
from getidx import getix
import matplotlib.pyplot as plt

def limbss(npnts_int, z_tp, mu0, caz, saz, r, aext, bext, asca, bsca, phasef):
    '''    
    Parameters
    ----------
    z_tp : TYPE
        DESCRIPTION.
    mu0 : TYPE
        DESCRIPTION.
    caz : TYPE
        DESCRIPTION.
    saz : TYPE
        DESCRIPTION.
    r : TYPE
        DESCRIPTION.
    aext : TYPE
        DESCRIPTION.
    bext : TYPE
        DESCRIPTION.
    asca : TYPE
        DESCRIPTION.
    bsca : TYPE
        DESCRIPTION.
    phasef : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    # number of shells & optical layers
    nr = len(r)
    nlr = nr - 1
    
    # extended
    rr, mask = refine_rgrid(r, npnts_int)
    nrr = len(rr)
    
    # number of shells above TP
    count = np.sum(rr > z_tp)
    
    # number of s-points for integration, including TP
    ns = 2*count + 1
    
    x_los = np.full(ns, -999.0)
    y_los = np.full_like(x_los, -999.0)
    z_los = np.full_like(x_los, -999.0)
    s_los = np.full_like(x_los, -999.0)
    r_los = np.full_like(x_los, -999.0)
    
    Jss = np.full_like(x_los, -999.0)

    mask_lr = np.full(ns, -999, dtype=np.int64)
    
    z_los[:] = z_tp # LOS is horizontal
    
    idx_cntr = ns // 2
    x_los[idx_cntr] = 0.0
    y_los[idx_cntr] = 0.0
    s_los[idx_cntr] = 0.0
    r_los[idx_cntr] = z_tp
    
    code, ir_below_tp, ir_above_tp = getix(r, z_tp)
    mask_lr[idx_cntr] = ir_above_tp - 1 # layer index is the "left" one
    
    idx_right = idx_cntr + 1
    idx_left = idx_cntr - 1
    
    code, irr_below_tp, irr_above_tp = getix(rr, z_tp)
    
    for ir in range(irr_above_tp, nrr):
        d = np.sqrt(rr[ir]**2 - z_tp**2)
        x = d * caz
        y = d * saz
        x_los[idx_right] =  x
        y_los[idx_right] =  y
        x_los[idx_left]  = -x
        y_los[idx_left]  = -y
        s_los[idx_right]  =  d
        s_los[idx_left]  = -d
        r_los[idx_right] = rr[ir]
        r_los[idx_left] = rr[ir]
        mask_lr[idx_right] = mask[ir]
        mask_lr[idx_left] = mask[ir]
        idx_right += 1
        idx_left  -= 1
        
    # loop over points 's' starting from the right (observation point)
    
    x_los_shell, y_los_shell = xy_los_shell(rr, z_tp, caz, saz)
    tau_los_plot = np.full(ns, -999.0)
    tau_sol_plot = np.full_like(tau_los_plot, -999.0)
    Jss_plot = np.full_like(tau_los_plot, -999.0)
    
    js = ns - 1
    ilr = mask_lr[js]
    tau_los = 0.0
    
    tau_sol = 0.0
    tau_los_plot[js] = tau_los
    tau_sol_plot[js] = tau_sol
    
    # Tsol = Tlos = 1.0
    Jss[js] = (asca[ilr]*r[nr-1] + bsca[ilr])*phasef[ilr]
    Jss_plot[js] = Jss[js]
    
    for js in range(ns - 2, -1, -1):
        ilr = mask_lr[js]
        xs = x_los[js]
        ys = y_los[js]
        xt = x_los[js + 1]
        yt = y_los[js + 1]
        tau_los = tau_los + tau_segment(xs, ys, z_tp, xt, yt, z_tp, aext[ilr], bext[ilr])
                
        # solar beam attenuation: always starts at the outermost shell ...
        tau_sol = 0.0
        ir = nr - 1
        xp, yp, zp = xyz_sol_shell(mu0, xs, ys, z_tp, r[ir])
        # ... and loop down through optical layers
        for jlr in range(nlr-1, ilr, -1): # from TOA to ilr+1, which must be included - hence ilr (without +1)
            ir = ir - 1
            xq, yq, zq = xyz_sol_shell(mu0, xs, ys, z_tp, r[ir])
    
            tau_sol = tau_sol + tau_segment(xp, yp, zp, xq, yq, zq, aext[jlr], bext[jlr]) # note jlr
            xp = xq
            yp = yq
            zp = zq

        # tiny segment along the solar beam between shell above TP and point S
        tau_sol = tau_sol + tau_segment(xp, yp, zp, xs, ys, z_tp, aext[ilr], bext[ilr])
        
        rs = r_los[js]
        Jss[js] = (asca[ilr]*rs + bsca[ilr])*phasef[ilr]*np.exp(-(tau_los + tau_sol))
 
        tau_los_plot[js] = tau_los
        tau_sol_plot[js] = tau_sol
        Jss_plot[js] = (asca[ilr]*rs + bsca[ilr])*phasef[ilr] # not T here
        
    Fo = 1.0
    Iss = simpson(Jss, s_los)*Fo/(4.0*np.pi)
        
    return Iss, s_los, tau_los_plot, tau_sol_plot, Jss_plot/(4.0*np.pi)
        
#------------------------------------------------------------------------------

if __name__ == "__main__":

    t0 = time.time()
    
    npnts_int = 2

    Re   = 100.0
    hkm  = np.arange(7) # to test vs. drawing
    h_tp = 2.5
    sza_tp = 0.0
    raz_tp = 0.0
    
    deg_to_rad = np.pi / 180.0
    nlr = len(hkm) - 1

    # (fake) linear approximation: ext(r) = aext*r + bext = sca(r)
    aext = np.zeros(nlr)
    bext = np.zeros(nlr)
    for ilr in range(nlr):
        aext[ilr] = 1.0e-4 * (ilr + 1)
        bext[ilr] = 1.0e-2 * (ilr + 1)
    asca = aext
    bsca = bext
      
    
    mu0 = np.cos(sza_tp * deg_to_rad)
    cos_raz = np.cos(raz_tp * deg_to_rad)
    sin_raz = np.sin(raz_tp * deg_to_rad)
    sin_sza = np.sqrt(1.0 - mu0**2)
    cos_scat_angle = sin_sza * cos_raz
    scat_angle_deg = np.arccos(cos_scat_angle) / deg_to_rad
    phase_function = 0.75 * (1.0 + cos_scat_angle**2)
    phasef = np.zeros(nlr)
    phasef[:] = phase_function
        
    
    Iss, s_los, tau_los_plot, tau_sol_plot, Jss_plot = \
        limbss(npnts_int, Re + h_tp, mu0, cos_raz, sin_raz, Re + hkm, aext, bext, asca, bsca, phasef) 
        
    print(f"Iss = {Iss: 8.4e}")
    
    # --- Print results ---
    print("Heights [km]:")
    for i, h in enumerate(hkm):
        print(f"{i:3d}  {h:6.2f}")

    print("\nTangent height [km]:", h_tp)
    
    fig, ax = plt.subplots(3, 1, figsize=(8, 6), sharex=True)

    ax[0].plot(s_los, Jss_plot, color='k', linewidth=2)
    ax[0].set_ylabel('Jss', fontsize=12)
    ax[0].grid(True)
    
    ax[1].plot(s_los, tau_sol_plot, color='red', linewidth=2)
    ax[1].set_ylabel('tau_sol', fontsize=12)
    ax[1].grid(True)

    ax[2].plot(s_los, tau_los_plot, color='blue', linewidth=2)
    ax[2].set_ylabel('tau_los', fontsize=12)
    ax[2].set_xlabel('s, km', fontsize=12)
    ax[2].grid(True)

    plt.tight_layout()
    plt.show()

    # --- Elapsed Time ---
    elapsed_time_sec = time.time() - t0
    print("total elapsed time:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time_sec)))
#-EOF-