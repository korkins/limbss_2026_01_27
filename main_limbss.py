import time
import numpy as np
import limbss
import linear_coefficients as lincoef
from scipy.integrate import simpson
import matplotlib.pyplot as plt

if __name__ == "__main__":
    t0 = time.time()
    # Zawada et al, AMT (2021)
    
    # number of extra points between shells (internal) for integration along LOS
    npnts_int = 1


    imodel = 0          # SASKTRAN - do not change
    itest_case = 0      # polarized single scattering no refraction - do not change
    isolar = 0
    icomposition = 0    # 0 or 1
    ialbedo = 0         # black - do not change
    iwavelength = 2
    istokes = 0         # Intensity - do not change
    
    itp_plot = 30
    
    
    Re = 6371.0
    cm_per_km = 1.0e+5
    path_bmrk = "./ort_limb_benchmarks/"
    
    if icomposition == 0:
        abs_on_off = 0.0
    else:
        abs_on_off = 1.0
    
    fname = "air_ozone_numden_profiles.txt"
    data = np.loadtxt(path_bmrk + fname, skiprows=1)
    hkm = data[:, 0] / 1000.0
    air_numden = data[:, 1]
    ozone_numden = data[:, 2]
    
    fname = "tangent_heights.txt"
    hkm_tp = np.loadtxt(path_bmrk + fname, skiprows=1)
    ntp = len(hkm_tp)
    
    fname = "xsection_rayleigh_ozone.txt"
    data = np.loadtxt(path_bmrk + fname, skiprows=1)
    wav_nm_all = data[:, 0]
    rayleigh_xsection = data[:, 1] 
    ozone_xsection = data[:, 2]
    nwav = len(wav_nm_all)

    print("Rayleigh OT, ozone OT, SSA:")
    for iwav in range(nwav):
        sca_coef = rayleigh_xsection[iwav] * air_numden * cm_per_km
        abs_coef = ozone_xsection[iwav] * ozone_numden * cm_per_km
        tau_ray = simpson(sca_coef, hkm)
        tau_o3 = simpson(abs_coef, hkm)
        ssa = tau_ray / (tau_ray + tau_o3)
        print(f"{wav_nm_all[iwav]: 8.1f} nm: {tau_ray:.4f}  {tau_o3:.4f}  {ssa: .4f}  ")
        
    fname = "solar_zenith_azimuth.txt"
    data = np.loadtxt(path_bmrk + fname, skiprows=1)
    sza_tp_all = data[:, 0]
    raz_tp_all = data[:, 1]


    wav_nm = wav_nm_all[iwavelength]
    sza_tp = sza_tp_all[isolar]
    raz_tp = raz_tp_all[isolar]


    print("Echo input:")
    print(f"  sza ={sza_tp: .1f}  saa ={raz_tp: .1f}")
    print(f"  wavelength ={wav_nm: .1f}")
    print(f"  rayleigh x-section = {rayleigh_xsection[iwavelength]: .8e}")
    print(f"  ozone x-section = {ozone_xsection[iwavelength]: .8e}")
    print(f"  abs_on_off = {abs_on_off: .1f}")
    print(f"  number of TPs: {ntp: 3d}")
     
    
    # units: cm2 * 1/cm3 * cm/km = 1/km for both
    sca_coef = rayleigh_xsection[iwavelength] * air_numden * cm_per_km
    abs_coef = ozone_xsection[iwavelength] * ozone_numden * cm_per_km
    ext_coef = sca_coef + abs_coef * abs_on_off
    
    nlr = len(hkm) - 1
    aext, bext, asca, bsca = \
        lincoef.linear_coefficients(Re + hkm, ext_coef, sca_coef)


    deg_to_rad = np.pi/180.0
    mu0 = np.cos(sza_tp * deg_to_rad)
    cos_raz = np.cos(raz_tp * deg_to_rad)
    sin_raz = np.sin(raz_tp * deg_to_rad)
    sin_sza = np.sqrt(1.0 - mu0**2)
    cos_scat_angle = sin_sza * cos_raz
    scat_angle_deg = np.arccos(cos_scat_angle) / deg_to_rad
    phase_function = 0.75 * (1.0 + cos_scat_angle**2)
    phasef = np.zeros(nlr)
    phasef[:] = phase_function
    
    Iss = np.full(ntp, -999.0)
    for itp in range(ntp):
        Iss[itp], s_los, tau_los, tau_sol, Jss = \
            limbss.limbss(npnts_int, Re + hkm_tp[itp], mu0, cos_raz, sin_raz, Re + hkm, aext, bext, asca, bsca, phasef)
        if itp == itp_plot:
            s_los_plot = s_los
            tau_los_plot = tau_los
            tau_sol_plot = tau_sol
            Jss_plot = Jss
    
    elapsed_time_sec = time.time() - t0
    
    fname_Iss = (
        f"Iss_imod{imodel}"
        f"_itst{itest_case}"
        f"_isol{isolar}"
        f"_icom{icomposition}"
        f"_ialb{ialbedo}"
        f"_iwav{iwavelength}"
        f"_istk{istokes}.txt"
        )
    
    out = np.zeros((ntp, 2))
    out[:, 0] = hkm_tp
    out[:, 1] = Iss
    np.savetxt(fname_Iss, out, fmt=['%4.1f', '%12.4e'], header='hkm_tp   Iss')
 
    if (True):
        # --- Print results ---
 
        fig, ax = plt.subplots(4, 1, figsize=(8, 12), sharex=True)

        ax[0].plot(s_los_plot, Jss_plot, color='k', linewidth=2)
        ax[0].set_ylabel('Jss', fontsize=12)
        ax[0].grid(True)
 
        T_sol = np.exp(-tau_sol_plot)
        ax[1].plot(s_los_plot, T_sol, color='red', linewidth=2)
        ax[1].set_ylabel('T_sol', fontsize=12)
        ax[1].grid(True)

        J_plot = Jss_plot * T_sol 
        ax[2].plot(s_los_plot, J_plot, color='magenta', linewidth=2)
        ax[2].set_ylabel('J', fontsize=12)
        ax[2].grid(True)

        T_los = np.exp(-tau_los_plot)
        ax[3].plot(s_los_plot, T_los, color='blue', linewidth=2)
        ax[3].set_ylabel('T_los', fontsize=12)
        ax[3].set_xlabel('s, km', fontsize=12)
        ax[3].grid(True)

        plt.tight_layout()
        plt.show()
        
        value = simpson(J_plot*T_los, s_los_plot)
        print(f"Iss = {value:.6e} at hkm_tp = {hkm_tp[itp_plot]: .1f} km")
        
        half = len(s_los_plot)//2
        value = simpson(J_plot[half:]*T_los[half:], s_los_plot[half:])
        print(f"Iss[half:] = {value:.6e} at hkm_tp = {hkm_tp[itp_plot]: .1f} km")
 
    # --- Elapsed Time ---
    print(f"total elapsed time: {elapsed_time_sec: .3f} s")