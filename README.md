# Summary
This repository contains a Python prototype of an RT code simulating single scattering of monochromatic sunlight at the limb. Prototype means the code is may not be thoroughly optimised numerically (e.g., it is in Python), twilight solar-view geometry is omitted, aerosol scenario is not yet included, comments in the functions are not yet thorough or could be inaccuarte. Optical parameteres of shells, constituing the atmosphere, change along the radius, `r` (km). Extinction and scattering coefficients change along `r`lenearly: ext = a + b * r ; a & b are constant within each shell (optical layer). The code was tested vs. all pure and absorbing Rayleigh scenarious from Zawada D., et al., *Systematic comparison of vectorial spherical radiative transfer models in limb scattering geometry*, Atmos. Meas. Tech., 14, 3953–3972, https://doi.org/10.5194/amt-14-3953-2021, 2021. Note, the paper also ignores the twilight condition. 

# Instructions
1. The GitHub repository contains 9 py-sources and 2 folders. The `tests` folder is discussed below, at step #6. The `ort_limb_benchmarks` folder contains input for the benchmarks (the word 'input' should have been added to the folder name...). This data was extracted from Zawada’s 'NetCDF' file and placed into 'txt' files for a) convenience – one does not have to read their NetCDF, and b) better data visibility – i.e., the sequence of the solar geometry in their NetCDF differs from that published in the paper.
2. In `main_limbss.py`, line 29, `path_bmrk = "./ort_limb_benchmarks/"` – this relative path should work as is, but if not, update as necessary.
3. Also, in the `main_limbss.py`:
-Line 13: `npnts_int = 1` – number of internal (auxiliary) shells in between the main shells; `npnts_int = 0, 1, 2` were tried;  
-Lines 16-22: 0-offset indices. Only `isolar`, `iwavelength`, and `icomposition` have effect to select different solar geometry, wavelength, and pure Rayleigh or Rayleigh + ozone case. Other indices are used only to create a unique output file name (see below) and quickly locate corresponding data in the NetCDF benchmark file.
4. Line 24: `itp_plot = 30` - 0-offset index corresponding to tangent height at which solution will be printed on screen (1 point; all points go to a txt file when finished). Due to Zawada’s input, its value is close to the actual tangent height, i.e. index 30 corresponds to 30.5 km
5. By default, one should get the following on the screen (almost immediately):
```
Rayleigh OT, ozone OT, SSA:  
300.0 nm: 1.2002  2.9165   0.2915
315.0 nm: 0.9747  0.3392   0.7418
351.0 nm: 0.6167  0.0010   0.9984
435.0 nm: 0.2522  0.0007   0.9973
442.0 nm: 0.2361  0.0014   0.9942
525.0 nm: 0.1165  0.0176   0.8690
600.0 nm: 0.0676  0.0419   0.6176
675.0 nm: 0.0419  0.0121   0.7759
943.0 nm: 0.0109  0.0000   1.0000
1020.0 nm: 0.0079  0.0000   1.0000
1700.0 nm: 0.0010  0.0000   1.0000

Echo input:  
sza = 70.0  saa = 30.0
wavelength = 351.0
rayleigh x-section =  2.87884699e-26
ozone x-section =  1.22579807e-22
abs_on_off =  0.0
number of TPs:  80
```

# Code Structure (LOC & Tree)
