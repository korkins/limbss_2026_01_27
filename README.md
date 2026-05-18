# Summary
This repository contains a Python prototype of an RT code simulating single scattering of monochromatic sunlight at the limb. Prototype means the code is may not be thoroughly optimised numerically (e.g., it is in Python), twilight solar-view geometry is omitted, aerosol scenario is not yet included, comments in the functions are not yet thorough or could be inaccuarte. Optical parameteres of shells, constituing the atmosphere, change along the radius, `r` (km). Extinction and scattering coefficients change along `r`lenearly: ext = a + b * r ; a & b are constant within each shell (optical layer). The code was tested vs. all pure and absorbing Rayleigh scenarious from Zawada D., et al., *Systematic comparison of vectorial spherical radiative transfer models in limb scattering geometry*, Atmos. Meas. Tech., 14, 3953–3972, https://doi.org/10.5194/amt-14-3953-2021, 2021.  

# Instructions

# Code Structure (LOC & Tree)
