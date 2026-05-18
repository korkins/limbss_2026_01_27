# Summary
This repository contains a Python prototype of an RT code simulating single scattering of monochromatic sunlight at the limb. Prototype means the code is may not be thoroughly optimised numerically (e.g., it is in Python), twilight solar-view geometry is omitted, aerosol scenario is not yet included, etc. However, comments in the functions are thorough. Optical parameteres of shells, constituing the atmosphere, change along the radius, `r` (km). Extinction and scattering coefficients change lenearly with `r`, e.g.: ext = a + b * r ; a & b are constant within each shell (optical layer). The code was tested vs. pure Rayleigh and Rayleigh with absorpbing gas scenarious defined in Zawada D., et al., *Systematic comparison of vectorial spherical radiative transfer models in limb scattering geometry*, Atmos. Meas. Tech., 14, 3953–3972, https://doi.org/10.5194/amt-14-3953-2021, 2021.  

# Instructions

# Code Structure (LOC & Tree)
