# Summary
This repository contains a Python prototype of an RT code simulating single scattering of monochromatic sunlight at the limb. Prototype means the code is may not be thoroughly optimised numerically (e.g., it is in Python), twilight solar-view geometry is omitted, aerosol scenario is not yet included, etc. However, comments in the functions are thorough.  
Optical parameteres of shells, constituing the atmosphere, change along the radius, `r` (km). Extinction and scattering coefficients change lenearly with `r`, e.g.: ext = a + b * r; a & b are constant within each shell (optical layer).  
# Instructions

# Code Structure (LOC & Tree)
