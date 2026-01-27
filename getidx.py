import numpy as np

def getix(x, x0):
    """
    Returns indices ixlo, ixhi of an interval x[ixlo]:x[ixhi] containing x0.
    Assumptions: 
        1) 0-offset for result
        2) xmin < x0 < xmax, i.e. both ixlo & ixhi exist (note strict inequality!)
        3) elements of x are sorted in ascending order
        4) elements of x are equidistant (regular grid)

    Parameters
    ----------
    x : ndarray
        Array of elements in ascending order
    x0 : float
        Value whose location is sought. 

    Returns
    -------
    code : integer
        -1 if no answer, 0 otherwise
    ixlo, ixhi : integers
        Desired indices
    """
    code = -1
    ixlo = None
    ixhi = None
    
    nx = len(x)
    xmin = x[0]
    xmax = x[nx - 1]
    dx = (xmax - xmin) / (nx - 1)
    if (xmin < x0 < xmax):
        ixlo = int((x0 - xmin) // dx) # floor division
        ixlo = min(ixlo, nx - 2)      # 2nd check in addition to xmin < x0 < xmax
        ixhi = ixlo + 1
        code = 0

    return code, ixlo, ixhi

# -------------------------
# Test example
# -------------------------
if __name__ == "__main__":

    nx = 6
    dx = 2.5
    xmin = 100.0   # grid starts at 100

    x = np.zeros(nx)
    for i in range(nx):
        x[i] = xmin + i * dx
    print("x grid :", x)
        
    print("\nTest 1: general case")

    x0 = 106.3   # between 105.0 and 107.5

    code, ixlo, ixhi = getix(x, x0)

    print("x0     :", x0)
    print("code   :", code)
    print("ixlo   :", ixlo)
    print("ixhi   :", ixhi)

    if code == 0:
        print("x[ixlo] =", x[ixlo])
        print("x[ixhi] =", x[ixhi])
        print("Check  :", x[ixlo] <= x0 <= x[ixhi])
    
    print("\nTest 2: x0 is at some node (returns x0 = x[ixlo]")

    x0 = 105.0

    code, ixlo, ixhi = getix(x, x0)

    print("x0     :", x0)
    print("code   :", code)
    print("ixlo   :", ixlo)
    print("ixhi   :", ixhi)

    if code == 0:
        print("x[ixlo] =", x[ixlo])
        print("x[ixhi] =", x[ixhi])
        print("Check  :", x[ixlo] <= x0 <= x[ixhi])
    
    print("\nTest 3: x0 = xmin (returns nones)")

    x0 = np.amin(x)

    code, ixlo, ixhi = getix(x, x0)

    print("x0     :", x0)
    print("code   :", code)
    print("ixlo   :", ixlo)
    print("ixhi   :", ixhi)

    if code == 0:
        print("x[ixlo] =", x[ixlo])
        print("x[ixhi] =", x[ixhi])
        print("Check  :", x[ixlo] <= x0 <= x[ixhi])
    
    print("\nTest 4: x0 = xmax (returns nones)")

    x0 = np.amax(x)

    code, ixlo, ixhi = getix(x, x0)

    print("x0     :", x0)
    print("code   :", code)
    print("ixlo   :", ixlo)
    print("ixhi   :", ixhi)

    if code == 0:
        print("x[ixlo] =", x[ixlo])
        print("x[ixhi] =", x[ixhi])
        print("Check  :", x[ixlo] <= x0 <= x[ixhi])