import numpy as np

def refine_rgrid(r, npnts):
    """
    Construct refined grid of radii r with equidistant subdivision of layers
    and return layer-membership mask.

    The input grid r defines nlr = len(r) - 1 layers.
    Each layer i spans the interval [r[i], r[i+1]) - note half-open interval
    Each r[i+1] belongs to the layer above, i.e., i+1.
    The very top boundary r = r_max belongs to the very top layer by definition.

    Parameters
    ----------
    r : ndarray
        Monotonically increasing shell radii
    npnts : int
        Number of interior points, equidistant over r, inserted per layer

    Returns
    -------
    rr : ndarray
        Refined height grid including all boundaries and interior points
    mask : ndarray
        Integer layer index (0-based) for each element of rr
    """
    
    nr = len(r)
    nlr = nr - 1

    if npnts == 0:
        rr = r.copy()
        mask = np.empty(nr, dtype=np.int64)
        mask[0] = 0
        for i in range(1, nr):
            mask[i] = i - 1
        return rr, mask

    # total points:
    # original boundaries + npnts interior points per layer
    nrr = nr + nlr * npnts
    rr = np.empty(nrr)
    mask = np.empty(nrr, dtype=np.int64)

    ix = 0

    # BOA
    rr[ix] = r[0]
    mask[ix] = 0
    ix += 1

    for ilr in range(nlr):
        r0 = r[ilr]
        r1 = r[ilr + 1]
        dr = (r1 - r0) / (npnts + 1)

        # interior points
        for ir in range(1, npnts+1):
            rr[ix] = r0 + ir * dr
            mask[ix] = ilr
            ix += 1
        
        # top boundary belongs to the highest layer
        rr[ix] = r1
        mask[ix] = ilr + 1
        ix += 1
    
    # by definition: the very top shell, r=r_max, belongs to the very top layer
    mask[nrr - 1] = nlr - 1
    
    return rr, mask

#------------------------------------------------------------------------------
if __name__ == "__main__":
    
    Re   = 100.0
    hkm  = np.arange(7)
    npnts = 3

    r = Re + hkm
    rr, mask = refine_rgrid(r, npnts)

    print("Input radii (km):")
    print(r)

    print("\nRefined grid with layer mask:")
    print(" irr      rr[km]    layer")
    print("-------------------------")
    for irr in range(len(rr)):
        print(f"{irr:4d}  {rr[irr]:10.4f}  {mask[irr]:5d}")
    
#-EOF-