import time
import numpy as np

def antideriv(s, rp, mu):
    """
    Computes integral{r(s) ds}, AKA antiderivative of r(s), in closed form;
    the s-independent constant is omitted as it will cancel out at calculation
    of the definite integral.
    
    Parameters
    ----------
    s : float
        Point of evaluation
    rp, mu : float
        Parameters

    Returns
    -------
    f : float
        Antiderivative of r(s); constants omitted
    
    Notes
    -----
    See [1: p.76] for Eq.380.201 or [2: p.320] for Eq.230.
    
    mu = +/- 1 is a special case:
        (1.0 - mu**2) == 0.0 but
        log(s + mu*rp + rs) = log(0) = inf
        
    This is handeled by takinf the limit mu --> +/- 1:
        integral{r(s)ds} = rp*s + (1/2)*mu*s**2 (note mu keeps sign here)
        
    Refs
    -----
    1. Dwight HB, Tables of Integrals and other Mathematical Data, 9th Ed, Macmillan, 1961
    2. Zwillinger D, CRC Standard Mathematical Tables and Formulas, 33d Ed, CRC Press, 2018 
        
    """
    
    tiny_mu = 1.0e-12
    
    if np.abs(abs(mu) - 1.0) < tiny_mu:
        antiderivative = rp*s + 0.5*mu*s**2 # note mu keeps sign
        
    else:
        rs = np.sqrt(rp**2 + 2.0*mu*rp*s + s**2)
        
        antiderivative = 0.5*( (s + mu*rp)*rs + (1.0 - mu**2)*rp**2*np.log(s + mu*rp + rs) )
    
    return antiderivative

if __name__ == "__main__":

    t0 = time.time()

    L  =   1.0
    mu =  -1.0    
    rp = 106.0
    
    antidir = antideriv(L, rp, mu)
    print(f"antiderivative = {antidir: 8.4f}")

    # --- Elapsed Time ---
    elapsed_time_sec = time.time() - t0
    print("total elapsed time:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time_sec)))
#-EOF-