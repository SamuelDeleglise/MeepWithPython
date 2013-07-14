##===============================================================
##       IMPORTS
##===============================================================
from pylab import *
from numpy import *
from meep import *
import time
from my_meep.structures.briques import make_fields,my_space,srcFreqCenter,srcPulseWidth,res,gridSizeX,gridSizeY,border
##===============================================================
##       NUMERICAL PARAMETERS
##===============================================================
#srcFreqCenter = 0.25 #gaussian source center frequency
#srcPulseWidth = 1.5 #gaussian source pulse width

def dispersion(objet,n_k = 25,f = 1.0,df = 1.5,**kwds):
    """calculates and displays the bloch bands for the object with n_k values for k
    f is the central frequency, and df is the width
    extra kwds are passed to the constructor of objet()"""
    


import h5py
f = h5py.File("modesBlochVert.h5","w")
for sym,col in zip([-1,+1],["r","b"]):
    results = []
    ks = arange(0.0,0.51,0.01)
    for k in ks:
        my_fields = make_fields(bloch_k = k, symmetry = sym)
        results.append(runWithHarminv(my_fields,my_space,Ez,vec(gridSizeX/2.12345,gridSizeY/2.12345),srcFreqCenter,1.0,15))#pHDF5OutputFile = my_file,pHDF5OutputFilePhase3 = other_file))
    
    #b = array((len(results),)
    k_s = []
    w_s = []
    q_s = []
    for i,k in zip(results,ks):
        for j in i:
            k_s.append(k)
            w_s.append(float(real(j[0])))
            q_s.append(j[2])
    scatter(k_s,w_s,s = sqrt(array(q_s))/3,c = col,label = "sym = " + str(sym))
    
    plot([0,0.5,0.2],[0,0.5,0.8],"black")
    a = gca()
    a.set_xbound(0,0.5)
    a.set_ybound(0,0.8)
    legend(loc = "best")
    
    
    f["k_s_sym_" + str(sym)] = k_s
    f["w_s_sym" + str(sym)] = w_s
    f["q_s_sym" + str(sym)] = q_s
    show()
f.close()