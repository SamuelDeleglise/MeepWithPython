##===============================================================
##       IMPORTS
##===============================================================
from pylab import *
from numpy import *
from meep import *
import time
from my_meep.structures.layer import make_fields,my_space,gridSizeX,gridSizeY,border,srcFreqCenter,srcPulseWidth
from utils_meep import run_and_gif

    
my_fields = make_fields(bloch_k = 0)
my_fields_cal = make_fields(bloch_k = 0,calibre = True)

inp = []
out = []
for fields,ext in zip((my_fields_cal,my_fields),("cal","run")):
    input_plane = volume(vec(0,2*border),vec(gridSizeX,2*border))
    output_plane = volume(vec(0,gridSizeY-border),vec(gridSizeX,gridSizeY-border))
    input_flux = fields.add_dft_flux_plane(input_plane,srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 1000)
    output_flux = fields.add_dft_flux_plane(output_plane,srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 1000)

    run_and_gif(fields,my_space,steps = 10000,frames = 200,filename = "layerMirror_" + ext + ".gif")
    
        
    inp.append(array(getFluxData(input_flux)[1]))
    out.append(array(getFluxData(output_flux)[1]))
    freq = array(getFluxData(output_flux)[0])

plot(freq,out[0]/inp[0],label = "transmission no structure")
plot(freq,out[1]/inp[0],label = "transmission with structure")
plot(freq,(inp[0] - inp[1])/inp[0],label = "reflexion with structure")
    
legend(loc = "best")
a = gca()
a.set_xbound(-0.1,1.1)

show()