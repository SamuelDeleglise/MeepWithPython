from pylab import *
from numpy import *
from meep import *

import time

res = 10
gridSizeX = 32.0
gridSizeY = 20.0
border = 1.0
srcFreqCenter = 1.0 #gaussian source center frequency
srcPulseWidth = 0.5 #gaussian source pulse width
v = vol2d(2, 2, 10)

s = structure(v, EPS, pml(0.1,Y) )
f = fields(s)

class epsilon(PolygonCallback2D):
    def __init__(self):
        PolygonCallback2D.__init__(self)
        master_printf("Callback function for epsilon contructed. Now defining the polygons...\n")
        #Points outside this polygon are assumed to have epsilon 1.0
        pol_points = numpy.zeros((5,2))
        pol_points[0] = [0.0, gridSizeY/2]
        pol_points[1] = [gridSizeX, gridSizeY/2]
        pol_points[2] = [gridSizeX, gridSizeY/2 + 2.0]
        pol_points[3] = [0.0, gridSizeY/2 + 2.0]
        pol_points[4] = [0.0, gridSizeY/2]
        
        
        self.add_polygon(pol_points, 1.0)
        master_printf("Polygons OK.\n")




class amplitudeFactor(Callback):
    def __init__(self):
        Callback.__init__(self)
        master_printf("Callback function for amplitude factor activated.\n")

    def complex_vec(self,vec):
        #BEWARE, these are coordinates RELATIVE to the source center !!!!
        x = vec.x()
        y = vec.y()
        master_printf("Fetching amplitude factor for x=%f - y=%f\n" %(x,y) )
        result = complex(exp(-(x/2)**2),0)
        return result

#we create a structure with PML of thickness = 1.0 on all boundaries,
#in all directions,
#using the material function EPS
my_space = voltwo(gridSizeX,gridSizeY,res)

material = epsilon()
set_EPS_Callback(material.__disown__())
#use_averaging(False)  #--> uncomment this line to disable eps-averaging!
my_structure = structure(my_space, EPS, pml(1.0))
my_fields = fields(my_structure)
#define a gaussian line source of length 'wgWidth' at X=wgLength/2, Y=padSize
srcGaussian = gaussian_src_time(srcFreqCenter, srcPulseWidth )
srcGeo = volume(vec(border,border),vec(gridSizeX-border,border))
af = amplitudeFactor()
set_AMPL_Callback(af.__disown__())
my_fields.add_volume_source(Ez, srcGaussian, srcGeo, AMPL)
master_printf("Field created...\n")

#my_file = prepareHDF5File("./coucou.h5")
#runUntilFieldsDecayed(my_fields, my_space, Ez, vec(2.0,2.0),pH5OutputIntervalSteps=1000,pHDF5OutputFile = my_file)


##===============================================================
##            CREATE FLUX PLANE
##===============================================================
fluxplanes = []
my_fluxes = []
for i in arange(0.1,15,0.5):
    fluxplanes.append(volume(vec(gridSizeX/2-i,gridSizeY-border),vec(gridSizeX+i,gridSizeY - border)))
    my_fluxes.append(my_fields.add_dft_flux_plane(fluxplanes[-1],srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 100))


##===============================================================
##            CREATE SAVE DIRECTORY
##===============================================================
savedir = str(__file__[:-3])
my_fields.outdir = savedir
import os
if os.path.exists(savedir):
    yn = raw_input("savedir allready exists, erase all in it?")
    if yn == "y":
        for f in os.listdir(savedir):
            os.remove(savedir + "/" + f)
else:
    os.mkdir(savedir)

##===============================================================
##            STEP THROUGH THE FIELD
##===============================================================
my_file = prepareHDF5File(savedir+"/images.h5")
steps_total = 800
n_frames = 20
steps_per_frame = steps_total/n_frames
for i in range(steps_total):
    my_fields.step()
    if i%steps_per_frame == 0 :
        my_fields.output_hdf5(Ez, my_space.surroundings(),my_file,1)

del my_file
##===============================================================
##            CONVERT H5 FILES TO PNG
##   see http://ab-initio.mit.edu/wiki/index.php/Meep_Tutorial
##===============================================================
import subprocess,glob
#for f in glob.glob(savedir +"/"+"*.h5"):
import os
curdir = os.getcwd()
os.chdir(savedir)

try:
    master_printf("building individual png...\n")
    t = time.time()
    sp = subprocess.Popen(["h5topng","-t 0:%i"%(n_frames-1),"-R","-Zc","dkbluered","images.h5"])
    sp.wait()
    master_printf("done in " + str(time.time()-t) + "\n")
    master_printf("making gif...\n")
    t = time.time()
    sp = subprocess.Popen(["convert","*.png", "animation.gif"])
    sp.wait()
    master_printf("done in " + str(time.time()-t) + "\n")
    master_printf("done\n")
except e:
    raise e
finally:
    os.chdir(curdir)

##===============================================================
##            PLOT FLUXES
##===============================================================
fig = figure()
center_freq = []
for index,i in enumerate(my_fluxes):
    fluxdata = getFluxData(i)
    center_freq.append(fluxdata[1][35])
plot(diff(center_freq),label = "low_freq")
side_freq = []
for index,i in enumerate(my_fluxes):
    fluxdata = getFluxData(i)
    side_freq.append(fluxdata[1][65])
plot(diff(side_freq),label = "high_freq")
legend()
show()
fig.savefig(savedir  +"/profiles.png")
