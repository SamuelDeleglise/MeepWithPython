'''
Conversion into Python of the 90 degree BENT WAVEGUIDE SAMPLE from Meep-Scheme tutorial
(see http://ab-initio.mit.edu/wiki/index.php/Meep_Tutorial )

Using the technique of POLYGONS to define the simulation geometry.

@author: Emmanuel.Lambert@intec.ugent.be
'''

from meep import *  # make it 'meep_mpi' for MPI-meep and 'meep' for non-MPI meep

from math import *
import numpy
import matplotlib.pyplot
import sys

res = 10.0
gridSizeX = 16.0
gridSizeY = 32.0
padSize = 4.0
wgLengthX = gridSizeX - padSize
wgLengthY = gridSizeY - padSize
wgWidth = 1.0 #width of the waveguide
wgHorYCen = padSize +  wgWidth/2.0 #horizontal waveguide center Y-pos
wgVerXCen = wgLengthX - wgWidth/2.0 #vertical waveguide center X-pos (in case there is a bend)
srcFreqCenter = 0.15 #gaussian source center frequency
srcPulseWidth = 0.1 #gaussian source pulse width
srcComp = Ez #gaussian source component


#this function plots the waveguide material as a function of a vector(X,Y)
class epsilon(PolygonCallback2D):
    def __init__(self, pIsWgBent):
        PolygonCallback2D.__init__(self)
        master_printf("Callback function for epsilon contructed. Now defining the polygons...\n")
        grid_step = 1.0 / float(res)
        if pIsWgBent:
           master_printf("Starting polygons for case without bend.\n")
           #Points outside this polygon are assumed to have epsilon 1.0
           pol_points = numpy.zeros((7,2))
           pol_points[0] = [0.0,padSize]
           pol_points[1] = [gridSizeX - padSize, padSize]
           pol_points[2] = [gridSizeX - padSize, gridSizeY]
           pol_points[3] = [gridSizeX - padSize - wgWidth, gridSizeY]
           pol_points[4] = [gridSizeX - padSize - wgWidth, padSize + wgWidth]
           pol_points[5] = [0.0, padSize + wgWidth]
           pol_points[6] = [0.0,padSize]
           self.add_polygon(pol_points, 12.0)
        else:
           master_printf("Starting polygons for case without bend.\n")
           #Points outside this polygon are assumed to have epsilon 1.0
           pol_points = numpy.zeros((5,2))
           pol_points[0] = [0.0, padSize]
           pol_points[1] = [gridSizeX, padSize]
           pol_points[2] = [gridSizeX, padSize+wgWidth]
           pol_points[3] = [0.0, padSize+wgWidth]
           pol_points[4] = [0.0, padSize]
           self.add_polygon(pol_points, 12.0)

           master_printf("Polygons OK for case without bend.\n")

class amplitudeFactor(Callback):
    def __init__(self):
        Callback.__init__(self)
        master_printf("Callback function for amplitude factor activated.\n")

    def complex_vec(self,vec):
        #BEWARE, these are coordinates RELATIVE to the source center !!!!
        x = vec.x()
        y = vec.y()
        master_printf("Fetching amplitude factor for x=%f - y=%f\n" %(x,y) )
        result = complex(1.0,0)
        return result

def createField(pCompVol, pWgLengthX, pWgWidth, pIsWgBent, pSrcFreqCenter, pSrcPulseWidth, pSrcComp):
        #we create a structure with PML of thickness = 1.0 on all boundaries,
        #in all directions,
        #using the material function EPS
        material = epsilon(pIsWgBent)
        set_EPS_Callback(material.__disown__())
        #use_averaging(False)  #--> uncomment this line to disable eps-averaging!
        s = structure(pCompVol, EPS, pml(1.0))
        f = fields(s)
        #define a gaussian line source of length 'wgWidth' at X=wgLength/2, Y=padSize
        srcGaussian = gaussian_src_time(pSrcFreqCenter, pSrcPulseWidth )
        srcGeo = volume(vec(1.0,padSize),vec(1.0,padSize+pWgWidth))
        af = amplitudeFactor()
        set_AMPL_Callback(af.__disown__())
        f.add_volume_source(pSrcComp, srcGaussian, srcGeo, AMPL)
        master_printf("Field created...\n")
        return f


master_printf("** Bent waveguide sample using POLYGONS to geometry, version 28-02-2011 **\n")

master_printf("Running on %d processor(s)...\n",count_processors())

#FIRST WE WORK OUT THE CASE WITH NO BEND
#----------------------------------------------------------------
master_printf("*1* Starting the case with no bend...\n")
#create the computational grid
noBendVol = voltwo(gridSizeX,gridSizeY,res)

#create the field
wgBent = 0 #no bend
noBendField = createField(noBendVol, wgLengthX, wgWidth, wgBent, srcFreqCenter, srcPulseWidth, srcComp)

bendFnEps = "./bentwgNB_Eps.h5"
bendFnEz = "./bentwgNB_Ez.h5"
#export the dielectric structure (so that we can visually verify the waveguide structure)
noBendDielectricFile =  prepareHDF5File(bendFnEps)
noBendField.output_hdf5(Dielectric, noBendVol.surroundings(), noBendDielectricFile)

#create the file for the field components
noBendFileOutputEz = prepareHDF5File(bendFnEz)

#define the flux plane for the reflected flux
noBendReflectedfluxPlaneXPos = 1.5 #the X-coordinate of our reflection flux plane
noBendReflectedFluxPlane = volume(vec(noBendReflectedfluxPlaneXPos,wgHorYCen-wgWidth),vec(noBendReflectedfluxPlaneXPos,wgHorYCen+wgWidth))
noBendReflectedFlux = noBendField.add_dft_flux_plane(noBendReflectedFluxPlane, srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 100)

#define the flux plane for the transmitted flux
noBendTransmfluxPlaneXPos = gridSizeX - 1.5;  #the X-coordinate of our transmission flux plane
noBendTransmFluxPlane = volume(vec(noBendTransmfluxPlaneXPos,wgHorYCen-wgWidth),vec(noBendTransmfluxPlaneXPos,wgHorYCen+wgWidth))
noBendTransmFlux = noBendField.add_dft_flux_plane(noBendTransmFluxPlane, srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 100 )

master_printf("Calculating...")
noBendProbingpoint = vec(noBendTransmfluxPlaneXPos,wgHorYCen) #the point at the end of the waveguide that we want to probe to check if source has decayed
runUntilFieldsDecayed(noBendField, noBendVol, srcComp, noBendProbingpoint, pHDF5OutputFile = noBendFileOutputEz, pH5OutputIntervalSteps=50)
master_printf("Done..!")

#construct 2-dimensional array with the flux plane data
#see python_meep.py
noBendReflFlux = getFluxData(noBendReflectedFlux)
noBendTransmFlux = getFluxData(noBendTransmFlux)

#save the reflection flux from the "no bend" case as minus flux in the temporary file 'minusflux.h5'
noBendReflectedFlux.scale_dfts(-1);
f = open("minusflux.h5", 'w') #truncate file if already exists
f.close()
noBendReflectedFlux.save_hdf5(noBendField, "minusflux")

del_EPS_Callback()


#AND SECONDLY FOR THE CASE WITH BEND
#----------------------------------------------------------------
master_printf("*2* Starting the case with bend...\n")
#create the computational grid
bendVol = voltwo(gridSizeX,gridSizeY,res)

#create the field
wgBent = 1 #there is a bend
bendField = createField(bendVol, wgLengthX, wgWidth, wgBent, srcFreqCenter, srcPulseWidth, srcComp)

#export the dielectric structure (so that we can visually verify the waveguide structure)
bendFnEps = "./bentwgB_Eps.h5"
bendFnEz = "./bentwgB_Ez.h5"
bendDielectricFile = prepareHDF5File(bendFnEps)
bendField.output_hdf5(Dielectric, bendVol.surroundings(), bendDielectricFile)

#create the file for the field components
bendFileOutputEz = prepareHDF5File(bendFnEz)

#define the flux plane for the reflected flux
bendReflectedfluxPlaneXPos = 1.5 #the X-coordinate of our reflection flux plane
bendReflectedFluxPlane = volume(vec(bendReflectedfluxPlaneXPos,wgHorYCen-wgWidth),vec(bendReflectedfluxPlaneXPos,wgHorYCen+wgWidth))
bendReflectedFlux = bendField.add_dft_flux_plane(bendReflectedFluxPlane, srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 100)

#load the minused reflection flux from the "no bend" case as initalization
bendReflectedFlux.load_hdf5(bendField, "minusflux")


#define the flux plane for the transmitted flux


#define the flux plane for the transmitted flux
bendTransmfluxPlaneYPos = padSize + wgLengthY - 1.5; #the Y-coordinate of our transmission flux plane
bendTransmFluxPlane = volume(vec(wgVerXCen - wgWidth,bendTransmfluxPlaneYPos),vec(wgVerXCen + wgWidth,bendTransmfluxPlaneYPos))
bendTransmFlux = bendField.add_dft_flux_plane(bendTransmFluxPlane, srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 100 )

master_printf("Calculating...\n")
bendProbingpoint = vec(wgVerXCen,bendTransmfluxPlaneYPos) #the point at the end of the waveguide that we want to probe to check if source has decayed
runUntilFieldsDecayed(bendField, bendVol, srcComp, bendProbingpoint, pHDF5OutputFile = bendFileOutputEz)
master_printf("Done..!")

#construct 2-dimensional array with the flux plane data
#see python_meep.py
bendReflFlux = getFluxData(bendReflectedFlux)
bendTransmFlux = getFluxData(bendTransmFlux)

del_EPS_Callback()

#SHOW THE RESULTS IN A PLOT
frequencies = bendReflFlux[0] #should be equal to bendTransmFlux.keys() or noBendTransmFlux.keys() or ...
Ptrans = [x / y for x,y in zip(bendTransmFlux[1], noBendTransmFlux[1])]
Prefl = [ abs(x / y) for x,y in zip(bendReflFlux[1], noBendTransmFlux[1]) ]
Ploss = [ 1-x-y for x,y in zip(Ptrans, Prefl)]

matplotlib.pyplot.plot(frequencies, Ptrans, 'bo')
matplotlib.pyplot.plot(frequencies, Prefl, 'ro')
matplotlib.pyplot.plot(frequencies, Ploss, 'go' )

matplotlib.pyplot.show()
