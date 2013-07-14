##===============================================================
##       IMPORTS
##===============================================================
from pylab import *
from numpy import *
from meep import *
import time
##===============================================================
##       NUMERICAL PARAMETERS
##===============================================================
res = 20
gridSizeX = 1.0
gridSizeY = 10.0
border = 1.0
srcFreqCenter = 1.5 #gaussian source center frequency
srcPulseWidth = 1.0 #gaussian source pulse width

##===============================================================
##       DEFINE SPACE
##===============================================================

#we create a structure with PML of thickness = 1.0 on all boundaries,
#in all directions,
#using the material function EPS
my_space = voltwo(gridSizeX,gridSizeY,res)

##===============================================================
##       DEFINE STRUCTURE
##===============================================================
class epsilon(PolygonCallback2D):
    def cube(self,x_center,y_center,width = 1.0,height = 1.0):
        pol_points = numpy.zeros((5,2))
        pol_points[0] = [x_center - width/2, y_center - height/2]
        pol_points[1] = [x_center + width/2,y_center - height/2]
        pol_points[2] = [x_center + width/2, y_center + height/2]
        pol_points[3] = [x_center - width/2, y_center + height/2]
        pol_points[4] = [x_center - width/2, y_center - height/2]
        return pol_points
    
    def __init__(self,normalize = False):
        PolygonCallback2D.__init__(self)
        master_printf("Callback function for epsilon contructed. Now defining the polygons...\n")
        #Points outside this polygon are assumed to have epsilon 1.0
        
        cubes = []
        fill_ratio = 0.56 #in units of lambda
        n_cubes = 1
        period = 0.7
        height = 0.35

        side = period*fill_ratio
        
        for i in range(n_cubes):
            cubes.append(self.cube(gridSizeX/2 - (n_cubes/2-i)*period,gridSizeY/2,side,height))
            
            
        if not normalize:
            for i in cubes:
                self.add_polygon(i, 3.5)
            master_printf("Polygons OK.\n")


##===============================================================
##       DEFINE FIELDS
##===============================================================

def make_fields(calibre = False,bloch_k = 0.0):
    material = epsilon(calibre)
    set_EPS_Callback(material.__disown__())
#use_averaging(False)  #--> uncomment this line to disable eps-averaging!
    my_structure = structure(my_space, EPS, pml(1.0,Y))
    the_fields = fields(my_structure)
    the_fields.set_boundary(Low,X,Periodic)
    the_fields.set_boundary(High,X,Periodic)
#    the_fields.set_boundary(Low,Y,Periodic)
#    the_fields.set_boundary(High,Y,Periodic)
    the_fields.use_bloch(vec(bloch_k,0))
    #the_fields.use_bloch()
##===============================================================
##       DEFINE SOURCE
##===============================================================
#define a gaussian line source of length 'wgWidth' at X=wgLength/2, Y=padSize
    srcGaussian = gaussian_src_time(srcFreqCenter, srcPulseWidth )
#    srcGeo = volume(vec(border,border),vec(gridSizeX-border,border))



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


#    af = amplitudeFactor()
#    set_AMPL_Callback(af.__disown__())
    the_fields.add_point_source(Ez, srcGaussian, vec(gridSizeX/2,gridSizeY/2))
    master_printf("Field created...\n")

    ##===============================================================
    ##            CREATE FLUX PLANE
    ##===============================================================
    y_low = gridSizeY/2-1.0
    y_up = gridSizeY/2+1.0
    
#    input_plane = volume(vec(border,y_low),vec(gridSizeX-border,y_low))
#    output_plane = volume(vec(border,y_up),vec(gridSizeX-border,y_up))
#    left_plane = volume(vec(border,y_low),vec(border,y_up))
#    right_plane = volume(vec(gridSizeX-border,y_low),vec(gridSizeX-border,y_up))

#    input_flux = the_fields.add_dft_flux_plane(input_plane,srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 100)
#    output_flux = the_fields.add_dft_flux_plane(output_plane,srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 100)
#    left_flux = the_fields.add_dft_flux_plane(left_plane,srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 100)
#    right_flux = the_fields.add_dft_flux_plane(right_plane,srcFreqCenter-(srcPulseWidth/2.0), srcFreqCenter+(srcPulseWidth/2.0), 100)
    

    return the_fields#,input_flux,output_flux,left_flux,right_flux
#my_file = prepareHDF5File("./coucou.h5")
#runUntilFieldsDecayed(my_fields, my_space, Ez, vec(2.0,2.0),pH5OutputIntervalSteps=1000,pHDF5OutputFile = my_file)


my_fields = make_fields()

#,input_flux,output_flux,left_flux,right_flux 
my_fields_cal = make_fields(calibre = True)
#,input_flux_cal,output_flux_cal,left_flux_cal,right_flux_cal

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
other_file = prepareHDF5File(savedir+"/other.h5")
steps_total = 500
n_frames = 200
steps_per_frame = steps_total/n_frames
#FileOutputComp = meep.prepareHDF5File(savedir + "/OutputComp.h5")
#ha = HarminvAnalysis(my_fields,0.5,2.0,10,500,Ez,vec(gridSizeX/2 + 0.1,gridSizeY/2 + 0.1),pOutputHDF5OnlyAfterAnalysis = True, pTimeAfterSources = 1000, pH5OutputIntervalSteps=6, pHDF5OutputFile = FileOutputComp))
results = []
ks = arange(0.0,1.0,0.05)
for k in ks:
    my_fields = make_fields(bloch_k = k)
    results.append(runWithHarminv(my_fields,my_space,Ez,vec(gridSizeX/2+1e-6,gridSizeY/2+1e-6),3.0,6.0,10, pH5OutputIntervalSteps = 6, pHDF5OutputFile = my_file,pHDF5OutputFilePhase3 = other_file))

if False:
    for i in range(steps_total):
        my_fields.step()
        if i%steps_per_frame == 0 :
            my_fields.output_hdf5(Ez, my_space.surroundings(),my_file,1)

del my_file
del other_file

#b = array((len(results),)
for i,k in zip(results,ks):
    for j in i:
        plot(k,float(real(j[0])),"o")
show()
##==============================================================
##      Harminv analysis
##==============================================================


##===============================================================
## IDEM FOR NORMALIZATION
##===============================================================
#material = epsilon(normalize = True)
#set_EPS_Callback(material.__disown__())

#my_structure_cal = structure(my_space, EPS, pml(1.0))
#my_fields_cal = fields(my_structure_cal)

#my_file = prepareHDF5File(savedir+"/imagesCalib.h5")
#my_fields.reset()

#for i in range(steps_total):
#    my_fields_cal.step()
#    if i%steps_per_frame == 0 :
#        my_fields_cal.output_hdf5(Ez, my_space.surroundings(),my_file,1)
#del my_file


##===============================================================
##            CONVERT H5 FILES TO PNG
##   see http://ab-initio.mit.edu/wiki/index.php/Meep_Tutorial
##===============================================================
my_file = prepareHDF5File(savedir+"/dielectric.h5")
my_fields.output_hdf5(Dielectric,my_space.surroundings(),my_file)
del my_file

##===============================================================
##            CONVERT H5 FILES TO PNG
##   see http://ab-initio.mit.edu/wiki/index.php/Meep_Tutorial
##===============================================================
import subprocess,glob
#for f in glob.glob(savedir +"/"+"*.h5"):
import os


sp = subprocess.Popen(["h5topng",savedir + "/dielectric.h5"])
sp.wait()


def make_gif(h5name,n_frames = n_frames,add_dielectric = False):
    curdir = os.getcwd()
    os.chdir(savedir)
    try:
        master_printf("building individual png...\n")
        t = time.time()
        sp = subprocess.Popen(["h5topng","-t 0:%i"%(n_frames-1),"-R","-Zc","dkbluered",h5name])
        sp.wait()
        
        name_pngs = h5name[:-3] + ".t*.png"
        if add_dielectric:
            name_pngs = "merged_" + name_pngs
            for f in glob.glob(h5name[:-3] + ".t*.png"):
                sp = subprocess.Popen(["convert",f,"dielectric.png","-compose", "Bumpmap","-composite","merged_" +f +".png"])
                sp.wait()
    

        master_printf("done in " + str(time.time()-t) + "\n")
        master_printf("making gif...\n")
        t = time.time()
        sp = subprocess.Popen(["convert",name_pngs, h5name[:-3] + ".gif"])
        sp.wait()
        master_printf("done in " + str(time.time()-t) + "\n")
        master_printf("done\n")
    except e:
        raise e
    finally:
        os.chdir(curdir)

make_gif("images.h5",add_dielectric = True)
make_gif("imagesCalib.h5")
make_gif("other.h5",7)

##===============================================================
##            PLOT FLUXES
##===============================================================
fig = figure()
down_cal = getFluxData(input_flux_cal)
up = getFluxData(output_flux)
down = getFluxData(input_flux)
left = getFluxData(left_flux)
right = getFluxData(right_flux)

freq = array(down[0])
trans = array(up[1])/array(down_cal[1])
refl = -(array(down[1])-array(down_cal[1]))/array(down_cal[1])
side = (array(right[1])-array(left[1]))/array(down_cal[1])
total = refl + trans + side
plot(freq,refl,label = "reflection")
plot(freq,trans,label = "transmission")
plot(freq,side,label = "side")
plot(freq,total,label = "total")

legend(loc = "best")

show()


#plot(diff(side_freq),label = "high_freq")
#legend()
#show()
#fig.savefig(savedir  +"/profiles.png")
