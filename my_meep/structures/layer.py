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
gridSizeY = 12.0
border = 1.5
srcFreqCenter = 0.45 #gaussian source center frequency
srcPulseWidth = 0.7 #gaussian source pulse width

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
class epsilon(Callback):
    def cube(self,x_center,y_center,width = 1.0,height = 1.0):
        pol_points = numpy.zeros((5,2))
        pol_points[0] = [x_center - width/2, y_center - height/2]
        pol_points[1] = [x_center + width/2,y_center - height/2]
        pol_points[2] = [x_center + width/2, y_center + height/2]
        pol_points[3] = [x_center - width/2, y_center + height/2]
        pol_points[4] = [x_center - width/2, y_center - height/2]
        return pol_points
    
    def __init__(self,normalize = False):
        Callback.__init__(self)
        self.normalize = normalize
                
    def double_vec(self, r):
        ## Reference option (wave passing through air)
        if self.normalize: return 1. 
        
        width = 1.2
        if r.y() > gridSizeY/2 + width/2:
            return 1.
        if r.y() < gridSizeY/2 - width/2:
            return 1.
        return 13.
        
##===============================================================
##       DEFINE FIELDS
##===============================================================

def make_fields(calibre = False,bloch_k = 0.0,f = 0.25,df = 1.5,symmetry = None):
    material = epsilon(calibre)
    set_EPS_Callback(material.__disown__())
    if symmetry == None:
        sym = identity()
    else:
        sym = mirror(Y,my_space)*complex(symmetry,0)
    
    
#use_averaging(False)  #--> uncomment this line to disable eps-averaging!
    my_structure = structure(my_space, EPS, pml(1.0,Y),sym)
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
    srcGaussian = gaussian_src_time(f, df)
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
            result = 1.0#complex(exp(-(x/2)**2),0)
            return result


#    af = amplitudeFactor()
#    set_AMPL_Callback(af.__disown__())
    srcGeo = volume(vec(0.,border),vec(gridSizeX,border))
    the_fields.add_volume_source(Ez, srcGaussian, srcGeo)
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
#my_fields_cal = make_fields(calibre = True)