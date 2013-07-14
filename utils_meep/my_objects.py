from meep import runWithHarminv,Dielectric,Ez,Ex,Ey,X,Y,Z,Low,High,Periodic,Callback,vec,getFluxData,pml,volume,set_EPS_Callback,identity,mirror,structure,voltwo,EPS,fields,gaussian_src_time,prepareHDF5File,vol3d
import h5py
from PIL import Image
import images2gif
from matplotlib import cm
from numpy import transpose,max,array
from pylab import imshow
from utils_meep import make_gif

class My_flux:
    def __init__(self,name,vol,f_start,f_stop,n_freqs):
        self.name = name
        self.vol = vol
        self.f_start = f_start
        self.f_stop = f_stop
        self.n_freqs = n_freqs
        
    def add_to_fields(self,fields):
        self.meep_flux = fields.add_dft_flux_plane(self.vol,self.f_start, self.f_stop, self.n_freqs)


class My_abstract_source:
    def __init__(self,f,df,comp = Ez):
        self.f = f
        self.df = df
        self.comp = comp
        
    def get_time_src(self):
        return gaussian_src_time(self.f,self.df)
        
class My_volume_source(My_abstract_source):
    def __init__(self,vol,*args):
        My_abstract_source.__init__(self, *args)
        self.vol = vol
    
    def get_source_geo(self):
        return self.vol
    
    def add_to_fields(self,fields):
        fields.add_volume_source(self.comp, self.get_time_src(), self.get_source_geo())
    
    
        
class My_gaussian_beam(My_volume_source):
    def __init__(self,waist,*args):
        My_volume_source.__init__(self, *args)
    
    
    def add_to_fields(self,fields):
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
        fields.add_volume_source(self.comp, self.get_time_src(), self.get_source_geo())
        
        

class My_point_source(My_abstract_source):
    def __init__(self,center,*args):
        My_point_source.__init__(self,*args)
        self.center = center
    def add_to_fields(self,fields):
        fields.add_point_source(self.comp, self.get_time_src(),self.center)
        
    
class My_flat_source(My_volume_source):
    pass
    
    

class My_space:
    def __init__(self,gridSizeX,gridSizeY,res,f,df,n_freqs = 1000,boundary_conditions = pml(1.0,Y),periodic_directions = [],bloch = vec(0.0,0.0),symmetry_direction = None,symmetry_val = complex(1.0),gridSizeZ = None):
        self.gridSizeX = gridSizeX
        self.gridSizeY = gridSizeY
        self.gridSizeZ = gridSizeZ
        self.res = res
        self.boundary_conditions = boundary_conditions
        self.fluxes = []
        self.f = f
        self.df = df
        self.n_freqs = n_freqs
        self.symmetry_direction = symmetry_direction
        self.symmetry_val = symmetry_val
        self.periodic_directions = periodic_directions
        self.bloch = bloch
        self.my_source = None

    def dim(self):
        if self.gridSizeZ is not None:
            return 3
        return 2
    
    def vol_cut(self):
        if self.dim() == 2:
            return self.meep_space.surroundings()
        if self.dim() == 3:
            return volume(vec(0,self.gridSizeY/2.,0),vec(self.gridSizeX,self.gridSizeY/2.,self.gridSizeZ)) 

    def add_structure(self,structure):
        self.my_structure = structure
        
    def add_flat_source(self,vol = volume(vec(0,1.0),vec(1.0,1.0)),comp = Ez):
        self.my_source = My_flat_source(vol,self.f,self.df,comp)
        
    def add_a_flux(self,name,vol):
        self.fluxes.append(My_flux(name,vol,self.f-self.df/2.0,self.f + self.df/2.0,self.n_freqs))
    
    def get_a_flux(self,name):
        for i in self.fluxes:
            if i.name == name:
                return array(getFluxData(i.meep_flux)[1])
    
    def show_dielectric(self):
        the_file = prepareHDF5File("coucou.h5")
        self.meep_fields.output_hdf5(Dielectric,self.meep_space.surroundings(),the_file)
        del the_file
        
        f = h5py.File("coucou.h5")
        d = f["eps"]
        imshow(transpose(d))
        f.close()
        
    def get_freqs(self):
        return array(getFluxData(self.fluxes[0].meep_flux)[0])
        
    def make_movie(self,steps = 5000,frames = 100,filename = "movie.gif"):
        fields = self.make_fields()
        h5filename = filename[:-4] + ".h5"
        the_file = prepareHDF5File(h5filename)
        for i in range(frames):
            for j in range(steps/frames):
                fields.step()
            fields.output_hdf5(self.my_source.comp,self.vol_cut(),the_file,1)
        del the_file
        self.save_dielectric()
        
        make_gif(h5filename,normalize = False)  
        
    def save_dielectric(self,filename = "dielectric.gif"):
        if filename[-4:]!=".gif":
            filename = filename + ".gif"
            
        the_file = prepareHDF5File("dielectric.h5")
        
        
        self.meep_fields.output_hdf5(Dielectric,self.vol_cut(),the_file)
        del the_file
        
        f = h5py.File("dielectric.h5")
        d = f["eps"]
        d = d/max(array(d).flatten()) 
        images2gif.writeGif(filename,[Image.fromarray(cm.jet(transpose(d[:,:]), bytes=True))])
        f.close()
    
    
    def make_fields(self):
        if self.gridSizeZ == None:
            self.meep_space = voltwo(self.gridSizeX,self.gridSizeY,self.res)
        else:
            self.meep_space = vol3d(self.gridSizeX,self.gridSizeY,self.gridSizeZ,self.res)
        material = epsilon(self.my_structure)
        set_EPS_Callback(material.__disown__())
        
        
        if self.symmetry_direction == None:
            sym = identity()
        else:
            sym = mirror(self.symmetry_direction,self.meep_space)*self.symmetry_val
            
        self.meep_structure = structure(self.meep_space, EPS, self.boundary_conditions,sym)
        
        
        the_fields = fields(self.meep_structure)
        
        for direc in self.periodic_directions:
            the_fields.set_boundary(Low,direc,Periodic)
            the_fields.set_boundary(High,direc,Periodic)
        
        
        the_fields.use_bloch(self.bloch)
        
        if self.my_source is not None:
            self.my_source.add_to_fields(the_fields)
     
        for f in self.fluxes:
            f.add_to_fields(the_fields)
        
        self.meep_fields = the_fields
        return the_fields
    

    def make_harminv(self,probing = vec(1.0,1.0),n_points = 15):
        the_field = self.make_fields()
        return runWithHarminv(the_field,self.meep_space,self.my_source.comp,probing,self.f,self.df,15)#pHDF5OutputFile = my_file,pHDF5OutputFilePhase3 = other_file))
    

class epsilon(Callback):
        def __init__(self,the_object):
            Callback.__init__(self)
            self.object = the_object
            
        def double_vec(self, r):
            ## Reference option (wave passing through air)
            if self.object.ghost:
                return 1. 
            else:
                return self.object.double_vec(r)
            

class My_structure:
    """an object without sources and without space"""
    def __init__(self,centerX,centerY,centerZ = 0.0,ghost = False):
        self.ghost = ghost
        self.centerX = centerX
        self.centerY = centerY
        self.centerZ = centerZ
        
