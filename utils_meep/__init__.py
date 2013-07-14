#from meep import *
import images2gif
import h5py
from numpy import shape,transpose,array,dstack,uint8,uint,ones,abs
from PIL import Image
from matplotlib import cm
import os
from meep import prepareHDF5File,Ez

def run_and_gif(fields,my_space,steps = 1000,frames = 40,filename = "default.gif"):
    filename = filename[:-3] + "h5"
    the_file = prepareHDF5File(filename)
    for i in range(frames):
        for j in range(steps/frames):
            fields.step()
        fields.output_hdf5(Ez,my_space.surroundings(),the_file,1)
    del the_file
    make_gif(filename,dielectric_bg = None,normalize = False)  

def make_gif(h5name,dielectric_bg = True,data_name = None,normalize = False):
    """makes a gif animation of the real part of the field in the file h5name.
    if dielectric_bg is a h5 name, will also plot the dielectric landscape. in case dielectric_bg = True
    will look for a file dielectric.h5 as well"""
    f = h5py.File(h5name)
    if data_name == None:
        data_name = f.keys()[-1]
    d = f[data_name]
    if normalize:
        d = d/abs(array(d).flatten()).max()
    if dielectric_bg == True:
        where = "./" + os.path.dirname(h5name) + "/dielectric.h5"
        if os.path.exists(where):
            dielectric_bg = where
        else:
            dielectric_bg = None
    if dielectric_bg != None:
        diel = h5py.File(dielectric_bg)
        diel_val = transpose(array(diel["eps"]))
        imdiel = diel_val
        diel.close()
    else:
        imdiel = transpose(ones(shape(d[:,:,0])))
        #m = max(imdiel.flatten())
    di = uint8(dstack((imdiel,imdiel,imdiel,imdiel)))
    images2gif.writeGif(h5name[:-3] + ".gif",[Image.fromarray(cm.jet(0.5+transpose(d[:,:,i]), bytes=True)*(di==1.0) + cm.gray(0.5+transpose(d[:,:,i]) , bytes=True)*(di!=1.0)) for i in range(shape(d)[-1])])
    f.close()


def gif_mode(make_fields,my_space,bloch_k,f,df,steps_per_frame = 5,symmetry = 1.0):
    """runs the field for some time, then record for one period.
    Make sure that the field was initialized with the right freq source k and df"""
#    from Meep.blochCristal import make_fields
    name = "mode_k_" + str(bloch_k) +"_f_"+str(f)+"_df_"+str(df)+".h5"
    
    the_file = prepareHDF5File(name)
    
    fields = make_fields(bloch_k = bloch_k,f = f,df = df,symmetry = symmetry)

    t = 0.0
    time_source = 1.0/df
    while(t<time_source * 10):
        fields.step()
        t = fields.time()
    t_end = t + 1./f
    i=0
    n_frames = 0
    while(t<t_end):
        t = fields.time()
        i = i+1
        fields.step()
        if i==steps_per_frame:
            fields.output_hdf5(Ez,my_space.surroundings(),the_file,1)
            i=0
#            n_frames+=1       
    del the_file
    make_gif(name,normalize = True)


