from utils_meep.my_objects import My_structure,My_space,My_flat_source
from meep import *
from HDnavigator import *
import h5py

class membrane1D(My_structure):
    def __init__(self,height = 1.2,width = 0.7,n = 3.0,centerX = 0.0,centerY = 0.0):
        My_structure.__init__(self,centerX,centerY)
        self.height = height
        self.width = width
        self.n = n
        
    def double_vec(self, r):
        if r.y() > self.centerY + self.height/2.:
            return 1.
        if r.y() < self.centerY - self.height/2.:
            return 1.

        #if r.x() > self.centerX + self.width/2.:
        #    return 1.
        #if r.x() < self.centerX - self.width/2.:
        #    return 1.
        
        return self.n**2
    
    
    
gridSizeX = 1
gridSizeY = 12
border = 1.0
res = 40
f = 2.0
df = 3.0
space = My_space(gridSizeX,gridSizeY, res, f, df)

membrane = membrane1D(height = 0.627, width = 0.62, n = 3.214, centerX = gridSizeX/2., centerY = gridSizeY/2.)
space.add_structure(membrane)


space.add_flat_source(comp = Ex)


space.add_a_flux("input",volume(vec(0,2*border),vec(gridSizeX,2*border)))
space.add_a_flux("output",volume(vec(0,gridSizeY-2*border),vec(gridSizeX,gridSizeY-2*border)))
##=============================================
##     RUN
##=============================================

#f = space.make_fields()



##=============================================
##   SAVEDIR
##=============================================
next_dir("just1layer")
space.make_movie(steps = 50000,filename = "run.gif")
inp = space.get_a_flux("input")
out = space.get_a_flux("output")

space.my_structure.ghost = True
space.make_movie(steps = 50000,filename = "cal.gif")

inp_cal = space.get_a_flux("input")
out_cal = space.get_a_flux("output")
freq = space.get_freqs()

refl = (inp_cal - inp)/inp_cal

plot(1.0/freq,out_cal/inp_cal,label = "transmission no structure")
plot(1.0/freq,out/inp_cal,label = "transmission with structure")
plot(1.0/freq,refl,label = "reflexion with structure")
        
plot(1.0/freq,(out+(inp_cal - inp))/inp_cal,label = "total")
space.my_structure.ghost = False

f = h5py.File("reflectivities.h5","w")

f["freq"] = freq
f["r"] = refl
f.close()
a = gca()
a.set_ybound(-0.1,1.2)
show()
gcf().savefig("reflectivities.pdf")
1/0
##================================================
##    Band diagram
##================================================
figure()
import h5py
f = h5py.File("modesBloch.h5","w")
for sym,col in zip([-1,+1],["r","b"]):
    results = []
    ks = arange(0.0,0.51,0.01)
    for k in ks:
        space.bloch = vec(k,0)
        space.symmetry_direction = Y
        space.symmetry_val = complex(sym,0)
        results.append(space.make_harminv(probing = vec(gridSizeX/2.12345,gridSizeY/2.12345),n_points = 15))
        
        #results.append(runWithHarminv(my_fields,space.meep_space,Ez,vec(gridSizeX/2.12345,gridSizeY/2.12345),1.0,1.0,15))#pHDF5OutputFile = my_file,pHDF5OutputFilePhase3 = other_file))
    
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
    
    plot([0,0.5,0.0],[0,0.5,1.0],"black")
    a = gca()
    a.set_xbound(0,0.5)
    a.set_ybound(0,1.8)
    legend(loc = "best")
    
    
    f["k_s_sym_" + str(sym)] = k_s
    f["w_s_sym" + str(sym)] = w_s
    f["q_s_sym" + str(sym)] = q_s
    show()
f.close()
gcf().savefig("bands.pdf")
##==============================================================
##    Add the modes on top
##==============================================================
figure()
f = h5py.File("reflectivities.h5")

freq = f["freq"]
r = f["r"]


plot(1.0/freq,r,"g")

f = h5py.File("modesBloch.h5")
ws = f["w_s_sym1"]
qs = f["q_s_sym1"]
ks = f["k_s_sym_1"]
was = f["w_s_sym-1"]
qas = f["q_s_sym-1"]
kas = f["k_s_sym_-1"]


ws  = [w for (w,k) in zip(ws , ks) if k==0]
qs  = [q for (q,k) in zip(qs , ks) if k==0]
was = [w for (w,k) in zip(was,kas) if k==0]
qas = [q for (q,k) in zip(qas , kas) if k==0]


def lorentz(center,gamma):
    x = linspace(center-5*g,center+5*g,100)
    return x,0.1*((g**2)/4)/((x-center)**2 + (g**2)/4)


for w,q in zip(was,qas):
    g = w/q
    x,y = lorentz(w,g)
    y = y + 1.0
    plot(1.0/x,y,"r")


for w,q in zip(ws,qs):
    g = w/q
    x,y = lorentz(w,g)
    y = y + 1.1
    plot(1.0/x,y,"b")

#scatter(ws,1.1*ones(len(ws)),s = sqrt(abs(array(qs)))/3,c = "b",label = "sym")
#scatter(was,1.1*ones(len(was)),s = sqrt(abs(array(qas)))/3,c = "r",label = "antisym")

a = gca()
a.set_ybound(-0.1,1.3)
legend(loc = "best")
show()

f.close()

gcf().savefig("mixture.pdf")

