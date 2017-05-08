import rebound
import numpy as np
from rebound import hash as h

f = open('foo.txt','w')
print(f)
f.write("index a    e    i    energy_error    n_final\n")


def rhill(m,mstar,a):
    return a*(m/(3.*mstar))**(1./3.)

def planet_distances(n,m,mstar):
    a = np.zeros(n)
    a[0] = 0.1 # the closest orbit is 0.1 AU
    r_hill_0 = rhill(m,mstar,a[0])
    #print(r_hill_0)
    for i in range(1,n):
        a[i] = a[i-1] + 3.3*rhill(m,mstar,a[i-1])
        #print(a[i])
    return a


def setupSimulation(n,m,mstar):
    sim = rebound.Simulation()
    #sim.integrator = "whfast"
    sim.add(m=mstar,hash="Sun",r=4.6e-3)
    
    # specify random positions along the orbit:
    # use l, the mean longitude

    # generate random angles from 0 to 2pi, for each planet:
    ls = 2.0*np.pi*np.random.rand(n)

    a = planet_distances(n,m,mstar)
    [sim.add(m=m,a=a[i],l=ls[i],hash=i+1,r=(50.0*4.66e-4)) for i in range(n)]

    
    #sim.add(m=0.001,a=0.5,hash=1)
    #sim.add(m=0.001,a=1.,hash=2)
    #sim.add(m=0.001,a=1.5,hash=3)
    #sim.add(m=0.001,a=2.,hash=4)
    #sim.add(m=0.001,a=2.5,hash=5)
    sim.move_to_com()
    return sim

def energy_offset():
    dE = abs((sim.calculate_energy()-E0)/E0)
    return dE

from rebound import hash as h


for sim_count in range(20):
    sim = setupSimulation(5,0.001,1.)
    #sim.status()
    #print(sim.particles[h(1)].x)
    
    infile = "/home/renata/rebound/examples/sim_bins/sa%d.bin" %sim_count
    print(infile)
    sim.initSimulationArchive(infile,interval=1000.)

    max_time = 1.0e2

    sim.exit_max_distance=30.
    Noutputs=1000
    times=np.linspace(0,max_time*2.*np.pi,Noutputs)
    #sim.dt = 1. # timestep

    sim.collision="direct"
    sim.collision_resolve="merge"

    sim.track_energy_offset = 1
    E0 = sim.calculate_energy()
    #print("Energy: {}".format(E0))


    x1,y1=np.zeros(Noutputs),np.zeros(Noutputs)
    x2,y2=np.zeros(Noutputs),np.zeros(Noutputs)
    x3,y3=np.zeros(Noutputs),np.zeros(Noutputs)
    x4,y4=np.zeros(Noutputs),np.zeros(Noutputs)
    x5,y5=np.zeros(Noutputs),np.zeros(Noutputs)

    a1 = np.zeros(Noutputs)
    a2 = np.zeros(Noutputs)
    a3 = np.zeros(Noutputs)
    a4 = np.zeros(Noutputs)
    a5 = np.zeros(Noutputs)

    semimajor_axes = [a1, a2, a3, a4, a5]

# note: this is how they do it in the examples, later maybe re-write it with purely hashes
# i.e. loop through the particles like: 
# ds = [(item.x**2.0+item.y**2.0+item.z**2.0) for item in sim.particles]
# mask = ds>sim.exit_max_distance**2
# for mask_index,count in enumerate(mask):
# if count:
# sim.remove(index=mask_index)
# (to do)

    for i,time in enumerate(times):
        try:
            sim.integrate(time)
        except rebound.Escape as error:
            #print(error)
            #print(sim.t/(2.0*np.pi))
            for j in range(sim.N):
                p = sim.particles[j]
                d2 = p.x*p.x+p.y*p.y+p.z*p.z
                if d2>sim.exit_max_distance**2:
                    index=j
                    #print(d2)
            sim.remove(index=index)
            sim.move_to_com()
    en_error = energy_offset()
    close_encounters = sim.ri_hermes._steps_miniactive
    [f.write("%d %2.3f %2.3f %2.3f %f %d %1.1f\n"% (sim_count,item.a,item.e,item.inc,en_error,sim.N,close_encounters)) for item in sim.particles[1:]]
    #[f.write("%d %2.3f %2.3f %2.3f %f %d\n"% (h(item.hash).value,item.a,item.e,item.inc,en_error,sim.N)) for item in sim.particles[1:]]
    print("{},{}\n".format(sim.particles[0].x,sim.particles[0].y))
    print("Particles left: {}".format(sim.N))
    #print(sim_count)

        #print(sim.calculate_energy())

        # loop through 5 particles
        # if they're still left, then the 'a' value is recorded
        # otherwise, a 0 is recorded
'''
        for k in range(5):
            try:
                semimajor_axes[k][i] = sim.particles[h(k)].a
            except rebound.ParticleNotFound as error:
                semimajor_axes[k][i] = 0
'''    

        #print("Number of particles remaining: {}".format(sim.N))
        #b = [item.hash for item in sim.particles]
        #print("Hashes of particles remaining: {}".format(b))
        #c = [item.a for item in sim.particles[1:]]
        #print("Semimajor axis:{}".format(c))
        #sim.status()


    #from matplotlib import pyplot as plt
f.close()
#plt.scatter(times,semimajor_axes[1])
#print(semimajor_axes[1])
#fig,ax = plt.subplots()


#[ax.plot(times/(2.0*np.pi),item,label=i) for i,item in enumerate(semimajor_axes)]
#plt.legend()

#ax.set_xlabel("Years")
#ax.set_ylabel("semimajor axis")
#ax.set_yscale('log')

#plt.savefig("distances.png")

