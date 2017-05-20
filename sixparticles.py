import rebound
import numpy as np
from rebound import hash as h

# change these values to specify simulation parameters
n_runs = 100
t_max = 1.0e3
#m_planet = 0.001
m_star = 1.
n_planets = 10
list_of_planet_ms = 0.001*np.linspace(0.5,10.,num=10)


# definition of hill radius
def rhill(m,mstar,a):
    return a*(m/(3.*mstar))**(1./3.)

# where to put the planets initially,
# could be random or separated by n*rhill
def planet_distances(n,m,mstar):

    a = 5.0*np.random.rand(n)
    print(a)
    
    '''
    a = np.zeros(n)
    a[0] = 0.1 # the closest orbit is 0.1 AU
    r_hill_0 = rhill(m,mstar,a[0])

    for i in range(1,n):
        a[i] = a[i-1] + 3.3*rhill(m,mstar,a[i-1])
    '''
    return a


# sets up simulation
# adds the Sun and the planets
# can choose # of planets(n), mass of planets(m), and
# mass of Sun

def setupSimulation(n,m,mstar):
    sim = rebound.Simulation()
    #sim.integrator = "whfast"
    sim.add(m=mstar,hash="Sun",r=4.6e-3)
    
    # specify random positions along the orbit:
    # use l, the mean longitude

    # generate random angles from 0 to 2pi, for each planet:
    ls = 2.0*np.pi*np.random.rand(n)

    a = planet_distances(n,m,mstar)

    # adding particles and identifying them with hashes
    # radius: 50 Jupiter radii (in AU)
    
    # random masses: next two lines. don't do this.
    #random_planet_masses = 0.0001+np.random.rand(n)*0.0019
    #[sim.add(m=random_planet_masses[i],a=a[i],l=ls[i],hash=i+1,r=(50.0*4.66e-4)) for i in range(n)]

    [sim.add(m=m,a=a[i],l=ls[i],hash=i+1,r=(50.0*4.66e-4)) for i in range(n)]

    # move to center of mass frame
    sim.move_to_com()
    return sim

# testing purposes: what is the change in energy from the initial value?
def energy_offset(sim,E0):
    dE = abs((sim.calculate_energy()-E0)/E0)
    return dE

def run_simulations(number_of_runs, t_max, number_of_particles, m_particles, m_stellar,sim_label):
    # open file to save data
    f = open('foo%d.txt' %sim_label,'w')
    print(f)
    f.write("index   a    e    i    energy_error    n_final   close_encounters  particle_mass\n")

    for sim_count in range(number_of_runs):
    
        sim = setupSimulation(number_of_particles,m_particles,m_stellar)
    

        # not necessary yet - this saves the simulations to binaries
        #infile = "/home/renata/rebound/examples/sim_bins/sa%d.bin" %sim_count
        #print(infile)
        #sim.initSimulationArchive(infile,interval=1000.)

        max_time = t_max

        sim.exit_max_distance=30.
        Noutputs=1000
        times=np.linspace(0,max_time*2.*np.pi,Noutputs)
        #sim.dt = 1. # timestep

        sim.collision="direct"
        sim.collision_resolve="merge"

        sim.track_energy_offset = 1
        E0 = sim.calculate_energy()
        #print("Energy: {}".format(E0))

        semimajor_axes = np.zeros((number_of_particles,Noutputs))

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

            # loop through the particles
            # if they're still left, then the 'a' value is recorded
            # otherwise, a 0 is recorded
            # record the semimajor axes of the particles
            # start with k=1, because we want planets and not the Sun
            # for the semimajor_axes, the particles start at 0
            for k in range(number_of_particles):
                try:
                    semimajor_axes[k][i] = sim.particles[h(k+1)].a
                except rebound.ParticleNotFound as error:
                    semimajor_axes[k][i] = 0  

        # for each simulation, output a plot of the distances (to check stability)
        from matplotlib import pyplot as plt
        fig,ax = plt.subplots()


        [ax.plot(times/(2.0*np.pi),item,label=i) for i,item in enumerate(semimajor_axes)]
        plt.legend()
        plt.title("{} particles".format(sim.N-1))
        ax.set_xlabel("Years")
        ax.set_ylabel("semimajor axis")
        ax.set_yscale('log')
        plt.ylim(0,10)

        plt.savefig("/home/renata/rebound/examples/a/a%d/distances%d.png" %(sim_label,sim_count))
        plt.close(fig)
        en_error = energy_offset(sim,E0)
        close_encounters = sim.ri_hermes._steps_miniactive

        # after each simulation ends, write final orbital elements to file
        [f.write("%d %2.3f %2.3f %2.3f %f %d %1.1f %2.4f\n"% (sim_count,item.a,item.e,item.inc,en_error,sim.N,close_encounters,item.m)) for item in sim.particles[1:]]
        #[f.write("%d %2.3f %2.3f %2.3f %f %d\n"% (h(item.hash).value,item.a,item.e,item.inc,en_error,sim.N)) for item in sim.particles[1:]]
        print("{},{}\n".format(sim.particles[0].x,sim.particles[0].y))
        #print("Particles left: {}".format(sim.N))
        print(sim_count)

        #print(sim.calculate_energy())



        print("Number of particles remaining: {}".format(sim.N))
        #b = [item.hash for item in sim.particles]
        #print("Hashes of particles remaining: {}".format(b))
        #c = [item.a for item in sim.particles[1:]]
        #print("Semimajor axis:{}".format(c))
        #sim.status()


    f.close()
    return 0

[run_simulations(n_runs, t_max, n_planets, m_planet, m_star,sim_label) for sim_label,m_planet in enumerate(list_of_planet_ms)]


#stability = [np.std(item[80:99]) for item in semimajor_axes]

#print(stability)
