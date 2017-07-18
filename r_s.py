import rebound
import numpy as np
from matplotlib import pyplot as plt
bin_dir = '/home/renata/rebound/examples/stability_tests'
indir = '/home/renata/rebound/examples'
stability_condition = 0.5
dividing_chunks = 50
repeats = 10

from rebound import hash as h

def planet_distances(n,m,mstar):

    a = 0.3*np.random.rand(n)
    print(a)
    return a
    
def setupSimulation(n,m,mstar):
    sim = rebound.Simulation()
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

def stability(semimajor_axis,intervals):
	step = Noutputs/intervals
	a_interval_spread = np.zeros((number_of_particles,intervals))

	for particle in range(number_of_particles):
		for i in range(intervals):
			begin = int(i*step)
			end = int((i+1)*step)
			#print(begin,end)
			#print("begin,end: %d %d" %(begin,end))
			#print("semimajor_axis")
			#print(semimajor_axis)
			a_interval = semimajor_axis[particle][begin:end]
			if (np.mean(a_interval) != 0):
				#a_interval_spread[i] = np.absolute((np.std(a_interval))/np.mean(a_interval))
				a_interval_spread[particle][i] = max(np.abs((np.amax(a_interval)-np.mean(a_interval))/np.mean(a_interval)),np.abs((np.amin(a_interval)-np.mean(a_interval))/np.mean(a_interval)))
				if (a_interval_spread[particle][i]>stability_condition):
					print("Unstable: %f" %a_interval_spread[particle][i])
				else:
					print("Stable %f" %a_interval_spread[particle][i])
			else:
				a_interval_spread[particle][i] = 0

	mask_sum = 0
	for particle in range(number_of_particles):
		errors = np.asarray(a_interval_spread[particle])
		print("Errors are:")
		print(errors)
		mask = errors>stability_condition
		mask_sum += np.sum(mask)
		print("Particle %d" %particle)
		print("Greater than %2.2f if unstable: %d" %(stability_condition,mask_sum))
	stability_value = (mask_sum==0)
	# true if stable, false if unstable
	return stability_value



number_of_particles = 3
Noutputs=1000

semimajor_axes = np.zeros((number_of_particles,Noutputs*repeats))
#eccentricites = np.zeros((number_of_particles,Noutputs))
#inclinations = np.zeros((number_of_particles,Noutputs))
#Omegas = np.zeros((number_of_particles,Noutputs))
#pomegas = np.zeros((number_of_particles,Noutputs))
#omegas = np.zeros((number_of_particles,Noutputs))
#fs = np.zeros((number_of_particles,Noutputs))
#Ms = np.zeros((number_of_particles,Noutputs))
#particle_masses = np.zeros((number_of_particles,Noutputs))
        
import time as t
def elapsed_time(start_time):
    return t.time() - start_time 



sim = setupSimulation(number_of_particles,0.001,1.0)
#sim.initSimulationArchive("bin_dir/archive.bin", interval=1e3)
    
max_time = 200
sim.exit_max_distance=30.0


initial_times =np.linspace(0,max_time*2.*np.pi,Noutputs)
all_times = np.linspace(0,max_time*2.*np.pi*repeats,Noutputs*repeats)
start_time = t.time()

def run_sim(sim,times):
	niter = 0


	while(True):
		print("Run number: %d" %niter)
		# integrate to this point in time, iterate
		# through the timesteps on the times array
		time = times[niter]
		#print("time")
		#print(time)
		try:
			sim.integrate(time)
		# exception if a particle escapes,
		# remove it, and repeat that integration
		# step without the escaping particle
		except rebound.Escape as error:
			print(error)
			print("Time_elapsed = %f s\n" %elapsed_time(start_time))
			for j in range(sim.N):
				p = sim.particles[j]
				d2 = p.x*p.x+p.y*p.y+p.z*p.z
				if d2>sim.exit_max_distance**2:
					index=j
					#print(d2)
			sim.remove(index=index)
			sim.move_to_com()

		# record the semimajor axes for the particles at each timestep
		for k in range(number_of_particles):
			try:
				semimajor_axes[k][niter] = sim.particles[h(k+1)].a
				#print(semimajor_axes[k][niter])
			except rebound.ParticleNotFound as error:
				semimajor_axes[k][niter] = 0  

		# check for stability, but only at intervals of 1000
		# record the semimajor axes at a particlular time
		
		if(niter>999 and niter%1000==0):
			new_a_values = [item[(niter-1000):niter] for item in semimajor_axes]
			#print(new_a_values)
			if(stability(new_a_values,dividing_chunks)):
				print("System is stable, exiting.")
				break
		niter += 1
	[plt.plot(all_times, item) for item in semimajor_axes]
	plt.savefig(indir+"/image.png")
	plt.clf()
		
	return 0


run_sim(sim,all_times)


    