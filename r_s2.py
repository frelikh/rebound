
import rebound
import numpy as np
#from matplotlib import pyplot as plt
#bin_dir = '/home/renata/rebound/examples/stability_tests'
# this for graymalkin:
#indir = '/home/renata/rebound/examples'
# for my computer:
indir = '/Users/renata/sftp_graymalkin'
stability_condition = 1.0
dividing_chunks = 50
repeats = 10

max_time = 1000
max_distance = 30.0

number_of_particles = 10
Noutputs=1000
max_niter = (repeats*max_time)

from rebound import hash as h

def planet_distances(n,m,mstar):
	#a = [0.5,1.0,1.5]

	a = 5.0*np.random.rand(n)
	#print(a)
	return a
	
def setupSimulation(n,m,mstar):
	sim = rebound.Simulation()
	sim.add(m=mstar,hash="Sun",r=4.6e-3)


	a = planet_distances(n,m,mstar)
	ls = 2.0*np.pi*np.random.rand(n)

	ms = [0.0003,0.0003,0.0003,0.0003,0.0003,0.003,0.003,0.003,0.003,0.003]
	[sim.add(m=ms[i],a=a[i],l=ls[i],hash=i+1,r=(50.0*4.66e-4)) for i in range(n)]

	sim.move_to_com()
	

	return sim

def stability(semimajor_axis,intervals):
	step = Noutputs/intervals
	a_interval_spread = np.zeros((number_of_particles,intervals))

	for particle in range(number_of_particles):
		for i in range(intervals):
			begin = int(i*step)
			end = int((i+1)*step)

			a_interval = semimajor_axis[particle][begin:end]
			if (np.mean(a_interval) != 0):
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
		#print("Errors are:")
		#print(errors)
		mask = errors>stability_condition
		mask_sum += np.sum(mask)
		#print("Particle %d" %particle)
		#print("Greater than %2.2f if unstable: %d" %(stability_condition,mask_sum))
	stability_value = (mask_sum==0)
	# true if stable, false if unstable
	return stability_value




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



initial_times =np.linspace(0,max_time*2.*np.pi,Noutputs)
all_times = np.linspace(0,max_time*2.*np.pi*repeats,Noutputs*repeats)
start_time = t.time()

# removes all particles from the simulation to clear
# values

def remove_all_particles(given_sim):
	while(len(given_sim.particles)):
		given_sim.remove(index=0)

def run_sim(sim,times,run_number,sim_label):
	niter = 0
	while(niter<(max_niter+1)):
		#print("Run number: %d" %niter)
		# integrate to this point in time, iterate
		# through the timesteps on the times array
		time = times[niter]
		print(run_number,niter)

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
				remove_all_particles(sim)
				break
				
			if(niter==10000):
				print("%d_%d never became stable" %(run_number,sim_label))
				remove_all_particles(sim)
		
		niter += 1
	t = range(Noutputs*repeats)

	#automatically closes the file
	with open(indir+"/semimajor_axis_file%d_%d.txt" %(run_number,sim_label),'w') as semimajor_axis_file:
		for row_number in range(Noutputs*repeats):
			for p_number in range(number_of_particles):
				semimajor_axis_file.write("%f %f %d\n" %(t[row_number],semimajor_axes[p_number][row_number],p_number))

	semimajor_axis_file.closed

	return 0

def run_simulations(number_of_runs, t_max, number_of_particles, m_particles, m_stellar,sim_label,repeats):
	# open file to save data
	f = open(indir+'/foo%d.txt' %sim_label,'w')
	print(f)
	f.write("index   a    e    i    energy_error    n_final   close_encounters  particle_mass\n")

	all_times=np.linspace(0,t_max*2.*np.pi*repeats,Noutputs*repeats)

	for sim_count in range(number_of_runs):
	
		sim = setupSimulation(number_of_particles,m_particles,m_stellar)


		sim.exit_max_distance=max_distance


		sim.collision="direct"
		sim.collision_resolve="merge"
		close_encounters = sim.ri_hermes._steps_miniactive

		run_sim(sim,all_times,sim_count,sim_label)

		# after each simulation ends, write final orbital elements to file
		[f.write("%d %2.3f %2.3f %2.3f %d %1.1f %2.4f\n"% (sim_count,item.a,item.e,item.inc,sim.N,close_encounters,item.m)) for item in sim.particles[1:]]

		print(sim_count)

		print("Number of particles remaining: {}".format(sim.N))
		


	f.close()
	return 0

run_simulations(100, max_time, number_of_particles , 0.001, 1.0, 0,repeats)
	
