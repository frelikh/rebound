import rebound
import numpy as np
from rebound import hash as h
import time as t


class wholeSimulation(object):

	def __init__(self,x,rplanet,ls,mplanet,t_max=1.0e3,n=10,mstar=1.,rstar=4.6e-3,d_max=30.,
		Noutputs=1000):
		'''

		:param n: The number of planets
		:param mstar: The stellar mass
		:param rstar: The stellar radius
		:param x: The distances of planets array
		:param rplanet: The planet physical radii array
		:param ls: The mean longitudes of the planets array
		:param mplanet: The planet masses array
		:param t_max: The maximum time for simulation
		:param d_max: The maximum distance for planet to escape
		:param Noutputs: The number of "timesteps" between 0 and t_max that are output 

		'''

		# User inputs
		self.n = n
		self.mstar = mstar
		self.x = x
		self.rplanet = rplanet
		self.ls = ls
		self.mplanet = mplanet
		self.t_max = t_max
		self.rstar = rstar
		self.d_max = d_max
		self.Noutputs = Noutputs

		self.start_time = t.time()
		self.sim = rebound.Simulation()


	def setupSimulation(self):
		#sim = rebound.Simulation()
		self.sim.add(m = self.mstar,hash="Sun",r=self.rstar)
		[self.sim.add(m=self.mplanet[i],a=self.x[i],hash=i+1,r=self.rplanet[i]) for i in range(self.n)]
		self.sim.move_to_com()
		return

	def stability(self,semimajor_axis,intervals=20,stability_condition = 0.5):
		step = self.Noutputs/intervals
		a_interval_spread = np.zeros((self.n,intervals))
		for particle in range(self.n):
			for i in range(intervals):
				begin = int(i*step)
				end = int((i+1)*step)
				a_interval = semimajor_axis[particle][begin:end]
				if (np.mean(a_interval) !=0 ):
					a_interval_spread[particle][i] = max(np.abs((np.amax(a_interval)-
						np.mean(a_interval))/np.mean(a_interval)),np.abs((np.amin(a_interval)-
						np.mean(a_interval))/np.mean(a_interval)))
					'''
					if (a_interval_spread[particle][i]>stability_condition):
						print("Unstable: %f" %a_interval_spread[particle][i])
					else:
						print("Stable %f" %a_interval_spread[particle][i])
					'''
				else:
					a_interval_spread[particle][i] = 0
				mask_sum = 0
		for particle in range(self.n):
			errors = np.asarray(a_interval_spread[particle])
			mask = errors>stability_condition
			mask_sum += np.sum(mask)
		stability_value = (mask_sum==0)
		return stability_value

	def elapsed_time(self):
		return t.time() - self.start_time

	def remove_all_particles(self):
		while(len(self.sim.particles)):
			self.sim.remove(index=0)


	def run_simulation(self,infile,repeats):

		f = open("foo%d.txt" %infile,'w')

		self.setupSimulation()
		max_time = self.t_max
		self.sim.exit_max_distance = self.d_max
		Noutputs = self.Noutputs

		times = np.linspace(0,max_time*2.*np.pi*repeats,Noutputs*repeats)
		#times = np.linspace(0,max_time*2.*np.pi,Noutputs)

		self.sim.collision = 'direct'
		self.sim.collision_resolve = 'merge'
		
		semimajor_axes = np.zeros((self.n,Noutputs*repeats))
		niter = 0
		while(niter<Noutputs*repeats):

			time = times[niter]
			try:
				self.sim.integrate(time)
			except rebound.Escape as error:
				for j in range(self.sim.N):
					p = self.sim.particles[j]
					d2 = p.x*p.x+p.y*p.y+p.z*p.z
					if d2>self.sim.exit_max_distance**2:
						index=j
				self.sim.remove(index=index)
				self.sim.move_to_com()
		
			for k in range(self.n):
				try:
					semimajor_axes[k][niter] = self.sim.particles[h(k+1)].a
				except rebound.ParticleNotFound as error:
					semimajor_axes[k][niter] = 0
			if((niter+1)>999 and (niter+1)%1000==0):
				new_a_values = [item[(niter-999):(niter+1)] for item in semimajor_axes]
				if(self.stability(new_a_values)):
					print("System is stable, exiting.")
					break
			if((niter+1)==(Noutputs*repeats)):
				print("System never became stable.")
			niter += 1

		[f.write("%d %2.3f %2.3f %2.3f %d %2.4f\n"% (infile, item.a,item.e,item.inc,
			self.sim.N,item.m)) for item in self.sim.particles[1:]]
		print("Number of particles remaining: {}".format(self.sim.N))
		f.close()
		return 0

