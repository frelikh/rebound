from dynamics import wholeSimulation
import rebound
import numpy as np
from rebound import hash as h
import time as t
from multiprocessing import Pool
NUMBER_OF_PROCESSES = 5
NPLANETS = 40
mplan = 0.0005

def planet_distances(n):
	#exponents = -0.3+1.3*np.random.rand(n)
	#a = [10.0**item for item in exponents]
	a = 5.0*np.random.rand(n)
	return a

def radii(n):
	#r = np.zeros(n)
	r = [50.0*4.66e-4 for i in range(n)]
	return r

def masses(n):
	m = [mplan for i in range(n)]
	return m

def l(n):
	return 2.0*np.pi*np.random.rand(n)

n_runs = 100

all_distances = np.zeros((n_runs,NPLANETS))
all_ls = np.zeros((n_runs,NPLANETS))

for j in range(n_runs):
	all_distances[j] = planet_distances(NPLANETS)
	all_ls[j] = l(NPLANETS)


def multi_run_simulation(infile):
	sim1 = wholeSimulation(n=NPLANETS,x = all_distances[infile],rplanet = radii(NPLANETS),
		ls = all_ls[infile], mplanet = masses(NPLANETS))
	sim1.run_simulation(infile,20)
	return 0




files = []
[files.append(i) for i in range(n_runs)]

if __name__ == "__main__":
	pool = Pool(NUMBER_OF_PROCESSES)
	results = pool.map(multi_run_simulation,files)
	print(results)

g = open("all.txt","w")
g.write("Semimajor axis Ecentricity Inclination  # Particles   Mass\n")
g.write("----------------------------------------------------------\n")
