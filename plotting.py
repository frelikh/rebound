infile = "/home/renata/rebound/examples/foo.txt"
import numpy as np
a,e = np.loadtxt(infile, usecols=(1,2),unpack=True,skiprows=1)
from matplotlib import pyplot as plt
plt.scatter(a,e)
plt.xscale('log')

#semimajor_axis = np.linspace(0.01,2.,20)
semimajor_axis = np.logspace(-2.0,0.3,30)

G = 6.67e-8
Mpl = 2.0e30
Msun = 2.0e33
AU = 1.5e13
Rj = 50.0*7.0e9

v_esc = (2.0*G*Mpl/Rj)**0.5
def v_kep(a):
	omega = (G*Msun/(a*AU)**3.0)**0.5
	return a*AU*omega
print(v_esc)

es = [v_esc/v_kep(item) for item in semimajor_axis]
plt.scatter(semimajor_axis,es,color='r',s=2)
plt.savefig("foo.png")
