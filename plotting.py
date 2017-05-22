import numpy as np

infiles = ["/home/renata/rebound/examples/foo%d.txt" %i for i in range(10)]
infile = infiles[9]

planet_masses = 0.001*np.linspace(0.5,10.,num=10)
planet_mass = planet_masses[9]

a,e,m = np.loadtxt(infile, usecols=(1,2,7),unpack=True,skiprows=1)
from matplotlib import pyplot as plt
mask = m>0.001
plt.scatter(a,e,s=5,label="m<MJ",color='r')
plt.scatter(a[mask],e[mask],color='b',s=5,label="m>MJ")
plt.xscale('log')
#plt.legend()

#semimajor_axis = np.linspace(0.01,2.,20)
semimajor_axis = np.logspace(-2.0,1.0,100)

G = 6.67e-8
Mpl = 2.0e30
Msun = 2.0e33
AU = 1.5e13
Rj = 50.0*7.0e9

#planet_radii = [7.0e9*item for item in range(5,56,10)]
#planet_radius = 7.0e9*50.0

def v_esc(mass):
        return (2.0*G*Msun*mass/Rj)**0.5

#v_esc = (2.0*G*Mpl/Rj)**0.5
def v_kep(a):
	omega = (G*Msun/(a*AU)**3.0)**0.5
	return a*AU*omega
print(v_esc)

es = [v_esc(planet_mass)/v_kep(item) for item in semimajor_axis]

plt.plot(semimajor_axis,es,label="%f M_sun" %planet_mass)

plt.legend()
plt.ylim(0,1.1)
plt.xlabel("Semimajor axis, AU")
plt.ylabel("Eccentricity")
plt.title("IC: Random a=0-5AU & true anomaly=0-2pi, m=%f Msun, 100 runs, n=10" %planet_mass)
plt.savefig("%fMsun.png" %planet_mass)
