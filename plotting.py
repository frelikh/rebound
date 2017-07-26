import numpy as np

G = 6.67e-8
Mpl = 2.0e30
Msun = 2.0e33
AU = 1.5e13
Rj = 50.0*7.0e9

infile = "/Users/renata/sftp_graymalkin/foo0.txt"
indeces,a,e,m = np.loadtxt(infile, usecols=(0,1,2,6),unpack=True,skiprows=1)

# group by max mass
for j in range(1,len(m)):
	if(indeces[j]==indeces[j-1]):
		m[j] = max(m[j],m[j-1])
print(m)

planet_mass = 0.001

def v_esc(mass):
        return (2.0*G*Msun*mass/Rj)**0.5

#v_esc = (2.0*G*Mpl/Rj)**0.5
def v_kep(a):
	omega = (G*Msun/(a*AU)**3.0)**0.5
	return a*AU*omega



from matplotlib import pyplot as plt
e_theoretical = [v_esc(m[i])/v_kep(a[i]) for i in range(len(m))]

# above the curve, colored blue
mask = e>e_theoretical

plt.scatter(a,e,s=5,label="e<e_theoretical",color='r')
plt.scatter(a[mask],e[mask],color='b',s=5,label="e>e_theoretical")
plt.xscale('log')
#plt.legend()

#semimajor_axis = np.linspace(0.01,2.,20)
semimajor_axis = np.logspace(-2.0,1.0,100)



#planet_radii = [7.0e9*item for item in range(5,56,10)]
#planet_radius = 7.0e9*50.0



es = [v_esc(planet_mass)/v_kep(item) for item in semimajor_axis]

plt.plot(semimajor_axis,es,label="%f M_sun" %planet_mass)

plt.legend()
plt.ylim(0,1.1)
plt.xlabel("Semimajor axis, AU")
plt.ylabel("Eccentricity")
plt.title("IC: Random a=0-5AU & true anomaly=0-2pi, m=%f Msun, 100 runs, n=10" %planet_mass)
plt.savefig("%fMsun.png" %planet_mass)
