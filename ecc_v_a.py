import numpy as n
indir = "/Users/renata/Eccentricities_vs_a/"
infile = indir+"sorted.txt"
data = n.genfromtxt(infile, skip_header = 2, delimiter = ",")

mask = n.isfinite(data) # test to see if number is not infinite and not a nan
# list of array indeces with rows that have both semimajor axis and eccentricity data 
# (and RV data - column 3 (i.e. fourth)) and star name
lst = []
for i in range(len(mask)):
    if (mask[i][0] and mask[i][1] and mask[i][2] and mask[i][3]) == True:
        lst.append(i)
# only keep data with both values for plotting
clean_data = [data[i] for i in lst]

G = 6.67e-8
Msun = 1.989e33
Mjup = 1.899e30
AU = 1.496e13
dens = 1.0 # solar density

def radius(m,rho):
    return (3.0*m/(4.0*n.pi*rho))**(1.0/3.0)

# m in solar masses, a in AU
def v_esc_sun(m,a):
    return (2.0*G*Msun*m/(a*AU))**0.5

# m in jupiter masses
def v_esc_planet(m,rho):
    mass = m*Mjup
    r = radius(mass,rho)
    return (2.0*G*mass/r)**0.5

from matplotlib import pyplot as plt
#%matplotlib inline
msini = [item[0] for item in clean_data]
a = [item[1] for item in clean_data]
ecc = [item[2] for item in clean_data]
star_name = [item[3] for item in clean_data]
mass_by_highest = n.zeros(len(msini))
mass_by_highest[0] = msini[0]

for i in range(1,len(star_name)):
    if(star_name[i]==star_name[i-1]):
        mass_by_highest[i]=msini[i-1]

a = n.asarray(a)
ecc = n.asarray(ecc)


mask2 = mass_by_highest>1.0

plt.scatter(a,ecc,label='Observed',s=2)
plt.scatter(a[mask2],ecc[mask2],color='r',s='2')


plt.xlabel("Semi-major axis")
plt.ylabel("Eccentricity")
plt.title("Exoplanet data (RV planets) from exoplanets.org")
plt.xscale("log")


a_theoretical = n.logspace(-2.0,1.0,num=100)
theoretical_e = [v_esc_planet(1.0,dens)/v_esc_sun(1.0,item) for item in a_theoretical]

plt.plot(a_theoretical,theoretical_e,color='red',label='Theoretical')
plt.xscale("log")
plt.ylim(0,1)

#mask2 = [msini>1.0]
#plt.scatter(a[mask2],ecc[mask2],label='>MJ',color='g')

plt.legend(loc=2)
plt.savefig("ecc_vs_a.png")