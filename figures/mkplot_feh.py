import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import q2
from scipy.optimize import curve_fit
from matplotlib.ticker import FormatStrFormatter

def linear(x, m, b):
     model = m*x + b
     return model
     
root_dir = '/Users/mbedell/Documents/Research/HARPSTwins/Abundances/All/'
a = genfromtxt(root_dir+'final_abundances.csv', delimiter=',', dtype=None, names=True)
#data = q2.Data(root_dir+"final_parameters.csv",root_dir+"harpstwins_lines.csv")
par = np.genfromtxt(root_dir+"final_parameters.csv", delimiter=',', dtype=None, names=True)
gce = genfromtxt(root_dir+'GCE/gce_linear.txt', delimiter=',', dtype=None, names=True)
ages = np.genfromtxt(root_dir+'final_ages_combination.csv', delimiter=',', dtype=None, names=True)

age = ages['age_mean']
age_err = ages['age_std']
xs = np.arange(11.)

feh = par['feh'][:-1]
feh_err = par['err_feh'][:-1]

fit = [i not in ['HIP19911', 'HIP108158', 'HIP109821', 'HIP115577', 'HIP14501', 'HIP28066', 'HIP30476',
        'HIP33094', 'HIP65708', 'HIP73241', 'HIP74432', 'HIP64150'] for i in a['id'][:-1]] # mask out SB2, thick-disk
inv = np.invert(fit)

c2 = '#003399' # blue
c3 = '#CC0033' # red
c4 = '#339900' # green
plt.rcParams["font.sans-serif"] = "Helvetica"

fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(age[inv], feh[inv], xerr=age_err[inv], yerr=feh_err[inv], fmt='D', c=c3, ecolor=c3, ms=7)
ax.errorbar(age[fit], feh[fit], xerr=age_err[fit], yerr=feh_err[fit], fmt='o', c='black', ecolor='black', mec='black', ms=7)
ax.annotate(r'$\odot$', xy=(4.6, 0.0), horizontalalignment='center', verticalalignment='center', color=c4, fontsize=24, weight='bold')

#popt, pcov = curve_fit(linear, feh, age, sigma=age_err)
#fehs = (xs - popt[1])/popt[0]
#ax.plot(xs, fehs, c=c2)

ax.set_xlabel('Age (Gyr)')
ax.set_ylabel('[Fe/H] (dex)')

#fig.text(0.5, 0.02, 'Age (Gyr)', size=20, ha='center')
#fig.text(0.03, 0.5, '[Fe/H] (dex)', rotation=90, size=20, va='center')
fig.savefig('feh.png')