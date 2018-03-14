import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def linear(x, m, b):
     model = m*x + b
     return model

root_dir = '../data/'
par = np.genfromtxt(root_dir+'final_parameters.csv', delimiter=',', dtype=None, names=True, encoding=None)
ages = np.genfromtxt(root_dir+'final_ages_combination.csv', delimiter=',', dtype=None, names=True, encoding=None)

a = np.genfromtxt(root_dir+'final_abundances_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)  

thin = [i not in ['HIP19911', 'HIP108158', 'HIP109821', 'HIP115577', 'HIP14501', 'HIP28066', 'HIP30476',
                    'HIP33094', 'HIP65708', 'HIP73241', 'HIP74432'] for i in a['id'][:-1]] # mask out SB2, thick-disk
thick = np.invert(thin)  

adibekyan = np.genfromtxt('/Users/mbedell/Documents/Research/Stars/HARPS_GTO/Adibekyan2012/abundances.csv', dtype=None,
            delimiter=',', usecols=(9, 10, 16, 17), names=('mg', 'mg_err', 'si', 'si_err'), encoding=None)
#brewer = np.genfromtxt('/Users/mbedell/Documents/Research/Stars/Brewer/brewer2016.csv', dtype=None,
#            names=True, delimiter=',')

brewer = np.genfromtxt('/Users/mbedell/Documents/Research/Stars/Brewer/table_9_full_bedell.csv', dtype=None,
            names=True, delimiter=',', encoding=None)
bmask = (brewer['MH'] <= 0.15) & (brewer['MH'] >= -0.15) & \
        (brewer['logg'] >= 3.5)

#        (brewer['Teff'] <= 6100.) & (brewer['Teff'] >= 5500.) & \
brewer = brewer[bmask]
print("{0} stars selected from Brewer".format(np.sum(bmask)))


age = ages['age_mean']
age_err = ages['age_std']

# plot Mg/Si:

c2 = '#003399' # blue
c3 = '#CC0033' # red
c4 = '#339900' # green
plt.rcParams["font.sans-serif"] = "Helvetica"

fig = plt.figure(figsize=(15,12))

mgsi_sun = 1.05
mgsi = 10**(a["MgI_1"][:-1] - a["SiI_1"][:-1]) * mgsi_sun
mgsi_err = mgsi * np.log(10.) * np.sqrt(a["err_MgI"][:-1]**2 + a["err_SiI"][:-1]**2)

brewer_mgsi = 10**(brewer['MgH'] - brewer['SiH']) * mgsi_sun

ax3 = fig.add_subplot(223)

ax3.scatter(brewer['Age'], brewer_mgsi, marker='o', c='black', alpha=0.1)

ax3.errorbar(age[thin], mgsi[thin], yerr=mgsi_err[thin], xerr=age_err[thin], fmt='o',  c=c2, ecolor=c2, mec=c2, markersize=8)
ax3.errorbar(age[thick], mgsi[thick], yerr=mgsi_err[thick], xerr=age_err[thick], fmt='^', c=c3, ecolor=c3, mec=c2, markersize=8)
ax3.annotate(r'$\odot$', xy=(4.6, mgsi_sun), horizontalalignment='center', verticalalignment='center', color=c4, fontsize=30, weight='bold')

ax3.set_xticks(np.arange(0,11,2))
ax3.set_xticks(np.arange(0,11,1), minor=True)
ax3.set_yticks(np.arange(0.6,1.5,0.2))
ax3.set_yticks(np.arange(0.5,1.5,0.05), minor=True)
ax3.set_xlim((-0.5, 10.5))
ax3.set_ylim((0.5, 1.5))
ax3.tick_params(axis='both', which='major', labelsize=20)


ax3.set_xlabel('Age (Gyr)', fontsize=24)
ax3.set_ylabel('Mg/Si', fontsize=24)

ax4 = fig.add_subplot(224)

ax4.scatter(brewer['FeH'], brewer_mgsi, marker='o', c='black', alpha=0.1)

ax4.errorbar(par['feh'][:-1][thin], mgsi[thin], yerr=mgsi_err[thin], xerr=par['err_feh'][:-1][thin], fmt='o',  c=c2, ecolor=c2, mec=c2, markersize=8)
ax4.errorbar(par['feh'][:-1][thick], mgsi[thick], yerr=mgsi_err[thick], xerr=par['err_feh'][:-1][thick], fmt='^', c=c3, ecolor=c3, mec=c2, markersize=8)
ax4.annotate(r'$\odot$', xy=(0.0, mgsi_sun), horizontalalignment='center', verticalalignment='center', color=c4, fontsize=30, weight='bold')

ax4.set_xticks(np.arange(-0.1,0.15,0.1))
ax4.set_xticks(np.arange(-0.15,0.18,0.02), minor=True)
ax4.tick_params(axis='both', which='major', labelsize=20)
ax4.set_xlim((-0.15, 0.15))
ax4.set_ylim((0.5, 1.5))
ax4.set_yticks(np.arange(0.6,1.5,0.2))
ax4.set_yticks(np.arange(0.5,1.5,0.05), minor=True)


ax4.set_xlabel('[Fe/H]', fontsize=24)

ax4.set_yticklabels('',visible=False)

ind = a['id'] != 'sun'

co_sun = 0.55
co = 10**(a["CI_1"][ind] - a["OI_1"][ind]) * co_sun
co_err = co * np.log(10.) * np.sqrt(a["err_CI"][ind]**2 + a["err_OI"][ind]**2)
feh = par['feh'][ind]
feh_err = par['err_feh'][ind]

no_oxygen = np.where(a['id'] == 'HIP114328')[0]
co = np.delete(co, no_oxygen)
co_err = np.delete(co_err, no_oxygen)
age = np.delete(age, no_oxygen)
age_err = np.delete(age_err, no_oxygen)
feh = np.delete(feh, no_oxygen)
feh_err = np.delete(feh_err, no_oxygen)
thin = np.delete(thin, no_oxygen)
thick = np.delete(thick, no_oxygen)


brewer_co = 10**(brewer['CH'] - brewer['OH']) * co_sun

ax = fig.add_subplot(221)

ax.scatter(brewer['Age'], brewer_co, marker='o', c='black', alpha=0.1)

ax.errorbar(age[thin], co[thin], yerr=co_err[thin], xerr=age_err[thin], fmt='o',  c=c2, ecolor=c2, mec=c2, markersize=8)
ax.errorbar(age[thick], co[thick], yerr=co_err[thick], xerr=age_err[thick], fmt='^', c=c3, ecolor=c3, mec=c2, markersize=8)
ax.annotate(r'$\odot$', xy=(4.6, co_sun), horizontalalignment='center', verticalalignment='center', color=c4, fontsize=30, weight='bold')

ax.set_xticks(np.arange(0,11,2))
ax.set_xticks(np.arange(0,11,1), minor=True)
ax.set_xlim((-0.5, 10.5))
ax.set_ylim((0.2, 0.8))

ax.set_xticklabels('',visible=False)
ax.set_yticks(np.arange(0.2,0.9,0.2))
ax.set_yticks(np.arange(0.2,0.8,0.05), minor=True)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.set_ylabel('C/O', fontsize=24)

ax2 = fig.add_subplot(222)

ax2.scatter(brewer['FeH'], brewer_co, marker='o', c='black', alpha=0.1)

ax2.errorbar(feh[thin], co[thin], yerr=co_err[thin], xerr=feh_err[thin], fmt='o',  c=c2, ecolor=c2, mec=c2, markersize=8)
ax2.errorbar(feh[thick], co[thick], yerr=co_err[thick], xerr=feh_err[thick], fmt='^', c=c3, ecolor=c3, mec=c2, markersize=8)
ax2.annotate(r'$\odot$', xy=(0.0, co_sun), horizontalalignment='center', verticalalignment='center', color=c4, fontsize=30, weight='bold')

ax2.set_xticks(np.arange(-0.1,0.18,0.1))
ax2.set_xticks(np.arange(-0.18,0.18,0.02), minor=True)
ax2.set_xlim((-0.15, 0.15))
ax2.set_ylim((0.2, 0.8))

ax2.set_yticks(np.arange(0.2,0.9,0.2))
ax2.set_yticks(np.arange(0.2,0.8,0.05), minor=True)
ax2.set_xticklabels('',visible=False)
ax2.set_yticklabels('',visible=False)


fig.subplots_adjust(wspace=.05, hspace=.05)


fig.savefig('ratios.pdf')
