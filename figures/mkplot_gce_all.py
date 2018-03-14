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
a = genfromtxt(root_dir+'final_abundances_w_ncapture.csv', delimiter=',', dtype=None, names=True)
#data = q2.Data(root_dir+"final_parameters.csv",root_dir+"harpstwins_lines.csv")
par = np.genfromtxt(root_dir+"final_parameters.csv", delimiter=',', dtype=None, names=True)
gce = genfromtxt(root_dir+'GCE/gce_linear_w_ncapture.txt', delimiter=',', dtype=None, names=True)
ages = np.genfromtxt(root_dir+'final_ages_combination.csv', delimiter=',', dtype=None, names=True)

age = ages['age_mean']
age_err = ages['age_std']
xs = np.arange(11.)
elements = ['C', 'O', 'Na', 'Mg', 'Al', 'Si', 'S', 'Ca', 'Sc',
       'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni',
       'Cu', 'Zn', 'Sr', 'Y', 'Zr', 'Ba', 'La', 'Ce', 'Pr',
       'Nd', 'Sm', 'Eu', 'Gd', 'Dy']
ms, bs = np.zeros(len(elements)), np.zeros(len(elements))
for i,el in enumerate(elements):
    if el in ['C', 'O', 'Na', 'Mg', 'Al', 'Si', 'S', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni', 'Cu', 'Zn', 'Sr']:
        ms[i] = gce['slope'][np.where(gce['element'] == el+'I')]
        bs[i] = gce['intercept'][np.where(gce['element'] == el+'I')]
    else:
        ms[i] = gce['slope'][np.where(gce['element'] == el+'II')]
        bs[i] = gce['intercept'][np.where(gce['element'] == el+'II')]

fit = [i not in ['HIP19911', 'HIP108158', 'HIP109821', 'HIP115577', 'HIP14501', 'HIP28066', 'HIP30476',
        'HIP33094', 'HIP65708', 'HIP73241', 'HIP74432', 'HIP64150'] for i in a['id'][:-1]] # mask out SB2, thick-disk
inv = np.invert(fit)
print("{0} stars used in fit; {1} stars excluded.".format(np.sum(fit), np.sum(inv)))

c2 = '#003399' # blue
c3 = '#CC0033' # red
c4 = '#339900' # green
plt.rcParams["font.sans-serif"] = "Helvetica"
plt.rcParams["axes.linewidth"] = 0.6

fig = plt.figure(figsize=(6,6), dpi=300)

for i,el in enumerate(elements):
    ax = fig.add_subplot(6,5,i+1)
    
    if el in ['C', 'O', 'Na', 'Mg', 'Al', 'Si', 'S', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni', 'Cu', 'Zn', 'Sr']:
        abund = a[el+"I_1"][:-1] - par['feh'][:-1] # exclude sun
        err = a["err_"+el+'I'][:-1]
    else:
        abund = a[el+"II_1"][:-1] - par['feh'][:-1] # exclude sun
        err = a["err_"+el+'II'][:-1]

    ax.errorbar(age[inv], abund[inv], xerr=age_err[inv], yerr=err[inv], fmt='D', c=c3, ecolor=c3, ms=2, elinewidth=0.6, alpha=0.8)
    ax.errorbar(age[fit], abund[fit], xerr=age_err[fit], yerr=err[fit], fmt='o', c='black', ecolor='black', mec='black', ms=2, elinewidth=0.6)
    ax.annotate(r'$\odot$', xy=(4.6, -0.01), horizontalalignment='center', verticalalignment='center', color=c4, fontsize=8, weight='bold')

    ax.plot(xs, ms[i]*xs + bs[i], color=c2, lw=1.0)
    
    ax.tick_params(length=3, width=0.6, labelsize=5, which='major', pad=4)
    ax.tick_params(length=2, width=0.4, which='minor')    
    ax.set_xticks(np.arange(0,11,2))
    ax.set_xticks(np.arange(0,11,1), minor=True)
    

    if i<20:
        ax.set_ylim([-0.2,0.25])
        ax.set_yticks(np.arange(-0.2,0.3,0.1))
        ax.set_yticks(np.arange(-0.2,0.3,0.05), minor=True)
        ax.text(1.0,0.15, el, size=8)
    else:
        ax.set_ylim([-0.15,0.45])
        ax.set_yticks(np.arange(-0.1,0.5,0.1))
        ax.set_yticks(np.arange(-0.15,0.5,0.05), minor=True)
        ax.text(1.0,0.3, el, size=8)
    
    if el in ['C', 'Cr', 'Sc', 'Ti']:
        if el == 'C':
            el2 = 'CH'
        else:
            el2 = el+'II'
        abund = a[el2+"_1"][:-1] - par['feh'][:-1] # exclude sun
        err = a["err_"+el2][:-1]
        ax.errorbar(age[inv], abund[inv], xerr=age_err[inv], yerr=err[inv], ls='None', marker=(4,0,0), mfc='None', 
                    ecolor=c3, mec=c3, mew=0.3, ms=2, elinewidth=0.4)
        ax.errorbar(age[fit], abund[fit], xerr=age_err[fit], yerr=err[fit], ls='None', marker=(1,3,0), mfc='None', 
                    ecolor='black', mec='black', mew=0.6, ms=2.2, elinewidth=0.6)
        

    if (i % 5) != 0:
        ax.set_yticklabels('',visible=False)

    if el not in elements[-5:]:
        ax.set_xticklabels('',visible=False)

fig.subplots_adjust(hspace=.05, wspace=.05)
fig.text(0.5, 0.07, 'Age (Gyr)', size=10, ha='center')
fig.text(0.06, 0.5, '[X/Fe] (dex)', rotation=90, size=10, va='center')
fig.savefig('gce_all.pdf')
