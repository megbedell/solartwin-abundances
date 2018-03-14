import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import q2

def linear(x, m, b):
     model = m*x + b
     return model
     
root_dir = '../data/'
par = np.genfromtxt(root_dir+'final_parameters.csv', delimiter=',', dtype=None, names=True, encoding=None)

a = np.genfromtxt(root_dir+'final_abundances_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None) 
gce = np.genfromtxt(root_dir+'GCE/gce_linear_w_ncapture.txt', delimiter=',', dtype=None, names=True, encoding=None)
a_gce = np.genfromtxt(root_dir+'GCE/harpstwins_gcecorrected_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)


fit = [i not in ['HIP19911', 'HIP108158', 'HIP109821', 'HIP115577', 'HIP14501', 'HIP28066', 'HIP30476',
            'HIP33094', 'HIP65708', 'HIP73241', 'HIP74432', 'HIP64150'] for i in a['id'][:-1]] # mask out SB2, thick-disk
inv = np.invert(fit)  
print("{0} stars used; excluding {1} thick-disk or otherwise atypical stars".format(np.sum(fit), np.sum(inv)))

sp_names = gce['element']
Tc = np.asarray([q2.abundances.gettc(x) for x in sp_names])
abund = np.empty(len(Tc))
err = np.empty(len(Tc))
abund_xh = np.empty(len(Tc))
abund_gce = np.empty(len(Tc))
err_gce = np.empty(len(Tc))
n_stars = len(a['CI_1']) - 1
n_stars_gce = np.sum(fit)
for i,el in enumerate(sp_names):
    abund_xh[i] = np.log10(np.nanmean(np.power(10., a[el+"_1"][:-1])))
    abund[i] = np.log10(np.nanmean(np.power(10., a[el+"_1"][:-1] - par['feh'][:-1])))
    err[i] = np.nanstd(a[el+"_1"][:-1] - par['feh'][:-1])/np.sqrt(n_stars - 1.) # hack
    abund_gce[i] = np.log10(np.nanmean(np.power(10., a_gce[el+"Fe"][fit])))
    err_gce[i] = np.nanstd(a_gce[el+"Fe"][fit])/np.sqrt(n_stars_gce - 1.) # hack
    
    print("species {0}: median error bar {1:.4f} dex".format(el, np.median(a["err_"+el][:-1])))

    
for t in set(Tc):
   # eliminate duplicate measurements of the same element
   ind = np.where(Tc == t)[0]
   if len(ind) == 2:  # won't take care of 3+ states of the same thing
       (abund[ind[0]], err[ind[0]]) = np.average(abund[ind], weights=err[ind], returned=True)
       (abund_xh[ind[0]], _) = np.average(abund_xh[ind], weights=err[ind], returned=True)
       (abund_gce[ind[0]], err_gce[ind[0]]) = np.average(abund_gce[ind], weights=err_gce[ind], returned=True)
       sp_names[ind[0]] = '<'+sp_names[ind[0]][:-1]+'>'
       abund = np.delete(abund, ind[1])
       abund_xh = np.delete(abund_xh, ind[1])
       err = np.delete(err, ind[1])
       abund_gce = np.delete(abund_gce, ind[1])
       err_gce = np.delete(err_gce, ind[1])
       Tc = np.delete(Tc, ind[1])
       sp_names = np.delete(sp_names, ind[1])

fig = plt.figure()
ax = fig.add_subplot(111)

c2 = '#003399' # blue
c3 = '#CC0033' # red
c4 = '#339900' # green
xs = np.arange(850., 1800.)
ref = Tc >= 800.

ax.errorbar(Tc, -abund, yerr=err, fmt='o', markersize=6, c=c2, label='before GCE corrections')

popt, pcov = curve_fit(linear, Tc[ref], -abund[ref], sigma=err[ref])
ax.plot(xs, popt[0]*xs+popt[1], c=c2)

ax.errorbar(Tc, -abund_gce, yerr=err_gce, fmt='^', markersize=8, c=c4, label='after GCE corrections')

popt, pcov = curve_fit(linear, Tc[ref], -abund_gce[ref], sigma=err_gce[ref])
ax.plot(xs, popt[0]*xs+popt[1], c=c4)

names = np.asarray(['C', 'O', 'Na', 'Mg', 'Al', 'Si', 'S', 'Ca', 'Sc',
       'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni', 'Cu', 'Zn', 'Sr',
       'Y', 'Zr', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Sm',
       'Eu', 'Gd', 'Dy'])

for i, txt in enumerate(names):
    #text above
    x_txt, y_txt = Tc[i], max(-abund[i], -abund_gce[i]) + 0.02
    if txt in ['C', 'Ni', 'Y', 'Zr']:
        y_txt += 0.004
    if txt in ['Ti', 'Ce', 'Cu']:
        y_txt += 0.008
    if txt in [ 'Gd', 'Sc']:
        y_txt += 0.012
    if txt in ['Sc']:
        x_txt -= 40
    if txt in ['Mg']:
        y_txt += 0.012
        x_txt -= 30
    if txt in ['Ni']:
        x_txt -= 40
    if txt in ['Co']:
        y_txt += 0.015
        x_txt += 60
    if txt in ['Sr']:
        y_txt += 0.008
        x_txt += 40
    if txt in ['Zr', 'Y']:
        x_txt += 40
    if txt in ['Pr']:
        x_txt -= 80
        y_txt -= 0.008
    if txt in ['Sm']:
        x_txt += 200
        y_txt -= 0.015
    if txt in ['La']:
        x_txt -= 50
    if txt in ['V']:
        x_txt -= 20
    #text below
    if txt in ['Ca', 'Ce', 'Mg', 'Eu']:
        y_txt = min(-abund[i], -abund_gce[i]) - 0.02
    if txt in ['Dy']:
        y_txt = min(-abund[i], -abund_gce[i]) - 0.015
        x_txt += 50
    if txt in ['Ba']:
        y_txt = min(-abund[i], -abund_gce[i]) - 0.01
        x_txt -= 50
    if txt in ['Cr']:
        y_txt = min(-abund[i], -abund_gce[i]) - 0.02
        x_txt -= 80
    if txt in [ 'Gd']:
        y_txt = min(-abund[i], -abund_gce[i]) - 0.018
        x_txt += 80
    if txt in ['Al']:
        y_txt = min(-abund[i], -abund_gce[i]) - 0.015
        x_txt += 80
    if txt in ['Nd']:
        y_txt = min(-abund[i], -abund_gce[i]) - 0.02
        x_txt += 30
    
    x_linept, x_linetxt = Tc[i], x_txt   
    y_linept1, y_linept2, y_linetxt = -abund[i], -abund_gce[i], y_txt 
    ax.plot([x_linept, x_linetxt], [y_linept1, y_linetxt], color='k', alpha=0.2, lw=1)
    ax.plot([x_linept, x_linetxt], [y_linept2, y_linetxt], color='k', alpha=0.2, lw=1)
    ax.text(x_txt, y_txt, txt, ha='center', va='center', fontsize='small', 
            bbox={'fc':'white', 'ec':'none', 'alpha':0.8, 'pad':1})


ax.set_ylim([-0.12,0.12])
ax.set_yticks(np.arange(-0.10,0.12,0.05))
ax.set_yticks(np.arange(-0.12,0.12,0.01), minor=True)
ax.set_xticks(np.arange(0,1800,250))
ax.set_xticks(np.arange(0,1850,50), minor=True)

ax.legend(loc='upper right', frameon=True)
ax.set_xlabel(r'$T_\mathrm{C}$ (K)')
plt.ylabel(r'Sun - [X/Fe]$_{<twin>}$')

fig.savefig('avgtwin_w_ncapture.pdf')
