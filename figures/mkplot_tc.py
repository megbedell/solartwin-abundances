import numpy as np
import matplotlib.pyplot as plt
import q2
from scipy.optimize import curve_fit


def linear(x, m, b):
     model = m*x + b
     return model
     
def broken_line(x, m2, m1, b2, b1):
    x_break = (b2 - b1) / (m1 - m2)
    return np.where(x < x_break, m1*x + b1, m2*x + b2)

root_dir = '../data/'
a = np.genfromtxt(root_dir+'final_abundances_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)
par = np.genfromtxt(root_dir+"final_parameters.csv", delimiter=',', dtype=None, names=True, encoding=None)
gce = np.genfromtxt(root_dir+'GCE/gce_linear_w_ncapture.txt', delimiter=',', dtype=None, names=True, encoding=None)
a_gce = np.genfromtxt(root_dir+'GCE/harpstwins_gcecorrected_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)

tc_all = np.genfromtxt(root_dir+'Tc/harpstwins_tcslopes_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)

c1 = 'black'
c2 = '#003399' # blue
c3 = '#CC0033' # red
c4 = '#339900' # green
plt.rcParams["font.sans-serif"] = "Helvetica"

stars = ['HIP54287', 'HIP74389', 'HIP4909']
slopes = []
slope_errs = []


fig = plt.figure(figsize=(6,6))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

for (s,ax) in zip(stars, [ax1, ax2, ax3]):
    i = np.where(a['id'] == s)[0]
    
    sp_names = gce['element']
    Tc = np.asarray([q2.abundances.gettc(x) for x in sp_names])
    abund = np.asarray([a[sp+"_1"][i] for sp in sp_names])
    err = np.asarray([a["err_"+sp][i] for sp in sp_names])
    abund_gce = np.asarray([a_gce[sp+"Fe"][i] + par['feh'][i] for sp in sp_names])
    err_gce = np.asarray([np.sqrt(a_gce["err_"+sp][i]**2 + par['err_feh'][i]**2) for sp in sp_names])
    
    sp_names = np.append(sp_names, 'Fe')
    Tc = np.append(Tc, q2.abundances.gettc('FeI'))
    abund = np.append(abund, par['feh'][i])
    err = np.append(err, par['err_feh_'][i])
    abund_gce = np.append(abund_gce, par['feh'][i])
    err_gce = np.append(err_gce, par['err_feh_'][i])
    
    for t in np.unique(Tc):
       # eliminate duplicate measurements of the same element
       ind = np.where(Tc == t)[0]
       if len(ind) == 2:  # won't take care of 3+ states of the same thing
           (abund[ind[0]], err[ind[0]]) = np.average(abund[ind], weights=err[ind], returned=True)       
           abund = np.delete(abund, ind[1])
           err = np.delete(err, ind[1])
           Tc = np.delete(Tc, ind[1])
           
    good = (Tc > 900.0) & np.isfinite(abund)
    xs = np.arange(900., 1800., 20) # plotting
    
    p0 = [5.e-4, -0.3]
    popt, pcov = curve_fit(linear, Tc[good], abund[good], sigma=err[good], p0=p0)
    perr = np.sqrt(np.diag(pcov))
    
    p0 = [5.e-4, 0., -0.3, 0.2] # x_break = 1000 K
    popt2, pcov2 = curve_fit(broken_line, Tc, abund, sigma=err, p0=p0)
    xs2 = np.arange(100., 1800., 10) # plotting
    
    ax.errorbar(Tc, abund, err, fmt='o', color='k', mec='k', markersize=5, lw=2)
    ax.plot(xs, linear(xs, *popt), color=c2, lw=2)
    ax.plot(xs2, broken_line(xs2, *popt2), color=c3, ls='--', lw=2)
    
    
    ax.set_xticks(np.arange(0,1800,500))
    ax.set_xticks(np.arange(0,1850,100), minor=True)
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    ax.set_yticks(np.arange(-0.1,0.25,0.1))
    ax.set_yticks(np.arange(-0.15,0.3,0.05), minor=True)
    ax.set_ylim([-0.17, 0.27])
    
    name = 'HIP '+s[3:]
    ax.text(0.05, 0.8, name, transform=ax.transAxes, fontsize=16)
    
    #slopes = np.append(slopes, popt[0])
    #slope_errs = np.append(slope_errs, perr[0])
    #ax.text(0.05, 0.9, "{0} \n $({1:.2f} \pm {2:.2f}) \\times 10^{{-4}}\ dex\ K^{{-1}}$".format(s, 
    #            popt[0]*1.e4, perr[1]*1.e4), transform=ax.transAxes, fontsize=16)
    

ax1.set_xticklabels('',visible=False)    
ax2.set_xticklabels('',visible=False) 
ax2.set_ylabel('[X/H] (dex)', fontsize=18)
ax3.set_xlabel(r'$T_\mathrm{C}$ (K)', fontsize=18)   
fig.subplots_adjust(hspace=.05) 
    
plt.savefig('tc.pdf')
