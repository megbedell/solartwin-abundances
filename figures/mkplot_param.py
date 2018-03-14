import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

root_dir = '/Users/mbedell/Documents/Research/HARPSTwins/Abundances/All/'
b = genfromtxt(root_dir+'final_parameters.csv', delimiter=',', dtype=None, names=True)
#ages = genfromtxt(root_dir+'final_ages_combination.csv', delimiter=',', dtype=None, names=True)

fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(16,16))
c2 = '#003399' # blue
c3 = '#CC0033' # red
c4 = '#339900' # green


bins1 = np.arange(5660.,5920.,20.)
ax1.hist(b['teff'][:-1], bins=bins1, histtype='stepfilled', alpha=0.5, color=c2)
ax1.hist(b['teff'][:-1], bins=bins1, color=c2, lw=3, histtype='step')
ax1.set_xticks(bins1[2::5])
ax1.set_xticks(bins1, minor=True)
ax1.axvline(5778., color=c3, linewidth=4, ls='dashed')
ax1.errorbar(5680, 10.5, xerr=np.median(b['err_teff']), fmt='None', c='k', ecolor='k', elinewidth=3, capsize=5, capthick=2)
ax1.set_xlabel(r'$T_{\rm eff}$ (K)')

bins2 = np.arange(4.1,4.6,0.05)
ax2.hist(b['logg'][:-1], bins=bins2, histtype='stepfilled', alpha=0.5, color=c2)
ax2.hist(b['logg'][:-1], bins=bins2, color=c2, lw=3, histtype='step')
ax2.set_xticks(bins2, minor=True)
ax2.axvline(4.44, color=c3, linewidth=4, ls='dashed')
ax2.errorbar(4.15, 16, xerr=np.median(b['err_logg']), fmt='None', c='k', ecolor='k', elinewidth=3, capsize=5, capthick=2)
ax2.set_xlabel(r'$\log g$')

bins3 = np.arange(-0.16,0.16,0.02)
ax3.hist(b['feh'][:-1], bins=bins3, histtype='stepfilled', alpha=0.5, color=c2)
ax3.hist(b['feh'][:-1], bins=bins3, color=c2, lw=3, histtype='step')
ax3.set_xticks(bins3, minor=True)
ax3.axvline(0.0, color=c3, linewidth=4, ls='dashed')
ax3.errorbar(-0.14, 8.5, xerr=np.median(b['err_feh']), fmt='None', c='k', ecolor='k', elinewidth=3, capsize=5, capthick=2)
ax3.set_xlabel(r'$\mathrm{[Fe/H]}$')

bins4 = np.arange(0.9,1.275,0.025)
ax4.hist(b['vt'][:-1], bins=bins4, histtype='stepfilled', alpha=0.5, color=c2)
ax4.hist(b['vt'][:-1], bins=bins4, color=c2, lw=3, histtype='step')
ax4.set_xticks(bins4, minor=True)
ax4.axvline(b['vt'][-1], color=c3, linewidth=4, ls='dashed')
ax4.errorbar(0.93, 13, xerr=np.median(b['err_vt']), fmt='None', c='k', ecolor='k', elinewidth=3, capsize=5, capthick=2)
ax4.set_xlabel(r'$v_t$ (km s$^{-1}$)')

fig.subplots_adjust(wspace=0.2, hspace=0.3)
fig.savefig('param_hist.pdf')

twins = (b['teff'] >= 5678.) & (b['teff'] <= 5878.) & (b['logg'] >= 4.34) & \
        (b['logg'] <= 4.54) & (b['feh'] >= -0.1) & (b['feh'] <= 0.1)
print np.sum(twins)