import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

root_dir = '../data/'
gce = genfromtxt(root_dir+'GCE/gce_linear_w_ncapture.txt', delimiter=',', dtype=None, names=True, encoding=None)

elements = ['C', 'O', 'Na', 'Mg', 'Al', 'Si', 'S', 'Ca', 'Sc',
       'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni',
       'Cu', 'Zn', 'Sr', 'Y', 'Zr', 'Ba', 'La', 'Ce', 'Pr',
       'Nd', 'Sm', 'Eu', 'Gd', 'Dy']
ms, m_errs = np.zeros(len(elements)), np.zeros((2,len(elements)))
for i,el in enumerate(elements):
    if el in ['C', 'O', 'Na', 'Mg', 'Al', 'Si', 'S', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni', 'Cu', 'Zn', 'Sr']:
        ms[i] = gce['slope'][np.where(gce['element'] == el+'I')]
        m_errs[0,i] = gce['slope_errp'][np.where(gce['element'] == el+'I')]
        m_errs[1,i] = gce['slope_errm'][np.where(gce['element'] == el+'I')]
    else:
        ms[i] = gce['slope'][np.where(gce['element'] == el+'II')]
        m_errs[0,i] = gce['slope_errp'][np.where(gce['element'] == el+'II')]
        m_errs[1,i] = gce['slope_errm'][np.where(gce['element'] == el+'II')]
        
c1 = 'black'
c2 = '#003399' # blue
c3 = '#CC0033' # red
c4 = '#339900' # green
plt.rcParams["font.sans-serif"] = "Helvetica"
plt.rcParams["axes.linewidth"] = 0.8


fig = plt.figure(figsize=(10,5), dpi=300)

ax = fig.add_subplot(111)

element_order = np.arange(len(elements))
labels = [elements[i] for i in element_order] # HACK
ax.set_xticks(element_order)
ax.set_xticklabels(labels, fontdict={'fontsize':12}, rotation=35, va="center", position=(0,-0.03))

ax.set_ylabel(r'$m$ (dex Gyr$^{-1}$)', fontsize=16)
ax.grid(True, alpha=0.3)

ax.errorbar(element_order, ms, yerr=m_errs, fmt='o', c=c2, ecolor=c2)


fig.tight_layout()
fig.savefig('gce_slopes.pdf') 

