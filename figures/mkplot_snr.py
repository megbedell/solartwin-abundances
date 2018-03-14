import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
    
root_dir = '/Users/mbedell/Documents/Research/HARPSTwins/Results/Bulk/'
a = genfromtxt(root_dir+'summary.csv', delimiter=',', dtype=None, names=True, usecols=(0,6))
star, snr = a['star'], a['snr']

snr = np.delete(snr, np.where(star == 'HIP19911'))
snr = np.delete(snr, np.where(star == 'HIP103983'))
snr = np.delete(snr, np.where(star == 'HIP67620'))

# plot it up:
c1 = 'black'
c2 = '#003399' # blue
c3 = '#CC0033' # red
c4 = '#339900' # green


fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1E'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%i'))

bins = np.arange(200, 2500., 100.)
ax1.hist(snr, bins=bins, histtype='stepfilled', alpha=0.5, color=c2)
ax1.hist(snr, bins=bins, color=c2, lw=3, histtype='step')

ax1.set_xlim([0,2700])
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

ax1.set_xlabel(r'SNR (pix$^{-1}$)')

fig.savefig('snr.pdf')


