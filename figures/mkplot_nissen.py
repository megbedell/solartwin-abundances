import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from starchive import identifiers

def linear(x, m, b):
     model = m*x + b
     return model
     
root_dir = '/Users/mbedell/Documents/Research/HARPSTwins/Abundances/All/'
a = np.genfromtxt(root_dir+'final_abundances.csv', delimiter=',', dtype=None, names=True)
par = np.genfromtxt(root_dir+"final_parameters.csv", delimiter=',', dtype=None, names=True)
nissen_oxygen = np.genfromtxt("nissen_oxygen.txt", dtype=None, names=True)

nissen = np.genfromtxt("nissen.txt", delimiter='\t', dtype=None, names=True)
nissen_mask = []
mask = []
conv = identifiers.Converter()
for i,hdname in enumerate(nissen['Star']):
    hdnumber = float(hdname[3:])
    hipnumber = conv.hdtohip(hdnumber)
    hipname = 'HIP' + str(hipnumber)
    if hipname in a['id']:
        print '{0} == {1}'.format(hdname, hipname)
        nissen_mask = np.append(nissen_mask, i)
        mask = np.append(mask, np.where(a['id'] == hipname)[0][0])
        
nissen_mask = [int(i) for i in nissen_mask] # idk why this is necessary but it is
mask = [int(i) for i in mask]


c2 = '#003399' # blue
c3 = '#CC0033' # red
c4 = '#339900' # green
plt.rcParams["font.sans-serif"] = "Helvetica"

fig = plt.figure()

xs = np.arange(-0.25,0.3,0.1)
elements = ['C', 'O', 'Si', 'Ni']
nissen_errors = [0.013, 0.0, 0.007, 0.006]

for i,el in enumerate(elements):
    ax = fig.add_subplot(2,2,i+1)

    abund = a[el+"I_1"][mask]
    err = a["err_"+el+"I"][mask]
    if el == 'C':
        abunds = [a['CI_1'][mask], a['CH_1'][mask]]
        errs = [a['err_CI'][mask], a['err_CH'][mask]]
        (abund, err) = np.average(abunds, weights=errs, returned=True, axis=0) 
    
    abund_nissen = nissen[el+'Fe'][nissen_mask] + nissen['FeH'][nissen_mask]
   
    if el == 'O':
        err_nissen = nissen_oxygen['sigma'][nissen_mask]
    else:
        err_nissen = np.zeros_like(abund_nissen) + nissen_errors[i]
    

    ax.errorbar(abund_nissen, abund, xerr=err_nissen, yerr=err, fmt='o', c='black', ecolor='black', mec='black', ms=7)

    ax.plot(xs, xs, color=c2, lw=2, ls='--')
    
    diff = abund - abund_nissen
    ax.text(-0.19, 0.12, '{0}: $\mu$ = {1:.3f} dex\n  $\sigma$ = {2:.3f} dex'.format(el, np.mean(diff), np.std(diff)), fontsize=18)

    ax.set_ylim([-0.22,0.22])
    ax.set_xlim([-0.22,0.22])

    ax.set_yticks(np.arange(-0.2,0.22,0.1))
    ax.set_yticks(np.arange(-0.2,0.22,0.05), minor=True)
    
    ax.set_xticks(np.arange(-0.2,0.22,0.1))
    ax.set_xticks(np.arange(-0.2,0.22,0.05), minor=True)


    ax.tick_params(axis='both', which='major', labelsize=16)

    if (i % 2) != 0:
        ax.set_yticklabels('',visible=False)

    if el not in elements[-2:]:
        ax.set_xticklabels('',visible=False)

fig.subplots_adjust(hspace=.05, wspace=.05)
fig.text(0.5, 0.015, '[X/H] (Nissen)', size=20, ha='center')
fig.text(0.015, 0.5, '[X/H] (this work)', rotation=90, size=20, va='center')
fig.savefig('nissen.pdf')
