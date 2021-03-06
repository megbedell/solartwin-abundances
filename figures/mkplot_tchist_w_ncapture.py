import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
from matplotlib.ticker import FormatStrFormatter

def linear(x, m, b):
     model = m*x + b
     return model
     
def percentile_uncertainty(value, xs, sigmas, Ntrial=int(2e4)):
    N = len(xs)
    assert len(sigmas) == N
    percentiles = np.zeros(Ntrial)
    for trial in range(Ntrial):
        inds = np.random.randint(0, N, size=N)
        thisxs = xs[inds] + sigmas[inds] * np.random.normal(size=N)
        percentiles[trial] = np.sum(thisxs <= value) * 100. / N
    return percentiles
     
root_dir = '../data/'
a = genfromtxt(root_dir+'final_abundances_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)
par = genfromtxt(root_dir+"final_parameters.csv", delimiter=',', dtype=None, names=True, encoding=None)
gce = genfromtxt(root_dir+'GCE/gce_linear_w_ncapture.txt', delimiter=',', dtype=None, names=True, encoding=None)
a_gce = genfromtxt(root_dir+'GCE/harpstwins_gcecorrected_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)
ages = genfromtxt(root_dir+'final_ages_combination.csv', delimiter=',', dtype=None, names=True, encoding=None)

tc_all = genfromtxt(root_dir+'Tc/harpstwins_tcslopes_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)
tc_all_pw = genfromtxt(root_dir+'Tc/harpstwins_tcslopes_w_ncapture_pw.csv', delimiter=',', dtype=None, names=True, encoding=None)


if True:
    fit = [i not in ['HIP19911', 'HIP108158', 'HIP109821', 'HIP115577', 'HIP14501', 'HIP28066', 'HIP30476',
                    'HIP33094', 'HIP65708', 'HIP73241', 'HIP74432', 'HIP64150'] for i in a['id'][:-1]] # mask out SB2, thick-disk
    inv = np.invert(fit)  
    print("{0} stars used; excluding {1} thick-disk or otherwise atypical stars".format(np.sum(fit), np.sum(inv)))

tc_slopes_nov = tc_all['tc_slope_nov']
tc_slope_errs_nov = tc_all['tc_slope_err_nov']
tc_slopes_gce_nov = tc_all['tc_slope_gce_nov'][fit]
tc_slope_errs_gce_nov = tc_all['tc_slope_err_gce_nov'][fit]
sun_age = (ages['age_mean'] > 3.5) & (ages['age_mean'] < 5.5)
sun_age_gce = (ages[fit]['age_mean'] > 3.5) & (ages[fit]['age_mean'] < 5.5)

print("The stars at solar age are:")
print(a['id'][:-1][sun_age])


if False: # use piecewise instead
    tc_slopes_nov = tc_all_pw['tc_slope']
    tc_slope_errs_nov = tc_all_pw['tc_slope_err']
    tc_slopes_gce_nov = tc_all_pw['tc_slope_gce'][fit]
    tc_slope_errs_gce_nov = tc_all_pw['tc_slope_err_gce'][fit]

percentiles = percentile_uncertainty(0.0, tc_slopes_gce_nov, tc_slope_errs_gce_nov)


print("{0} stars of {1} total have Tc slopes at or below the Solar value after GCE correction.".format(np.sum(tc_slopes_nov <= 0.0), len(tc_slopes_nov)))
print("The Sun lies at the {0:.1f} percentile before GCE corrections".format(np.sum(tc_slopes_nov <= 0.0)*100. / len(tc_slopes_nov)))
print("{0} stars of {1} total have Tc slopes at or below the Solar value after GCE correction.".format(np.sum(tc_slopes_gce_nov <= 0.0), len(tc_slopes_gce_nov)))
print("The Sun lies at the {0:.1f} percentile after GCE corrections".format(np.sum(tc_slopes_gce_nov <= 0.0)*100. / len(tc_slopes_gce_nov)))
print("The uncertainty on this is: {0:.1f} +- {1:.1f}".format(np.mean(percentiles), np.std(percentiles)))
print("At 95% confidence, the Sun is below the {0:.1f}th percentile.".format(np.percentile(percentiles, 95)))

c1 = 'black'
c2 = '#003399' # blue
c3 = '#CC0033' # red
c4 = '#339900' # green
plt.rcParams["font.sans-serif"] = "Helvetica"


fig = plt.figure(figsize=(8,12))

ax1 = fig.add_subplot(211)

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1E'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%i'))

ax1.hist(tc_slopes_nov, bins=np.arange(-8.e-5, 40.e-5, 4.e-5),histtype='stepfilled', alpha=0.4, color=c2, label='before GCE corrections')
ax1.hist(tc_slopes_nov, bins=np.arange(-8.e-5, 40.e-5, 4.e-5), color=c2, lw=2, histtype='step')
#ax1.errorbar(-5.e-5, 15., xerr=np.median(tc_slope_errs_nov), fmt='None', c='k', ecolor='k', elinewidth=3, capsize=5, capthick=2)

ax1.hist(tc_slopes_gce_nov, bins=np.arange(-8.e-5, 40.e-5, 4.e-5),histtype='stepfilled', alpha=0.4, color=c4, label='after GCE corrections')
ax1.hist(tc_slopes_gce_nov, bins=np.arange(-8.e-5, 40.e-5, 4.e-5), color=c4, lw=2, histtype='step')

ax1.tick_params(axis='x', labelsize=14)
ax1.set_yticks(np.arange(0,19,5))
ax1.set_yticks(np.arange(0,20,1), minor=True)
#ax1.set_xticks(np.arange(-8.e-5,40.e-5,4.e-5), minor=True)
 
ax1.set_ylabel(r'Number of stars (all)', fontsize=18)
ax1.legend(loc='upper right', frameon=True, fontsize=18)

ax1.axvline(0.0, lw=2, color=c1, ls='dotted')
ax1.text(0.6e-5, 8., 'Sun', color=c1, size=18)

ax2 = fig.add_subplot(212)

ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1E'))

ax2.hist(tc_slopes_nov[sun_age], bins=np.arange(-8.e-5, 40.e-5, 4.e-5),histtype='stepfilled', alpha=0.4, color=c2, label='before GCE corrections')
ax2.hist(tc_slopes_nov[sun_age], bins=np.arange(-8.e-5, 40.e-5, 4.e-5), color=c2, lw=2, histtype='step')
#ax1.errorbar(-5.e-5, 15., xerr=np.median(tc_slope_errs_nov), fmt='None', c='k', ecolor='k', elinewidth=3, capsize=5, capthick=2)

ax2.hist(tc_slopes_gce_nov[sun_age_gce], bins=np.arange(-8.e-5, 40.e-5, 4.e-5),histtype='stepfilled', alpha=0.4, color=c4, label='after GCE corrections')
ax2.hist(tc_slopes_gce_nov[sun_age_gce], bins=np.arange(-8.e-5, 40.e-5, 4.e-5), color=c4, lw=2, histtype='step')

ax2.tick_params(axis='x', labelsize=14)
ax2.set_yticks(np.arange(0,9,5))
ax2.set_yticks(np.arange(0,9,1), minor=True)
#ax1.set_xticks(np.arange(-8.e-5,40.e-5,4.e-5), minor=True)
 
ax2.set_xlabel(r'$T_{c}$ slope (dex K$^{-1}$)', fontsize=18)
ax2.set_ylabel(r'Number of stars (Solar age)', fontsize=18)

ax2.legend(loc='upper right', frameon=True, fontsize=18)

ax2.axvline(0.0, lw=2, color=c1, ls='dotted')
ax2.text(0.6e-5, 3., 'Sun', color=c1, size=18)


fig.tight_layout()
fig.savefig('tc_histogram_w_ncapture.pdf') 
fig.clear()
