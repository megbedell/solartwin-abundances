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
         
def fit_slope(abund, err, Tc, include_volatiles=True, std=False, plot=False, starname='star', func=broken_line):    
    for t in np.unique(Tc):
       # eliminate duplicate measurements of the same element
       ind = np.where(Tc == t)[0]
       if len(ind) == 2:  # won't take care of 3+ states of the same thing
           (abund[ind[0]], err[ind[0]]) = np.average(abund[ind], weights=err[ind], returned=True)       
           abund = np.delete(abund, ind[1])
           err = np.delete(err, ind[1])
           Tc = np.delete(Tc, ind[1])
                  
    if (include_volatiles):
        good = (Tc > 0.0) & np.isfinite(abund)
        xs = np.arange(0., 1800., 20) # plotting
    else:
        good = (Tc > 900.0) & np.isfinite(abund)
        xs = np.arange(900., 1800., 20) # plotting

    if func==broken_line:
        p0 = [5.e-4, 0., -0.3, 0.2] # x_break = 1000 K
    else:
        p0 = [5.e-4, -0.3]
    popt, pcov = curve_fit(func, Tc[good], abund[good], sigma=err[good], p0=p0)
    perr = np.sqrt(np.diag(pcov))
        
    if plot:
        plt.errorbar(Tc, abund, err, fmt='o')
        plt.plot(xs, func(xs, *popt))
        plt.xlabel(r'$T_\mathrm{C}$')
        plt.ylabel('[X/H]')
        plt.title(starname)
        plt.savefig('plots/'+starname+'_Tcfit_'+func.__name__+'.png')
        plt.clf()
    if std:
        stdev = np.std(abund[good] - func(Tc[good], *popt))
        return popt[0], perr[0], stdev
    else:
        return popt[0], perr[0]
    
    
 
if __name__ == "__main__":
    root_dir = '../'
    a = np.genfromtxt(root_dir+'final_abundances_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)
    gce = np.genfromtxt(root_dir+'GCE/gce_linear_w_ncapture.txt', delimiter=',', dtype=None, names=True, encoding=None)
    a_gce = np.genfromtxt(root_dir+'GCE/harpstwins_gcecorrected_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)
    par = np.genfromtxt(root_dir+"final_parameters.csv", delimiter=',', dtype=None, names=True, encoding=None)
    
    f = open('harpstwins_tcslopes_w_ncapture.csv', 'w')
    f.write('id, tc_slope, tc_slope_err, tc_slope_nov, tc_slope_err_nov, tc_std_nov, tc_slope_gce, tc_slope_err_gce, tc_slope_gce_nov, tc_slope_err_gce_nov, tc_std_gce_nov \n')

    f2 = open('../../aas-journals/tc_table.txt', 'w')
    
    f3 = open('harpstwins_tcslopes_w_ncapture_pw.csv', 'w')
    f3.write('id, tc_slope, tc_slope_err, tc_slope_gce, tc_slope_err_gce \n')
    
    for i,s in enumerate(a['id']):
        print(s)
        if (s == 'sun'):
            continue
        assert a['id'][i] == par['id'][i]
        assert a['id'][i] == a_gce['id'][i]
        
        sp_names = gce['element']
        Tc = np.asarray([q2.abundances.gettc(x) for x in sp_names])
        abund = np.asarray([a[sp+"_1"][i] for sp in sp_names])
        err = np.asarray([a["err_"+sp][i] for sp in sp_names])
        abund_gce = np.asarray([a_gce[sp+"Fe"][i] + par['feh'][i] for sp in sp_names])
        err_gce = np.asarray([np.sqrt(a_gce["err_"+sp][i]**2 + par['err_feh'][i]**2) for sp in sp_names])
                
        # append iron:
        sp_names = np.append(sp_names, 'Fe')
        Tc = np.append(Tc, q2.abundances.gettc('FeI'))
        abund = np.append(abund, par['feh'][i])
        err = np.append(err, par['err_feh_'][i])
        abund_gce = np.append(abund_gce, par['feh'][i])
        err_gce = np.append(err_gce, par['err_feh_'][i])
        
            
        # fit Tc slopes without GCE corrections:
        tc_slope, tc_slope_err = fit_slope(abund, err, Tc, include_volatiles = True, func=linear)
        tc_slope_nov, tc_slope_err_nov, tc_std_nov = fit_slope(abund, err, Tc, include_volatiles = False, std=True, func=linear)
        
        # fit Tc slopes with GCE corrections:
        tc_slope_gce, tc_slope_err_gce = fit_slope(abund_gce, err_gce, Tc, include_volatiles = True, func=linear)
        tc_slope_gce_nov, tc_slope_err_gce_nov, tc_std_gce_nov = fit_slope(abund_gce, err_gce, Tc, include_volatiles = False, std=True, plot=True, starname=s, func=linear)


        f.write('{id}, {s1:8e}, {ds1:8e}, {s2:8e}, {ds2:8e}, {std2:8e}, {s3:8e}, {ds3:8e}, {s4:8e}, {ds4:8e}, {std4:8e} \n'.format(id=s, 
                s1=tc_slope, ds1=tc_slope_err, s2=tc_slope_nov, ds2=tc_slope_err_nov, s3=tc_slope_gce, ds3=tc_slope_err_gce,
                s4=tc_slope_gce_nov, ds4=tc_slope_err_gce_nov, std2=tc_std_nov, std4=tc_std_gce_nov))
                
        tc_slope_pw, tc_slope_err_pw = fit_slope(abund, err, Tc, include_volatiles = True, func=broken_line)
        tc_slope_gce_pw, tc_slope_err_gce_pw = fit_slope(abund_gce, err_gce, Tc, include_volatiles = True, func=broken_line)
        
        f3.write('{id}, {s1:8e}, {ds1:8e}, {s3:8e}, {ds3:8e} \n'.format(id=s, 
                s1=tc_slope_pw, ds1=tc_slope_err_pw, s3=tc_slope_gce_pw, ds3=tc_slope_err_gce_pw))
        
        hipno = s[3:]
        f2.write('HIP {id} & {s1:.2f} & {s2:.2f} & {s3:.2f} & {s4:.2f} \\\ \n'.format(id=hipno, 
                s1=tc_slope_nov*1.e4, s2=tc_slope_gce_nov*1.e4,
                s3=tc_slope_pw*1.e4, s4=tc_slope_gce_pw*1.e4))
                
    f.close()
    f2.close()
    f3.close()
