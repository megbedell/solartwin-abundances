import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    root_dir = '../'
    a = np.genfromtxt(root_dir+'final_abundances_w_ncapture.csv', delimiter=',', dtype=None, names=True, encoding=None)
    par = np.genfromtxt(root_dir+"final_parameters.csv", delimiter=',', dtype=None, names=True, encoding=None)
    ages = np.genfromtxt(root_dir+'final_ages_combination.csv', delimiter=',', dtype=None, names=True, encoding=None)
    
    age = ages['age_mean']
    age_err = ages['age_std']
    
    gce = np.genfromtxt('gce_linear_w_ncapture.txt', delimiter=',', dtype=None, names=True, encoding=None)
    
    abund_all = np.zeros((len(age), len(gce['element'])*2+1))
    #abund_all[:,0] = ages['id']
    names_all = ['id']
    
    for i,el in enumerate(gce['element']):
        abund = a["{0}_1".format(el)][:-1] - par['feh'][:-1] # exclude sun
        err = a["err_{0}".format(el)][:-1]
        
        slope = gce['slope'][i]
        slope_errp = gce['slope_errp'][i]
        slope_errm = gce['slope_errm'][i]
        gce_errp = np.abs((age - 4.5) * (slope+slope_errp) - (age - 4.5) * slope)
        gce_errm = np.abs((age - 4.5) * (slope-slope_errm) - (age - 4.5) * slope)
        gce_err = np.mean([gce_errp, gce_errm], axis=0)  # really I should save abund_errp and abund_errm separately
        
        # apply corrections:
        
        abund -= (age - 4.5) * slope # differential correction
        err = np.sqrt(err**2 + gce_err**2)
        
        # save it:
        
        abund_all[:,2*i+1] = abund        
        names_all = np.append(names_all, '[{0}/Fe]'.format(el))
        abund_all[:,2*i+2] = err
        names_all = np.append(names_all, 'err_{0}'.format(el))
    
    header = np.array2string(names_all, separator=',', max_line_width=5000)    
    header = str.replace(header, '\'', '')
    header = str.lstrip(header, '[')
    header = str.rstrip(header, ']')
    np.savetxt('harpstwins_gcecorrected_w_ncapture.csv', abund_all, fmt='%s', delimiter=',', header=header)
    
    # NOTE: this does not save perfectly! must go in and reformat header + paste in star IDs