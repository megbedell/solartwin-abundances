import numpy as np

# equivalent width table headers

root_dir = '/Users/mbedell/Documents/Research/HARPSTwins/Abundances/All/'
a = np.genfromtxt(root_dir+'final_abundances_w_ncapture.csv', delimiter=',', dtype=None, names=True)
gce = np.genfromtxt(root_dir+'GCE/gce_linear_w_ncapture.txt', delimiter=',', dtype=None, names=True)
linelist = np.genfromtxt(root_dir+'final_lines.csv', delimiter=',', dtype=None, names=True)
f = open('../aas-journals/ews_table_header.txt', 'w')
f_lbls = open('mrt_ews_labels.txt', 'w')
f_units = open('mrt_ews_units.txt', 'w')
f_exps = open('mrt_ews_explanations.txt', 'w')

f.write(r'1 & Wavelength & \r{A} & Rest wavelength of line \\')
f.write('\n')
f.write(r'2 & Species &  & Species identifier \\')
f.write('\n')
f.write(r'3 & EP & eV & Excitation potential \\')
f.write('\n')
f.write(r'4 & log($gf$) & dex & Log of the oscillator strength \\')
f.write('\n')

f_lbls.write('Wavelength\n')
f_lbls.write('Species\n')
f_lbls.write('EP\n')
f_lbls.write('log(gf)\n')

f_units.write('0.1nm\n')
f_units.write('---\n')
f_units.write('eV\n')
f_units.write('[-]\n')

f_exps.write('Rest wavelength of line in Angstroms\n')
f_exps.write('Species identifier\n')
f_exps.write('Excitation potential\n')
f_exps.write('Log of the oscillator strength\n')


for i,s in enumerate(a['id']):
    if s == 'sun':
        name = s
    else:
        hipno = s[3:]
        name = '\hip{'+hipno+'}'
    f.write('{0} & {1} & m\\r{{A}} & \\\ \n'.format(i+5, name))
    f_lbls.write('{0}\n'.format(s))
    f_units.write('10-13m\n')
    no_nans_check = np.isfinite(linelist[s]).all()
    if not no_nans_check:
        f_exps.write('?')
    f_exps.write('Measured equivalent width for {0}\n'.format(s))

f.close()
f_lbls.close()
f_units.close()
f_exps.close()

# equivalent widths table contents
#f = open('mrt_ews.txt', 'w')

#f.close()


# abundances table headers

f = open('../aas-journals/abund_table_header.txt', 'w')
f_lbls = open('mrt_abunds_labels.txt', 'w')
f_units = open('mrt_abunds_units.txt', 'w')
f_exps = open('mrt_abunds_explanations.txt', 'w')

f.write(r'1 & Star & & Star identifier \\')
f.write('\n')

f_lbls.write('Star\n')
f_units.write('---\n')
f_exps.write('Star identifier\n')

for i,sp in enumerate(gce['element'][:21]):
    name = sp.replace('I', ' \I', 1)
    f.write('{0} & [{1} /H] & dex & \\\ \n'.format(2 + i*2, name))
    f.write('{0} & u\\_[{1} /H] & dex & Estimated uncertainty \\\ \n'.format(3 + i*2, sp))
    f_lbls.write('[{0}/H]\n'.format(sp))
    f_lbls.write('u_[{0}/H]\n'.format(sp))
    f_units.write('[-]\n')
    f_units.write('[-]\n')
    if sp=='OI':
        f_exps.write('?')
    f_exps.write('Differential abundance of species {0} relative to Solar\n'.format(sp))
    if sp=='OI':
        f_exps.write('?')
    f_exps.write('Estimated uncertainty\n')
    
f.close()
f_lbls.close()
f_units.close()
f_exps.close()
