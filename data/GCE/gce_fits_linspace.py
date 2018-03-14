import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import emcee
import corner

def vee(par):
    # Hogg+2010 eqn 29
    (m, b, lnjitter) = par
    return 1./np.sqrt(1. + m**2) * np.asarray([-m, 1.])
    
def ortho_displacement(par, ys, xs):
    # Hogg+2010 eqn 30
    (m, b, lnjitter) = par
    disp = np.zeros_like(ys)
    for i, (y, x) in enumerate(zip(ys, xs)):
        z0 = np.asarray([0.0, b])
        zi = np.asarray([x, y])
        disp[i] = np.dot( vee(par), zi - z0 )
    return disp
    
def ortho_variance(par, dys, dxs):
    # Hogg+2010 eqn 31
    #(m, b, jitter) = par
    var = np.zeros_like(dys)
    for i, (dy, dx) in enumerate(zip(dys, dxs)):
        cov = np.eye(2)
        cov[0,0] = dx**2
        cov[1,1] = dy**2 #+ jitter**2
        var[i] = np.dot( np.dot(vee(par), cov), vee(par) )
    return var

def twodlike(par, y, dy, x, dx):
    # log(likelihood) considering errors in both x and y
    # Hogg+2010 eqn 31 with jitter
    (m, b, lnjitter) = par
    delta = ortho_displacement(par, y, x)
    sigmasq = ortho_variance(par, dy, dx) + np.exp(2.*lnjitter)
    return -0.5 * np.sum(delta**2/sigmasq + np.log(sigmasq))
    
def lnprior(par):
    (m, b, lnjitter) = par
    if (-5. < m < 5.) and (-1. < b < 1.) and (-20. < lnjitter < 1.):
        return 0.0
    return -np.inf

def lnprob_1sp(par, y, dy, x, dx):
    lp = lnprior(par)
    if not np.isfinite(lp):
        return -np.inf
    return lp + twodlike(par, y, dy, x, dx)

def bestfit_1sp(abund, err, age, age_err):   
    # fit
    nll = lambda *args: -lnprob_1sp(*args)
    par0 = np.asarray([0.0, 0.0, 0.0])
    result = op.minimize(nll, par0, args=(abund, err, age, age_err))
    (m, b, lnjitter) = result['x']    
    return m, b, lnjitter
    
    
def lnprob_2sp(par, y1, dy1, y2, dy2, x, dx):
    par1 = par[0:3]
    par2 = np.append(par[0], par[3:])
    lp1 = lnprior(par1)
    lp2 = lnprior(par2)
    if not (np.isfinite(lp1) and np.isfinite(lp1)):
        return -np.inf
    return lp1 + lp2 + twodlike(par1, y1, dy1, x, dx) + twodlike(par2, y2, dy2, x, dx)

def bestfit_2sp(abund1, err1, abund2, err2, age, age_err):   
    # fit
    nll = lambda *args: -lnprob_2sp(*args)
    par0 = np.asarray([0.0, 0.0, 0.0, 0.0, 0.0])
    result = op.minimize(nll, par0, args=(abund1, err1, abund2, err2, age, age_err))
    (m, b1, lnjitter1, b2, lnjitter2) = result['x']    
    return m, b1, lnjitter1, b2, lnjitter2
    


if __name__ == "__main__":
    
    root_dir = '/Users/mbedell/Documents/Research/HARPSTwins/Abundances/All/'
    a = np.genfromtxt(root_dir+'final_abundances_w_ncapture.csv', delimiter=',', dtype=None, names=True)
    par = np.genfromtxt(root_dir+"final_parameters.csv", delimiter=',', dtype=None, names=True)
    ages = np.genfromtxt(root_dir+'final_ages_combination.csv', delimiter=',', dtype=None, names=True)
    
    age = ages['age_mean'] - 4.6
    age_err = ages['age_std']
    
    #fit_keep = np.where(age < 6.5)[0] 
    

    fit_keep = [i not in ['HIP19911', 'HIP108158', 'HIP109821', 'HIP115577', 'HIP14501', 'HIP28066', 'HIP30476',
            'HIP33094', 'HIP65708', 'HIP73241', 'HIP74432', 'HIP64150'] for i in a['id'][:-1]] # mask out SB2, thick-disk
    inv = np.invert(fit_keep)   
    
    # get best-fit GCE relations:
    
    elements = ['CI','CH', 'OI','NaI','MgI','AlI','SiI','SI', 'CaI', 'ScI','ScII', 'TiI', 'TiII', 
                'VI','CrI', 'CrII', 'MnI', 'CoI','NiI','CuI','ZnI',
                'SrI', 'YII', 'ZrII', 'BaII', 'LaII', 'CeII', 'PrII', 'NdII', 'SmII', 'EuII', 
                'GdII', 'DyII']
    ms = np.zeros(len(elements)) # slopes
    m_errp = np.zeros_like(ms)
    m_errm = np.zeros_like(ms)
    bs = np.zeros(len(elements)) # intercepts
    b_errp = np.zeros_like(ms)
    b_errm = np.zeros_like(ms)
    lnjs = np.zeros(len(elements)) # ln(jitter)s
    j_errp = np.zeros_like(ms)
    j_errm = np.zeros_like(ms)
    
    
    for i,el in enumerate(elements):
        if el in ['CI', 'ScI', 'TiI', 'CrI']:
            # do a two-species fit
            el2 = elements[i+1]
            abund1 = a[el+"_1"][:-1] - par['feh'][:-1] # exclude sun
            err1 = a["err_"+el][:-1]
            abund2 = a[el2+"_1"][:-1] - par['feh'][:-1] # exclude sun
            err2 = a["err_"+el2][:-1]
            
            abund1 = np.power(10., abund1)
            err1 = err1 * abund1 * np.log(10.)
            
            abund2 = np.power(10., abund2)
            err2 = err2 * abund2 * np.log(10.)
            
            # mask missing values
            fit = np.copy(fit_keep)
            bad = np.where(np.isnan(abund1))[0]
            if len(bad) > 0:
                fit[bad] = False
                
                
            bf = bestfit_2sp(abund1[fit], err1[fit], abund2[fit], err2[fit], age[fit], age_err[fit])
            ms[i], bs[i], lnjs[i] = bf[0:3]
            ms[i+1], bs[i+1], lnjs[i+1] = np.append(bf[0], bf[3:])
        
            # mcmc for slope errors:
            ndim, nwalkers = 5, 20
            spread = [1e-3, 1e-2, 1e-4, 1e-2, 1e-4]
            pos = [bf + spread*np.random.randn(ndim) for j in range(nwalkers)]
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_2sp, args=(abund1[fit], err1[fit], abund2[fit], err2[fit], age[fit], age_err[fit]))
            sampler.run_mcmc(pos, 1000)

            # save mcmc results:
            samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
            samples[:,2] = np.exp(samples[:,2])
            samples[:,-1] = np.exp(samples[:,-1])
            fig = corner.corner(samples, labels=["$m$", "$b_1$", "$j_1$", "$b_2$", "$j_2$"])        
            fig.savefig('{0}_emcee.png'.format(el))
            plt.close()
            
            np.savetxt('{0}_emcee.txt'.format(el), samples, delimiter=',',
                    header='slope, intercept1, jitter1, intercept2, jitter2')
            m_16, m_50, m_84 = np.percentile(samples[:,0], [16,50,84])
            b1_16, b1_50, b1_84 = np.percentile(samples[:,1], [16,50,84])
            j1_16, j1_50, j1_84 = np.exp(np.percentile(samples[:,2], [16,50,84]))
            b2_16, b2_50, b2_84 = np.percentile(samples[:,3], [16,50,84])
            j2_16, j2_50, j2_84 = np.exp(np.percentile(samples[:,4], [16,50,84]))
            m_errp[i] = m_84 - m_50
            m_errm[i] = m_50 - m_16  
            m_errp[i+1] = m_84 - m_50
            m_errm[i+1] = m_50 - m_16  
            b_errp[i] =   b1_84 - b1_50  
            b_errm[i] =   b1_50 - b1_16  
            b_errp[i+1] =   b2_84 - b2_50  
            b_errm[i+1] =   b2_50 - b2_16
            j_errp[i] =   j1_84 - j1_50 
            j_errm[i] =   j1_50 - j1_16  
            j_errp[i+1] =   j2_84 - j2_50 
            j_errm[i+1] =   j2_50 - j2_16
            
            
            print "elements {0} and {1} completed with slope = {2:.3e} +{3:.3e} -{4:.3e}".format(el, el2, ms[i], m_errp[i], m_errm[i])
            print "(best-fit slope - median slope) = {0:.5e}".format(ms[i] - m_50)
            print "(intercept2 - intercept1) = {0:.3e} dex".format(bs[i+1] - bs[i])
            print "jitter1 = {0:.3e} dex, jitter2 = {1:.3e} dex".format(np.exp(lnjs[i]), np.exp(lnjs[i+1]))
            print "---------------------------------------------"
            
        elif el in ['CH', 'ScII', 'TiII', 'CrII']:
            continue
            
        else:
            # do a one-species fit
            abund = a[el+"_1"][:-1] - par['feh'][:-1] # exclude sun
            err = a["err_"+el][:-1]
            
            abund = np.power(10., abund)
            err = err * abund * np.log(10.)
            
            # mask missing values
            fit = np.copy(fit_keep)
            bad = np.where(np.isnan(abund))[0]
            if len(bad) > 0:
                fit[bad] = False
            
            bf = bestfit_1sp(abund[fit], err[fit], age[fit], age_err[fit])
            ms[i], bs[i], lnjs[i] = bf
        
            # mcmc for slope errors:
            ndim, nwalkers = 3, 20
            spread = [1e-3, 1e-2, 1e-4]
            pos = [bf + spread*np.random.randn(ndim) for j in range(nwalkers)]
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_1sp, args=(abund[fit], err[fit], age[fit], age_err[fit]))
            sampler.run_mcmc(pos, 1000)

            # save mcmc results:
            samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
            samples[:,-1] = np.exp(samples[:,-1])
            fig = corner.corner(samples, labels=["$m$", "$b$", "$j$"])        
            fig.savefig('{0}_emcee.png'.format(el))
            plt.close()
            
            np.savetxt('{0}_emcee.txt'.format(el), samples, delimiter=',',
                    header='slope, intercept, jitter')
            m_16, m_50, m_84 = np.percentile(samples[:,0], [16,50,84])
            b_16, b_50, b_84 = np.percentile(samples[:,1], [16,50,84])
            j_16, j_50, j_84 = np.exp(np.percentile(samples[:,2], [16,50,84]))
            m_errp[i] = m_84 - m_50
            m_errm[i] = m_50 - m_16 
            b_errp[i] =   b_84 - b_50  
            b_errm[i] =   b_50 - b_16  
            j_errp[i] =   j_84 - j_50 
            j_errm[i] =   j_50 - j_16     
            print "element {0} completed with slope = {1:.3e} +{2:.3e} -{3:.3e}".format(el, ms[i], m_errp[i], m_errm[i])
            print "(best-fit slope - median slope) = {0:.5e}".format(ms[i] - m_50)
            print "---------------------------------------------"
            

        
    # save:
    js = np.exp(lnjs)
    np.savetxt('gce_logspace.txt', np.transpose([elements, ms, m_errp, m_errm, bs, b_errp, b_errm, js, j_errp, j_errm]), fmt='%s', delimiter=',',
                header='element, slope, slope_errp, slope_errm, intercept, intercept_errp, intercept_errm, jitter, jitter_errp, jitter_errm')
    
    # plot it:
    if True:
    
        c2 = '#003399' # blue
        c3 = '#CC0033' # red
        c4 = '#339900' # green
        plt.rcParams["font.sans-serif"] = "Helvetica"
        xs = np.arange(-5., 6., 1.0)

        fig = plt.figure(figsize=(20,30))

        for i,el in enumerate(elements):
            ax = fig.add_subplot(9,4,i+1)

            abund = a[el+"_1"][:-1] - par['feh'][:-1] # exclude sun
            err = a["err_"+el][:-1]
            
            abund = np.power(10., abund)
            err = err * abund * np.log(10.)

            ax.errorbar(age[inv], abund[inv], xerr=age_err[inv], yerr=err[inv], fmt='^', c=c3, ecolor=c3)
            ax.errorbar(age[fit_keep], abund[fit_keep], xerr=age_err[fit_keep], yerr=err[fit_keep], fmt='o', c='black', ecolor='black', mec='black')
            ax.annotate(r'$\odot$', xy=(0.0, 1.0), horizontalalignment='center', verticalalignment='center', color=c4, fontsize=24, weight='bold')


            ax.plot(xs, ms[i]*xs + bs[i], color=c2, lw=2)

            #ax.set_ylim([-0.2,0.3])
            ax.text(-2.0,1.4,el)

            #ax.set_yticks(np.arange(-0.2,0.3,0.1))
            #ax.set_yticks(np.arange(-0.2,0.3,0.05), minor=True)
            
            ax.set_ylim([0.6,1.6])
            ax.set_yticks(np.arange(0.5,1.6,0.5))
            ax.set_yticks(np.arange(0.6,1.6,0.1), minor=True)

            ax.set_xticks(np.arange(-4,6,2))
            ax.set_xticks(np.arange(-5,6,0.5), minor=True)

            ax.tick_params(axis='both', which='major', labelsize=14)

            if (i % 4) != 0:
                ax.set_yticklabels('',visible=False)

            if el not in elements[-4:]:
                ax.set_xticklabels('',visible=False)

        fig.subplots_adjust(hspace=.05, wspace=.05)
        fig.text(0.5, 0.07, r'Age - Age$_{\odot}$ (Gyr)', size=28, ha='center')
        fig.text(0.05, 0.5, r'$N_{X} - N_{X, \odot}$', rotation=90, size=28, va='center')
        fig.savefig('gce_linspace.pdf')
    