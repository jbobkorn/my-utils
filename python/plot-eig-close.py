from __future__ import print_function, division
from numpy import *
from matplotlib.pyplot import *
import analyse_eig as eig

### INPUT
filename = sys.argv[-1]
#extract important eigenvlaues from all eigenvalues
#short form (without naming all objects inside observables):
observables = eig.observe(filename, PRINT=False)
analyzed = eig.analyse(filename, PRINT=False) #analyzed is just an extended version of observe (used for bar plots, see below)


### PLOT
params = {'legend.fontsize': '7',
             'figure.figsize': (6, 3),
           'axes.labelsize': 10,
            'axes.titlesize':'xx-large',
           'xtick.labelsize': 10,
           'ytick.labelsize': 10}
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rcParams.update(params)


### Example: Plot values of one EIG file, around fermi energy, spin resolved, with information on the picture 
#requires: observables
fig0 = eig.plotclose(observables, params, 0)
fig0.savefig('../Eig_files_plots/{}.pdf'.format(filename[:-4] + '-EIGzoom'), bbox_inches='tight')
