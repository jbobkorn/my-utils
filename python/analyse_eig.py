from __future__ import print_function, division
from numpy import *
from matplotlib.pyplot import *


CONST_close_lower = 6 #must be twice the same number here, close_up_range and close_down_range must have same length
CONST_close_upper = 6

### FUNCTIONS
def read_data( filename ):
    """ Function that read the results from a data file and returns the
    all data as an array of lists which are lines in data file."""

    # read first word at first line
    data = []
    with open( filename, 'r' ) as f:
        for line in f:
            line = line.split()
            if line:
                #important, float() to convert each string number to a float number
                line = [float(i) for i in line]
                data.append(line)
    data = array(data)
    return data

# function to output fermi energy
def fermi_energy(data):
    return data[0]
# function to output all eigen energies in an array
def eigs(data):
    eigs = []
    N = len(data)
    L = len(data[2])
    for k in arange(1,L):
        eigs.append(data[2][k])
    for i in arange(3,N):
        M = len(data[i])
        for j in arange(M):
            eigs.append(data[i][j])
    eigs = array(eigs)

    return eigs
# function to put all the eigen energies above/below fermi together in an array
# and return the minimal eig of CB & the maximal eig of VB
def eig_bonds(eigs,e_f_value):
    eigs_above = []
    eigs_below = []
    for e in eigs:
        if e > e_f_value:
            eigs_above.append(e)
        if e < e_f_value:
            eigs_below.append(e)

    eigs_above = array(eigs_above)
    eigs_below = array(eigs_below)
    return min(eigs_above), max(eigs_below),eigs_above,eigs_below

def degeneracy_check(index, array):
    """
    checks if energy levels below fermi energy are degenerate and returns the biggest index, whose energy is below E_f. (only the index below E_f must be checked, since eig_bonds always returns the smallest matching index, which is unwished below E_f but wished above E_f.
    """
    while array[index] == array[index+1]:
        index +=1
    return(index)


def observe(filename, PRINT=False):
    """Groups the eigenvalues into observables"""
    #call functions defined above
    #aquire data from file
    data = read_data(filename)
    #extract eigenvalues from data
    eigens = eigs(data)
    # split the total eig array into two
    #spin up eigs array
    eig_up = split(eigens,2)[0]
    #spin down eigs array
    eig_down = split(eigens,2)[1]
    N = len(eig_up)
    e_up = eig_up.tolist()
    e_down = eig_down.tolist()
    #fermi energy
    e_f_value = fermi_energy(data)
    
    #minimal energy in conduction band in spin up:
    e_up_min = eig_bonds(eig_up,e_f_value)[0]
    index_up_min = e_up.index(e_up_min)
    ##maximal energy in valence band in spin up:
    e_up_max = eig_bonds(eig_up,e_f_value)[1]
    index_up_max = degeneracy_check(e_up.index(e_up_max), e_up)
    #band gap = e_LUMO - e_HOMO
    e_g_up = e_up_min - e_up_max
    
    #minimal energy in conduction band in spin down:
    e_down_min = eig_bonds(eig_down,e_f_value)[0]
    index_down_min = e_down.index(e_down_min)
    ##maximal energy in valence band in spin down:
    e_down_max = eig_bonds(eig_down,e_f_value)[1]
    index_down_max = degeneracy_check(e_down.index(e_down_max), e_down)
    #band gap = e_LUMO - e_HOMO
    e_g_down = e_down_min - e_down_max
    
    close_up_range = range(index_up_min-CONST_close_lower, index_up_max+CONST_close_upper+1) #must have same length as close_down_range
    close_down_range = range(index_down_min-CONST_close_lower, index_down_max+CONST_close_upper+1) #must have same length as close_up_range
    
    if PRINT == True:
        print('Total indices per spin: N = {}'.format(N))
        print('')
        print('Fermi Energy: E_F = {}'.format(e_f_value[0]))
        print('')
        print('E_up_LUMO = {}, index = {}'.format(e_up_min,index_up_min+1))
        print('E_up_HOMO = {}, index = {}'.format(e_up_max,index_up_max+1))
        print('E_g_up = {}'.format(e_g_up))
        print('')
        print('E_down_LUMO = {}, index:{}'.format(e_down_min,index_down_min+N+1))
        print('E_down_HOMO = {}, index:{}'.format(e_down_max,index_down_max+N+1))
        print('E_g_down = {}'.format(e_g_down))
        print('')
        print('Values around the fermi energy:')
        print('')
        print('UP: (i, e_i); i counting from 1')
        for i in close_up_range:
            print((i+1, e_up[i]))
            if i+1 == index_up_min:
                print("-----------------------")
        print('')
        print('DOWN: (i, e_i); i counting from 1')
        for i in close_down_range:
            print((i+1+N, e_down[i])) #just PRINT the index counting from 1, but give the right energy for the index from 0.
            if i+1 == index_down_min:
                print("-----------------------")
    else:
        pass
    
    return(eig_up, eig_down, N, e_f_value[0], [e_up_min,index_up_min+1], [e_up_max,index_up_max+1], e_g_up, [e_down_min,index_down_min+N+1], [e_down_max,index_down_max+N+1], e_g_down, close_up_range, close_down_range)


def analyse(f, PRINT=False): #f is a filename
    "An advaned version of the observe function above."
    print(f)
    #call functions defined above
    data = read_data(f)
    eigens = eigs(data)

    # split the total eig array into two
    #spin up eigs array
    eig_up = split(eigens,2)[0]
    #spin down eigs array
    eig_down = split(eigens,2)[1]
    N = len(eig_up)
    e_up = eig_up.tolist()
    e_down = eig_down.tolist()
    #fermi energy
    e_f_value = fermi_energy(data)

    #minimal energy in conduction band in spin up:
    e_up_min = eig_bonds(eig_up,e_f_value)[0]
    index_up_min = e_up.index(e_up_min)
    ##maximal energy in valence band in spin up:
    e_up_max = eig_bonds(eig_up,e_f_value)[1]
    index_up_max = degeneracy_check(e_up.index(e_up_max), e_up)
    #band gap = e_LUMO - e_HOMO
    e_g_up = e_up_min - e_up_max

    #minimal energy in conduction band in spin down:
    e_down_min = eig_bonds(eig_down,e_f_value)[0]
    index_down_min = e_down.index(e_down_min)
    ##maximal energy in valence band in spin down:
    e_down_max = eig_bonds(eig_down,e_f_value)[1]
    index_down_max = degeneracy_check(e_down.index(e_down_max), e_down)
    #band gap = e_LUMO - e_HOMO
    e_g_down = e_down_min - e_down_max


    close_up_range = range(index_up_min-CONST_close_lower, index_up_max+CONST_close_upper+1) #must have same length as close_down_range
    close_down_range = range(index_down_min-CONST_close_lower, index_down_max+CONST_close_upper+1) #must have same length as close_up_range


    x = linspace(1,N,N)
    #here just for plotting the fermi line
    fermi = []
    for i in arange(N):
        fermi.append(e_f_value)
    fermi = array(fermi)

    x_up_close = x[close_up_range[0]:close_up_range[-1]+1]
    eig_up_close = eig_up[close_up_range[0]:close_up_range[-1]+1]
    x_down_close = x[close_down_range[0]:close_down_range[-1]+1]
    eig_down_close = eig_down[close_down_range[0]:close_down_range[-1]+1]
    x_fs = x[min(close_up_range[0], close_down_range[0]):max(close_up_range[-1], close_down_range[-1])]
    e_fs = fermi[min(close_up_range[0], close_down_range[0]):max(close_up_range[-1], close_down_range[-1])] 
    
    if PRINT == True:
        print('Total indices per spin: N = {}'.format(N))
        print('')
        print('Fermi Energy: E_F = {}'.format(e_f_value[0]))
        print('')
        print('E_up_LUMO = {}, index = {}'.format(e_up_min,index_up_min+1))
        print('E_up_HOMO = {}, index = {}'.format(e_up_max,index_up_max+1))
        print('E_g_up = {}'.format(e_g_up))
        print('')
        print('E_down_LUMO = {}, index:{}'.format(e_down_min,index_down_min+N+1))
        print('E_down_HOMO = {}, index:{}'.format(e_down_max,index_down_max+N+1))
        print('E_g_down = {}'.format(e_g_down))
        print('')
        print('Values around the fermi energy:')
        print('')
        print('UP: (i, e_i); i counting from 1')
        for i in close_up_range:
            print((i+1, e_up[i]))
            if i+1 == index_up_min:
                print("-----------------------")
        print('')
        print('DOWN: (i, e_i); i counting from 1')
        for i in close_down_range:
            print((i+1+N, e_down[i])) #just PRINT the index counting from 1, but give the right energy for the index from 0.
            if i+1 == index_down_min:
                print("-----------------------")
    else:
        pass
    
    return([x, eig_up, eig_down, index_up_min, index_up_max, e_up_min, e_up_max, index_down_min, index_down_max, e_down_min, e_down_max, fermi, e_f_value, N, x_fs, e_fs, x_up_close, x_down_close, eig_up_close, eig_down_close, e_g_up, e_g_down]) #those are 22 output values.


def get_analysis_dict(filenames):
    """Needs a list of strings (filename paths) as input"""
    dct = {} #dct to hold all parameters for the filenames
    for filename in filenames:
        par = analyse(filename)
        e_f = par[-10][0]
        e_g_up = par[-2]
        e_g_down = par[-1]
        e_close_up, e_close_down = par[-4], par[-3]
        e_close_up_shifted, e_close_down_shifted = fermi_shifter(e_close_up, e_f), fermi_shifter(e_close_down, e_f) #energy array, shifted by the fermi energy
        dct['{}'.format(filename)] = [e_f, e_close_up_shifted, e_close_down_shifted, e_g_up, e_g_down]
    return(dct)



def fermi_shifter(array, e_f):return(array-e_f)

def plot_energy_bars(E_array, linestyle, color, label):
    N = len(E_array)
    for i in range(N):
        bar_x = array([i-0.25, i+0.25])
        Ei = E_array[i]
        bar_E = zeros(len(bar_x)) + Ei
        if i == 1: #just add one label. 
            plot(bar_x, bar_E, linestyle=linestyle, color=color, label=label)
        else: 
            plot(bar_x, bar_E, linestyle=linestyle, color=color, label='')
    return(0)

    
def plotclose(observe_function_output, plot_parameters, fignum):
    params = plot_parameters
    [E_up, E_down, N, E_F, [E_up_LUMO, index_up_LUMO], [E_up_HOMO, index_up_HOMO], E_g_up, [E_down_LUMO,index_down_LUMO], [E_down_HOMO,index_down_HOMO], E_g_down, close_up_range, close_down_range] = observe_function_output
    x = linspace(1,N,N) #here just for plotting the fermi line
    fermi = []
    for i in arange(N):
        fermi.append(E_F)
    fermi = array(fermi)
    fig = figure(fignum)
    xlabel(r'$i_{\mathrm{up}}-1;\ i_{\mathrm{down}}-N-1$')
    ylabel(r'$E_i$')
    #Plot Eig values around fermi energy
    #SPIN UP:
    plot(x[close_up_range[0]:close_up_range[-1]+1],E_up[close_up_range[0]:close_up_range[-1]+1], "^",markersize=4.0, color='blue',label=r'$E_\mathrm{} \ (E^\mathrm{}_\mathrm{} = {})$'.format('{\uparrow}', '{\uparrow}', '{g}', E_g_up))
    #Plot HOMO level
    plot(index_up_HOMO, E_up_HOMO, "^",markersize=4.0, color='blue',label=r'$E_\mathrm{}={}$, $i=${}'.format('{\uparrow,HOMO}',E_up_HOMO,index_up_HOMO))
    #Plot LUMO level
    plot(index_up_LUMO, E_up_LUMO, "^",markersize=4.0, color='blue',label=r'$E_\mathrm{}={}$, $i=${}'.format('{\uparrow,LUMO}', E_up_LUMO, index_up_LUMO))
    
    #SPIN DOWN:
    plot(x[close_down_range[0]:close_down_range[-1]+1],E_down[close_down_range[0]:close_down_range[-1]+1], "v",markersize=4.0, color='red',label=r'$E_\mathrm{} \ (E^\mathrm{}_\mathrm{} = {})$'.format('{\downarrow}', '{\downarrow}', '{g}', E_g_down))
    #Plot HOMO level (and wrap to the same index as spin UP)
    plot(index_down_HOMO-N, E_down_HOMO, "v",markersize=4.0, color='red',label=r'$E_\mathrm{}={}$, $i=${}'.format('{\downarrow,HOMO}', E_down_HOMO, index_down_HOMO))
    #Plot LUMO level (and wrap to the same index as spin UP)
    plot(index_down_LUMO-N, E_down_LUMO, "v",markersize=4.0, color='red',label=r'$E_\mathrm{}={}$, $i=${}'.format('{\downarrow,LUMO}', E_down_LUMO, index_down_LUMO))
    #Plot a line for the Fermi level
    plot(x[min(close_up_range[0], close_down_range[0]):max(close_up_range[-1], close_down_range[-1])+1],fermi[min(close_up_range[0], close_down_range[0]):max(close_up_range[-1], close_down_range[-1])+1], "--", linewidth=0.5,color='black',label=r'$E_F={}$'.format(E_F))
    legend(loc='best')
    #show()
    return(fig)
    

def plotwide(analyse_function_output, plot_parameters, fignum): 
    '''call by using plotwide(analyse(filename))'''
    params = plot_parameters
    [x, eig_up, eig_down, index_up_min, index_up_max, e_up_min, e_up_max, index_down_min, index_down_max, e_down_min, e_down_max, fermi, e_f_value, N, x_fs, e_fs, x_up_close, x_down_close, eig_up_close, eig_down_close, e_g_up, e_g_down] = analyse_function_output
    fig = figure(fignum)
    xlabel(r'$i_{\mathrm{up}}-1;\ i_{\mathrm{down}}-N-1$')
    ylabel(r'$E_i$')
    plot(x,eig_up, "^",markersize=4.0, color='red',label=r'$E_\mathrm{UP}$')
    plot(x,eig_down, "v",markersize=4.0, color='blue',label=r'$E_\mathrm{DOWN}$')
    plot(index_up_min+1,e_up_min, "^",markersize=4.0, color='cyan',label=r'$E_\mathrm{}={}$ in CB, $i=${}'.format('{upmin}',e_up_min,index_up_min+1))
    plot(index_up_max+1,e_up_max, "^",markersize=4.0, color='cyan',label=r'$E_\mathrm{}={}$ in VB, $i=${}'.format('{upmax}',e_up_max,index_up_max+1))
    plot(index_down_min+1,e_down_min, "v",markersize=4.0, color='magenta',label=r'$E_\mathrm{}={}$ in CB, $i=${}'.format('{downmin}',e_down_min,index_down_min+N+1))
    plot(index_down_max+1,e_down_max, "v",markersize=4.0, color='magenta',label=r'$E_\mathrm{}={}$ in VB, $i=${}'.format('{downmax}',e_down_max,index_down_max+N+1))
    plot(x,fermi, "--", linewidth=0.5,color='black',label=r'$E_f={}$'.format(e_f_value[0]))
    legend(loc='best')
    show()
    return(fig)


def mysort(array):
    a = array
    flat_a = [item for sublist in a for item in sublist]
    print('flat_a: {}'.format( flat_a))
    mergepoints = []
    for j in range(len(a)-1):
        if not a[j][-1] == a[j+1][0]:
            first_flat_index = flat_a.index(a[j+1][0])
            last_flat_index = len(flat_a)-1 - flat_a[::-1].index(a[j+1][0])
            #print('first flat index {}'.format(first_flat_index))
            #print('last flat index {}'.format(last_flat_index))
            #mergepoints.append(first_flat_index)
            mergepoints.append(last_flat_index-1)
    mergepoints.insert(0,-1) #to sort from beginning to first mergepoint
    mergepoints.append(len(flat_a)) #to sort from last mergepoint to end
    print(mergepoints)
    res = []
    for i in range(len(mergepoints)-1):
        mi = [mergepoints[i], mergepoints[i+1]] #array, from-to
        #print(i,mi)
        #print(flat_a[mi[0]+1:mi[1]+1]) #array, from-to of the input array (with doubles)
        #print(sorted(list(set(flat_a[mi[0]+1:mi[1]+1])))) #remove doubles and sort the set of numbers
        res.append(sorted(list(set(flat_a[mi[0]+1:mi[1]+1]))))
    #mergepoints fails, if there is only 2-fold degeneracy. in that case, a is already the array we want.
    if flat_a == res[0]:
        res = a
    return(res)

def degeneracy_with_tolerance(E_array, tolerance):
    deg = []
    N = len(E_array)
    for i in range(N-1):
        one_deg = []
        if abs(E_array[i] - E_array[i+1]) < tolerance:
            one_deg.append(i)
            one_deg.append(i+1)
        if len(one_deg) != 0:
            deg.append(one_deg)
    #now deg will have degenerate eigenenergies as a list of lists, where sublist always list pairs of two degenerate energies. e.g.: [[1,2],[4,5],[5,6],...] here. we want to have [[1,2],[4,5,6],...] to get a degeneracy of two and a degeneracy of 3, represented as lists
    print('deg: {}'.format(deg))
    formatted_deg = mysort(deg)
    print('formatted_deg: {}'.format(formatted_deg))
    return(formatted_deg)


def get_degenerates(E_array, tolerance):
    N = len(E_array)
    degenerates = degeneracy_with_tolerance(E_array, tolerance)
    degenerate_energies = []
    for i in range(len(degenerates)):
        degenerate_energies.append(E_array[degenerates[i][0]])
    #print(degenerates)
    #print(degenerate_energies)
    return(degenerates, degenerate_energies)

def plot_eig_diagram(xpos, E_array, linestyle, color, label='', markersize=4.0):
    N = len(E_array)
    degeneracy_tolerance = 0.01 #set tolerance at 0.01eV
    deg, deg_e = get_degenerates(E_array, degeneracy_tolerance)     #print('Degenerate indices: {} \n With energies: {}'.format(deg,deg_e))
    #plot the degenerate energies
    for one_deg in deg:
        #print(one_deg)
        n_fold_degeneracy = len(one_deg)
        start_from = xpos - (n_fold_degeneracy - 1)/2
        for i in one_deg:
            bar_x = array([start_from-0.25, start_from+0.25]) #define constant bar -> will plot bars above each other
            Ei = E_array[i]
            bar_E = zeros(len(bar_x)) + Ei
            if i == 1: #just add one label. 
                plot(bar_x, bar_E, linestyle=linestyle, markersize=markersize, color=color, label=label)
            else: 
                plot(bar_x, bar_E, linestyle=linestyle, markersize=markersize, color=color, label='')
            start_from += 1
    #plot the non degenerate energies
    flat_deg = [item for sublist in deg for item in sublist]
    non_deg = []

    for i in range(N):
        if i not in flat_deg:
            non_deg.append(i)
    #print('Non degenerate indices: {}'.format(non_deg))
    for i in non_deg:
        bar_x = array([xpos-0.25, xpos+0.25]) #define constant bar -> will plot bars above each other
        Ei = E_array[i]
        bar_E = zeros(len(bar_x)) + Ei
        if i == 1: #just add one label. 
            plot(bar_x, bar_E, linestyle=linestyle, markersize=markersize, color=color, label=label)
        else: 
            plot(bar_x, bar_E, linestyle=linestyle, markersize=markersize, color=color, label='')
    return(0)



