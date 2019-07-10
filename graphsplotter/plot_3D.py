'''docstring'''
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import


def plot_3D(table, export_directory, sample_code, first_allele, second_allele,
            xlabel, zlabel, color_code='black'):

    color_codes = {"red":"r", "green":"g" , "yellow":"orange", "black":"k"}
    xvals = []; yvals = []; zvals=[]

    for key in table.keys():
        xvals.append(key[0])
        yvals.append(key[1])
        zvals.append(table[key])
    
    max_value = (max(max(xvals),max(yvals))) #to create am empty cube
    fig = plt.figure()
    ax = fig.add_subplot( projection='3d')
    
    for i in range(max_value): #create an empty square
        ax.bar3d(i, i, 0, 0.00, 0.00, 0,color='w', shade=True,edgecolor='w',alpha=0.0)



    for i in range(len(xvals)): #plot all values
        ax.bar3d(xvals[i]-0.5, yvals[i]-0.5, 0, 0.5, 0.5, zvals[i],color='#2E86C1', shade=True,alpha=1)
    
    first_allele_index = xvals.index(first_allele.x_count) 
    second_allele_index = xvals.index(second_allele.x_count)  

    ax.bar3d(xvals[second_allele_index]-0.5, yvals[second_allele_index]-0.5,
            0, 0.5, 0.5, zvals[second_allele_index],color='#EC7063', shade=True,alpha=1)

    ax.bar3d(xvals[first_allele_index]-0.5, yvals[first_allele_index]-0.5,
            0, 0.5, 0.5, zvals[first_allele_index],color='#c0392b', shade=True,alpha=1)
    
    first_proxy = plt.Rectangle((0, 0), 1, 1, fc="#c0392b")
    second_proxy = plt.Rectangle((0, 0), 1, 1, fc="#EC7063")
    frst_legnd = first_allele.sequence_string + "  ("+xlabel+ ": "+ str(xvals[first_allele_index]) +\
                ") , ("+zlabel+": "+str(zvals[first_allele_index])+")" 

    scnd_legnd = second_allele.sequence_string + "  ("+xlabel+ ": "+ str(xvals[second_allele_index]) +\
                ") , ("+zlabel+": "+str(zvals[second_allele_index])+")" 
    if first_allele == second_allele:
        frst_legnd += str(first_allele.abundance/first_allele.abundance+second_allele.abundance)
        scnd_legnd += str(second_allele.abundance/first_allele.abundance+second_allele.abundance)
    
    ax.legend((first_proxy,second_proxy),
           (frst_legnd,scnd_legnd),
           fontsize=6, borderaxespad=0, frameon=False, loc='upper left')

    ax.set_title(sample_code,loc='center',pad=35, color=color_codes[color_code])

    ax.set_xlabel(xlabel)
    ax.set_ylabel(zlabel)
    ax.set_zlabel('Abundance')

    plt.autoscale(enable=True, axis='both', tight=True)
    #plt.show()
    plt.savefig(export_directory, dpi=600)
    plt.clf()
