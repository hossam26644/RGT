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
        if xvals[i] == first_allele.x_count and yvals[i] == first_allele.z_count:
            ax.bar3d(first_allele.x_count, first_allele.z_count-0.5,
                0, 0.5, 0.5, table[(first_allele.x_count,first_allele.z_count)],
                color='#c0392b', shade=True,alpha=1)
        elif xvals[i] == second_allele.x_count and yvals[i] == second_allele.z_count:
            ax.bar3d(second_allele.x_count, second_allele.z_count-0.5,
                0, 0.5, 0.5, table[(second_allele.x_count,second_allele.z_count)],
                color='#EC7063', shade=True,alpha=1)
        else:
            ax.bar3d(xvals[i], yvals[i]-0.5, 0, 0.5, 0.5, zvals[i],color='#2E86C1', shade=True,alpha=1)

    
    first_proxy = plt.Rectangle((0, 0), 1, 1, fc="#c0392b")
    second_proxy = plt.Rectangle((0, 0), 1, 1, fc="#EC7063")
    frst_legnd = first_allele.sequence_string + "  ("+xlabel+ ": "+ str(first_allele.x_count) +\
                ") , ("+zlabel+": "+str(first_allele.z_count)+")" 

    scnd_legnd = second_allele.sequence_string + "  ("+xlabel+ ": "+ str(second_allele.x_count) +\
                ") , ("+zlabel+": "+str(second_allele.z_count)+")" 
    
    if (first_allele != second_allele and first_allele.x_count == second_allele.x_count
        and first_allele.z_count == second_allele.z_count): #two allels same counts
        frst_legnd += " "+str((first_allele.abundance/(first_allele.abundance+second_allele.abundance))*100)+"%"
        scnd_legnd += " "+str((second_allele.abundance/(first_allele.abundance+second_allele.abundance))*100)+"%"
        second_proxy = first_proxy

    if first_allele == second_allele: #homozygous
        ax.legend((first_proxy,),
           (frst_legnd,),
           fontsize=6, borderaxespad=0, frameon=False, loc='upper left')
    else:    
        ax.legend((first_proxy,second_proxy),
           (frst_legnd,scnd_legnd),
           fontsize=6, borderaxespad=0, frameon=False, loc='upper left')
    ticks_scaling_factor = (max_value//15)+1
    plt.xticks(list(range(1,max_value,ticks_scaling_factor)),
               list(range(1,max_value,ticks_scaling_factor)),
               fontsize=6 , fontweight='medium' ) #4.5

    plt.yticks(list(range(1,max_value,ticks_scaling_factor)),
               list(range(1,max_value,ticks_scaling_factor)),
               fontsize=6 , fontweight='medium' ) #4.5

    ax.set_title(sample_code,loc='center',pad=35, color=color_codes[color_code])

    ax.set_xlabel(xlabel)
    ax.set_ylabel(zlabel)
    ax.set_zlabel('Abundance')

    plt.autoscale(enable=True, axis='both', tight=True)
    #plt.show()
    plt.savefig(export_directory, dpi=600)
    plt.clf()
