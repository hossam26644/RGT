'''docstring'''
import matplotlib.pyplot as plt

def plot_counts_table(counts_table, export_directory, sample_code, first_allele, second_allele, color_code='black'):
    #print(first_allele.repeat_units_count)
    color_codes = {"red":"r", "green":"g" , "yellow":"orange", "black":"k"}
    sorted_table = sorted(counts_table.items(),reverse=True)
    try:
        x, y = zip(*sorted_table)

    except ValueError as e:
        x = [50]
        y = [0]
    
    x_ticks_scaling_factor = (max(x)//40)+1


    graph = plt.bar(x,y, align='center', alpha=0.65)
    try:
        first_allele_index = x.index(first_allele.repeat_units_count)
        second_allele_index = x.index(second_allele.repeat_units_count)
        graph[first_allele_index].set_facecolor('r')
        if first_allele != second_allele:
            graph[second_allele_index].set_facecolor('#EC7063')
    except:
        pass

    plt.xticks(list(range(1,max(x),x_ticks_scaling_factor)),
               list(range(1,max(x),x_ticks_scaling_factor)),
               fontsize=6 , rotation=30, fontweight='medium' ) #4.5
    

    plt.yticks(fontsize=6, fontweight='medium')

    plt.title(sample_code, color=color_codes[color_code])
    plt.xlabel("Total number of repeat units" )
    plt.ylabel("Number of reads" )
    try:
        ax = plt.gca()
        if first_allele != second_allele :
            ax.legend((graph[first_allele_index], graph[second_allele_index] ),
                ([first_allele.sequence_string, second_allele.sequence_string]),
                 fontsize=4, borderaxespad=0, frameon=False, loc='upper left')
        else:
            ax.legend((graph[first_allele_index],),
               (first_allele.sequence_string,),fontsize=4, borderaxespad=0, frameon=False, loc='upper left')
    except:
        pass

    plt.savefig(export_directory, dpi=600)
    plt.clf()
