'''docstring'''
import matplotlib.pyplot as plt

def plot_2d_table(counts_table, export_directory, sample_code, first_allele, second_allele, color_code='black'):

    color_codes = {"red":"r", "green":"g" , "yellow":"orange", "black":"k"}
    sorted_table = sorted(counts_table.items(),reverse=True)
    try:
        x, y = zip(*sorted_table)

    except ValueError as e: #when the counts table is empty, plot an empty plot
        x = [50]
        y = [0]
    
    x_ticks_scaling_factor = (max(x)//40)+1 # scale the number of x-ticks to avoid their overlap 


    graph = plt.bar(x,y, align='center', alpha=0.65)
    
    try:
        #get the index of the alleles in the x data
        first_allele_index = x.index(first_allele.repeat_units_count) 
        second_allele_index = x.index(second_allele.repeat_units_count)
        
        #color the bars that correspond to the alleles
        graph[second_allele_index].set_facecolor('#EC7063')
        graph[first_allele_index].set_facecolor('r')
        
    except:
        pass

    # manipulating the xticks, font, rotation,and step
    plt.xticks(list(range(1,max(x),x_ticks_scaling_factor)),
               list(range(1,max(x),x_ticks_scaling_factor)),
               fontsize=6 , rotation=30, fontweight='medium' ) #4.5
    

    plt.yticks(fontsize=6, fontweight='medium')

    plt.title(sample_code, color=color_codes[color_code])
    plt.xlabel("Total number of repeat units" )
    plt.ylabel("Number of reads" )
    
    try: #all the hassle in this try is to put the legend
        ax = plt.gca()
        if first_allele != second_allele : #heterozygous
            if first_allele_index == second_allele_index:
                
                two_alleles_counts = first_allele.abundance + second_allele.abundance
                first_allele_percentage = first_allele.abundance / two_alleles_counts
                second_allele_percentage = second_allele.abundance / two_alleles_counts

                first_legend_string = first_allele.sequence_string +" ("+ str(round(first_allele_percentage, 1)*100) + "%)" 
                second_legend_string = second_allele.sequence_string +" ("+ str(round(second_allele_percentage, 1)*100) + "%)"
                
                ax.legend((graph[first_allele_index], graph[second_allele_index] ),
                    ([first_legend_string, second_legend_string]),
                     fontsize=4, borderaxespad=0, frameon=False, loc='upper left')
            else:
                ax.legend((graph[first_allele_index], graph[second_allele_index] ),
                    ([first_allele.sequence_string, second_allele.sequence_string]),
                     fontsize=4, borderaxespad=0, frameon=False, loc='upper left')
        else: #homozygous
            ax.legend((graph[first_allele_index],),
               (first_allele.sequence_string,),fontsize=4, borderaxespad=0, frameon=False, loc='upper left')
    except:
        pass

    plt.savefig(export_directory, dpi=600)
    plt.clf()
