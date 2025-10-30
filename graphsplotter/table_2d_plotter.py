'''docstring'''
import matplotlib.pyplot as plt

def plot_2d_table(table, export_directory, sample_code, first_allele, second_allele,
                first_allele_count, second_allele_count,xlabel, color_code='black'):

    color_codes = {"red":"r", "green":"g" , "yellow":"orange", "black":"k"}
    sorted_table = sorted(table.items(),reverse=True)
    try:
        x, y = zip(*sorted_table)

    except ValueError as e: #when the counts table is empty, plot an empty plot
        x = [50]
        y = [0]
    
    x_ticks_scaling_factor = (max(x)//40)+1 # scale the number of x-ticks to avoid their overlap 


    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    bars = ax.bar(x, y, align="center", alpha=0.65)
    
    try:
        #get the index of the alleles in the x data
        first_allele_index = x.index(first_allele_count) 
        second_allele_index = x.index(second_allele_count)
        
        #color the bars that correspond to the alleles
        bars[second_allele_index].set_facecolor('#EC7063')
        bars[first_allele_index].set_facecolor('r')
        
    except:
        pass

    # manipulating the xticks, font, rotation,and step
    ticks = list(range(1, max(x), x_ticks_scaling_factor))
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks, fontsize=6, rotation=30, fontweight="medium")    
    ax.tick_params(axis="y", labelsize=6)
    ax.set_title(sample_code, color=color_codes.get(color_code, "k"))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of reads")
    
    try: #all the hassle in this try is to put the legend
        #ax = plt.gca()
        if first_allele != second_allele : #heterozygous
            if first_allele_index == second_allele_index:
                
                two_alleles_counts = first_allele.abundance + second_allele.abundance
                first_allele_percentage = first_allele.abundance / two_alleles_counts
                second_allele_percentage = second_allele.abundance / two_alleles_counts

                first_legend_string = first_allele.sequence_string +" ("+ str(round(first_allele_percentage, 1)*100) + "%)" 
                second_legend_string = second_allele.sequence_string +" ("+ str(round(second_allele_percentage, 1)*100) + "%)"
                
                ax.legend((bars[first_allele_index], bars[second_allele_index] ),
                    ([first_legend_string, second_legend_string]),
                     fontsize=4, borderaxespad=0, frameon=False, loc='upper left')
            else:
                ax.legend((bars[first_allele_index], bars[second_allele_index] ),
                    ([first_allele.sequence_string, second_allele.sequence_string]),
                     fontsize=4, borderaxespad=0, frameon=False, loc='upper left')
        else: #homozygous
            ax.legend((bars[first_allele_index],),
               (first_allele.sequence_string,),fontsize=4, borderaxespad=0, frameon=False, loc='upper left')
    except:
        pass

    fig.tight_layout()
    fig.savefig(export_directory, dpi=600)
    plt.close(fig)
