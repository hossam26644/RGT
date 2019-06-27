'''docstring'''
import matplotlib.pyplot as plt

class CountsPlotter():

    @staticmethod
    def plot_counts_table(counts_table, export_directory, sample_code):
        sorted_table = sorted(counts_table.items(),reverse=True)
        try:
            x, y = zip(*sorted_table)

        except ValueError as e:
            x = [50]
            y = [0]
        
        x_ticks_scaling_factor = (max(x)//40)+1


        graph = plt.bar(x,y, align='center', alpha=0.5)
        plt.xticks(list(range(1,max(x),x_ticks_scaling_factor)),
                   list(range(1,max(x),x_ticks_scaling_factor)),
                   fontsize=6 , rotation=30, fontweight='medium' ) #4.5
        

        plt.yticks(fontsize=6, fontweight='medium')

        plt.title(sample_code)
        plt.xlabel("Total number of repeat units" )
        plt.ylabel("Number of reads" )


        #plt.show()
        plt.savefig(export_directory, dpi=600)
        plt.clf()
