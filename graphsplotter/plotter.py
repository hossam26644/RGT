'''docstring'''
from .2d_table_plotter.py import plot_2d_table

def plot_graphs(genotype,output_directory, sample_code, first_allele, second_allele, color_code):
    
    table = genotype.get_counts_table()
    plot_directory = output_directory+ "/Plots/"+sample_code+".png"
        
    plot_2d_table(table, plot_directory, sample_code,
        first_allele, second_allele, color_code=color_code)
