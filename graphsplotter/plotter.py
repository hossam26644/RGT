'''docstring'''
from .table_2d_plotter import plot_2d_table
from .plot_3D import plot_3D

def plot_graphs(settings, genotype,output_directory, sample_code, first_allele, second_allele, color_code):
    
    #plot counts table
    total_counts_table = genotype.get_counts_table()
    plot_directory = output_directory+ "/Plots/total-counts-plots/"+sample_code+".png"
    first_allele_count = first_allele.repeat_units_count
    second_allele_count = second_allele.repeat_units_count
    
    plot_2d_table(total_counts_table, plot_directory, sample_code,
        first_allele, second_allele,
        first_allele_count, second_allele_count,"Total number of repeat units", color_code=color_code)

    #plot specific units count
    specified_units_counts_table = genotype.get_unique_counts_table()
    plot_directory = output_directory+ "/Plots/specified-units-count-plots/"+sample_code+".png"
    first_allele_count = first_allele.unique_units_count
    second_allele_count = second_allele.unique_units_count
    xlabel = ' , '.join(settings["unique_repeat_units"])
    xlabel += " count"
    plot_2d_table(specified_units_counts_table, plot_directory, sample_code,
        first_allele, second_allele,
        first_allele_count, second_allele_count,xlabel, color_code=color_code)

    #plot the 3D plot
    if settings["3D_plot_parameters"] != None:
        table_3d = genotype.table_3d
        plot_directory = output_directory+ "/Plots/3d_plots/"+sample_code+".png"

        xlabel =' , '.join(settings["3D_plot_parameters"]["x_units"]) + " count"
        zlabel =' , '.join(settings["3D_plot_parameters"]["z_units"]) + " count"

        if xlabel == zlabel:
            xlabel += "(1)"
            zlabel += "(2)"
        plot_3D(table_3d, plot_directory, sample_code, first_allele, second_allele,
                xlabel, zlabel, color_code=color_code)
