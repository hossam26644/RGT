''' RGT class, does analysis on a single sample (fastq file)
'''
import glob
from joblib import Parallel, delayed, cpu_count

from filereader.ReadFile import ReadFile
from genotyper.GenoType import Genotype, get_rev_complementry
from excelexporter.ExcelExport import ExcelWriter
from allelesdetector.AllelesDetector import AllelesDetector
from graphsplotter.plotter import plot_graphs




class RGT():
    
    def __init__(self, settings, input_directory, output_directory):
        self.settings = settings
        self.input_directory = input_directory
        self.output_directory = output_directory

    def rgt(self, sample):

        output_table = {} #dictionary to export result
        color_table = {}
        #sample code chopping from file + directory name
        sample_code = sample.split("/")[-1]  #gets the last item after backslashes
        sample_code = sample_code.split(".")[0] #gets file name without extension
        print(sample_code)

        try:
            #read file and extract sequence from between flanks
            file = ReadFile(sample ,start_flank=self.settings["start_flank"],
                            end_flank=self.settings["end_flank"],
                            discard_reads_with_no_end_flank=self.settings["discard_reads_with_no_end_flank"])
            reads = file.reads #extracted reads from between flanks
            #genotype the reads (create counts table and repeat sequence abundance table)
            genotype = Genotype(reads,self.settings)

            geno_table = genotype.get_geno_table() #the repeat sequence abundance table
            counts_table = genotype.get_counts_table() 
            unique_counts_table = genotype.get_unique_counts_table()

            #sort the three tables by abundance
            sorted_geno_table = dict(sorted(geno_table.items(), key=lambda x: x[1], reverse=True))
            sorted_counts_table = dict(sorted(counts_table.items(), key=lambda x: x[1], reverse=True))
            sorted_unique_counts_table = dict(sorted(unique_counts_table.items(), key=lambda x: x[1], reverse=True))

            #write three tabels to excel 
            excel_writer = ExcelWriter()
            if self.settings["3D_plot_parameters"] != None:
                xlabel =' , '.join(self.settings["3D_plot_parameters"]["x_units"]) + " count"
                zlabel =' , '.join(self.settings["3D_plot_parameters"]["z_units"]) + " count"
            else:
                xlabel = "x axis units count for 3D plots, not applicable"
                zlabel = "z axis units count for 3D plots, not applicable"
            geno_sheet_titles = ["sequence structure", "Abundance",
                                "Number of repeat units",xlabel , zlabel, "Number of unique repeat units", "Raw sequence structure"]
            counts_table_titles = ["Number of repeat units", "Abundance"]

            excel_writer.add_table_to_sheet(sorted_geno_table,"genotype", geno_sheet_titles)
            excel_writer.add_table_to_sheet(sorted_counts_table,"counts", counts_table_titles)
            excel_writer.add_table_to_sheet(sorted_unique_counts_table,"unique counts", counts_table_titles)
            excel_writer.save_file(self.output_directory + "/FilesSpecificResults/"+sample_code+".xlsx")

            #Automaticly detect allels from counts table and geom table
            a = AllelesDetector(sorted_counts_table, sorted_geno_table, self.settings["minimum_no_of_reads"])
            output_table[sample_code] = a.result_summery
            discarded_reads_percentage = file.get_discarded_reads_percentage()
            output_table[sample_code].append(str(round(discarded_reads_percentage, 1)))

            #export plot
            plot_graphs(self.settings, genotype,self.output_directory, sample_code, a.first_allele, a.second_allele, a.color_code)

            color_table[sample_code] = {4:a.color_code}
            self.color_code_discarded_reads_percntg(color_table, discarded_reads_percentage,sample_code)

        except Exception as e:
            #print(e)
            print("can not genotype " + sample_code)
            output_table[sample_code] = ["can not get allele","can not get allele","Error","[]","!"]
            color_table[sample_code] = {1:"red", 4:"red", 6:"red"}
        return[output_table, color_table]

    def color_code_discarded_reads_percntg(self, color_table, prcntg,sample_code):
        if prcntg >= int(self.settings["discarded_reads_flag_percentage"]):
            this_sample_dict = color_table[sample_code]
            this_sample_dict[6] = "red"

