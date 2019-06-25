from FileReader.ReadFile import ReadFile
from GenoTyper.GenoType import Genotype
from ExcelExporter.ExcelExport import ExcelWriter
from AllelesDetector.AllelesDetector import AllelesDetector
from CountsPlotter.CountsPlotter import CountsPlotter

import glob
from joblib import Parallel, delayed, cpu_count


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


        #read file and extract sequence from between flanks
        file = ReadFile(sample ,start_flank=self.settings["start_flank"],
                        end_flank=self.settings["end_flank"],
                        discard_reads_with_no_end_flank=self.settings["discard_reads_with_no_end_flank"])
        reads = file.reads #extracted reads from between flanks

        #genotype the reads (create counts table and repeat sequence abundance table)
        genotype = Genotype(reads,repeat_units=self.settings["repeat_units"],
                            unique_repeat_units=self.settings["unique_repeat_units"],
                            min_size_repeate=self.settings["min_size_repeate"],
                            max_interrupt_tract=self.settings["max_interrupt_tract"],
                            grouping_repeat_units=self.settings["grouping_repeat_units"])

        geno_table = genotype.get_geno_table() #the repeat sequence abundance table
        counts_table = genotype.get_counts_table() 
        unique_counts_table = genotype.get_unique_counts_table()

        #sort the three tables by abundance
        sorted_geno_table = dict(sorted(geno_table.items(), key=lambda x: x[1], reverse=True))
        sorted_counts_table = dict(sorted(counts_table.items(), key=lambda x: x[1], reverse=True))
        sorted_unique_counts_table = dict(sorted(unique_counts_table.items(), key=lambda x: x[1], reverse=True))

        #write three tabels to excel 
        excel_writer = ExcelWriter()
        geno_sheet_titles = ["sequence structure", "Abundance",
                            "Number of repeat units", "Number of unique repeat units", "Raw sequence structure"]
        counts_table_titles = ["Number of repeat units", "Abundance"]

        excel_writer.add_table_to_sheet(sorted_geno_table,"genotype", geno_sheet_titles)
        excel_writer.add_table_to_sheet(sorted_counts_table,"counts", counts_table_titles)
        excel_writer.add_table_to_sheet(sorted_unique_counts_table,"unique counts", counts_table_titles)
        excel_writer.save_file(self.output_directory + "/FilesSpecificResults/"+sample_code+".xlsx")

        #export plot
        table = genotype.get_counts_table()
        CountsPlotter.plot_counts_table(table, self.output_directory+ "/Plots/"+sample_code+".png" , sample_code)

        #Automaticly detect allels from counts table and geom table
        a = AllelesDetector(sorted_counts_table,sorted_geno_table)
        output_table[sample_code] = a.get_alleles()
        output_table[sample_code].append(file.get_discarded_reads_percentage())

        color_table[sample_code] = a.color_code

        return[output_table, color_table]

