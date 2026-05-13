''' RGT class, does analysis on a single sample (fastq file)
'''
import glob
from joblib import Parallel, delayed, cpu_count

from filereader.ReadFile import ReadFile
from genotyper.GenoType import Genotype, get_rev_complementry
from exporter.ExcelExport import ExcelWriter
from exporter.csv_exporter import export_csv
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
            file = ReadFile(sample, self.settings)

            #genotype the reads (create counts table and repeat sequence abundance table)
            genotype = Genotype(file.reads, self.settings)
            ad = AllelesDetector(genotype, self.settings)

            if self.settings['discard_reads_with_no_end_flank'] == 'smart':
                # reinclude reads with no end flank if the have more repeats than the longer allele
                max_detected_expanded = 0
                for read in file.reads:
                    for repeat in read.repeats:
                        if not read.end_flank_found and repeat.number_of_units and repeat.number_of_units > ad.longer_allele_repeats_count:
                            read.successfully_extracted = True
                            max_detected_expanded = max(max_detected_expanded, repeat.number_of_units)
                if max_detected_expanded:
                    genotype = Genotype(file.reads, self.settings)
                    ad = AllelesDetector(genotype, self.settings)
                    ad.color_code = "red"
                    ad.result_summery[2] += f" ,Expanded allele detected with no end flank"
                    ad.result_summery[2] += f" ,max detected expanded allele has {max_detected_expanded} repeat units, please check manually"                          
                        
            geno_table = genotype.get_geno_table() #the repeat sequence abundance table
            counts_table = genotype.get_counts_table() 
            unique_counts_table = genotype.get_unique_counts_table()

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

            excel_writer.add_table_to_sheet(geno_table,"genotype", geno_sheet_titles)
            if self.settings["match_singltons"]:
                excel_writer.add_table_to_sheet(genotype.before_matching_table,"before_matching_singltons", geno_sheet_titles)
            excel_writer.add_table_to_sheet(counts_table,"counts", counts_table_titles)
            excel_writer.add_table_to_sheet(unique_counts_table,"unique counts", counts_table_titles)
            excel_writer.save_file(self.output_directory + "/FilesSpecificResults/"+sample_code+".xlsx")

            if self.settings["additional_csv_export"]:
                csv_geno_table_path = f"{self.output_directory}/FilesSpecificResults/csv_exports/genotype_table/{sample_code}.csv"
                csv_counts_table_path = f"{self.output_directory}/FilesSpecificResults/csv_exports/counts_table/{sample_code}.csv"
                csv_unique_counts_table_path = f"{self.output_directory}/FilesSpecificResults/csv_exports/unique_counts_table/{sample_code}.csv"
                export_csv(geno_table, csv_geno_table_path, header=geno_sheet_titles)
                export_csv(counts_table, csv_counts_table_path, header=counts_table_titles)
                export_csv(unique_counts_table, csv_unique_counts_table_path, header=counts_table_titles)
                if self.settings["match_singltons"]:
                    csv_before_matching_table_path = f"{self.output_directory}/FilesSpecificResults/csv_exports/before_matching_singltons/{sample_code}.csv"
                    export_csv(genotype.before_matching_table, csv_before_matching_table_path, header=geno_sheet_titles)


            #Automaticly detect allels from counts table and geom table
            output_table[sample_code] = ad.result_summery
            discarded_reads_percentage = file.get_discarded_reads_percentage()
            output_table[sample_code].append(str(round(discarded_reads_percentage, 1)))

            #export plot
            plot_graphs(self.settings, genotype, self.output_directory, sample_code, ad.first_allele, ad.second_allele, ad.color_code)

            color_table[sample_code] = {4:ad.color_code}
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

