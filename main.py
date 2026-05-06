''' main module, run this file to start RGT

This module gets user inputs, assigns different processes to different files,
run them in parallel, gets the collective result from all samples, then
exports the summative result in an excel file.

this is done by creating multible instances of the RGT class that does the
analysis for each file

'''
import sys
import glob
from joblib import Parallel, delayed, cpu_count

from rgt import RGT
from interface.interface import get_user_inputs
from exporter.ExcelExport import ExcelWriter
from exporter.csv_exporter import export_csv


def get_collective_dictionary_from_list_of_output_dictionaries(list_of_output_dictionaries):
    ''' Convert the output of the parallel processing to a single dictionary''' 
    output_dictionary = {}
    for dictionary in list_of_output_dictionaries:
        key = list(dictionary.keys())[0] #get the first key(the only one)
        output_dictionary[key] = dictionary[key]
    return(output_dictionary)


def main():
    ''' gets user inputs, creats multible instances of the RGT class to analyse each sample
    and exports the summative result in an excel file
    '''
    input_directory, output_directory, settings, number_of_threads = get_user_inputs(sys.argv[1:])
    samples = glob.glob(input_directory + "/*.fastq")
    rgt_ = RGT(settings, input_directory, output_directory)

    if number_of_threads == None:
        number_of_threads = min(cpu_count(), len(samples))
    result = Parallel(n_jobs=number_of_threads, verbose=1)(map(delayed(rgt_.rgt),(samples)))
    
    automated_genotyope = [i[0] for i in result]
    color_table = [i[1] for i in result]
    
    output_dictionary = get_collective_dictionary_from_list_of_output_dictionaries(automated_genotyope)
    color_code_dictionary = get_collective_dictionary_from_list_of_output_dictionaries(color_table)
    collective_excel_writer = ExcelWriter()
    results_headers = ["sample ID", "First allele structure", "Second allele structure",
                        "Comments and Flags", "Identified peaks", "Discarded reads percentage %"]
    collective_excel_writer.add_table_to_sheet(output_dictionary,"results", results_headers,
                            color_table=color_code_dictionary, colored_cell_index=4)
    collective_excel_writer.save_file(output_directory + "/ResultsSummary.xlsx")

    if settings["additional_csv_export"]:
        export_csv(output_dictionary, output_directory + "/ResultsSummary.csv", header=results_headers)


if __name__== "__main__":
    main()