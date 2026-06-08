''' main module, run this file to start RGT

This module gets user inputs, assigns different processes to different files,
run them in parallel, gets the collective result from all samples, then
exports the summative result in an excel file.

this is done by creating multible instances of the RGT class that does the
analysis for each file

'''
import sys
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed
from os import cpu_count

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
    compressed = glob.glob(f"{input_directory}/*.fastq.gz")
    samples = compressed + [
        f for f in glob.glob(f"{input_directory}/*.fastq")
        if f + ".gz" not in compressed
    ]

    rgt_ = RGT(settings, input_directory, output_directory)

    if number_of_threads == None:
        number_of_threads = min(cpu_count(), len(samples))
    else:
        number_of_threads = min(number_of_threads, len(samples))
   
    with ProcessPoolExecutor(max_workers=number_of_threads) as executor:
        futures = {executor.submit(rgt_.rgt, s): i for i, s in enumerate(samples)}
        result = [None] * len(samples)
        for done_count, future in enumerate(as_completed(futures), 1):
            idx = futures[future]
            result[idx] = future.result()
            if done_count % 10 == 0 or done_count == len(samples):
                print(f"[{done_count}/{len(samples)}] samples processed")

    automated_genotyope = [i[0] for i in result]
    color_table = [i[1] for i in result]
    
    output_dictionary = get_collective_dictionary_from_list_of_output_dictionaries(automated_genotyope)
    color_code_dictionary = get_collective_dictionary_from_list_of_output_dictionaries(color_table)
    collective_excel_writer = ExcelWriter()
    
    if settings["report_consensus_flanking_sequence"]: header_flank_text = [f"consensus {f} flanks structure" for f in ["start","end"]]
    else: header_flank_text = [f"most common {f} flank" for f in ["start","end"]]
    
    results_headers = ["sample ID", "First allele structure", "Second allele structure",
                        "Comments and Flags", "Identified peaks",
                        f"Allele 1 {header_flank_text[0]}", f"Allele 1 {header_flank_text[1]}",
                        f"Allele 2 {header_flank_text[0]}", f"Allele 2 {header_flank_text[1]}",
                        "Reads with no identified flanks %", "Discarded reads percentage %"]
    
    collective_excel_writer.add_table_to_sheet(output_dictionary,"results", results_headers,
                            color_table=color_code_dictionary, colored_cell_index=4, align_right_index=[6,8],
                            color_odd_element_index=[6,7,8,9])

    collective_excel_writer.save_file(output_directory + "/ResultsSummary.xlsx")

    if settings["additional_csv_export"]:
        export_csv(output_dictionary, output_directory + "/ResultsSummary.csv", header=results_headers)


if __name__== "__main__":
    main()