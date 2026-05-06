"""
Exporting a tavle to a csv file 
"""
import csv


def export_csv(table:dict, file_name:str, header:list=""):

    """ table is a dictionary with keys as the first column and values the other columns"""
    with open(file_name, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)


        for key in sorted(table.keys()):
            values = table[key]
            if not isinstance(values, list):
                values = [values]
            writer.writerow([key] + values)
            
        