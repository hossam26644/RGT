import getopt
import sys, os, glob



def check_input_directory(input_directory):
    if not(os.path.isdir(input_directory)):
        print('invalid input directory')
        sys.exit()
    
    compressed = glob.glob(f"{input_directory}/*.fastq.gz")
    samples = compressed + [
        f for f in glob.glob(f"{input_directory}/*.fastq")
        if f + ".gz" not in compressed
    ]

    if len(samples)==0:
        print('no fastq files found in the input directory')
        sys.exit()

def check_number_of_threads(number_of_threads):
    try:
        int(number_of_threads)
    except Exception as e:
        print("Please enter a valid number of threads or leave blank") 
        sys.exit()

def check_or_create_output_directory(output_directory, settings):
    if not os.path.isdir(output_directory):
        try:
            os.mkdir(output_directory)
            os.mkdir(output_directory +"/Plots")
            print("output directory created")
        except:
            print("cannot acces or create output directory") 
            sys.exit()

    if not os.path.isdir(output_directory + "/Plots"):
        os.mkdir(output_directory +"/Plots")
        
    if not os.path.isdir(output_directory + "/Plots/total-counts-plots/"):
        os.mkdir(output_directory +"/Plots/total-counts-plots/")
    if not os.path.isdir(output_directory + "/Plots/specified-units-count-plots/"):
        os.mkdir(output_directory +"/Plots/specified-units-count-plots/")

    if not os.path.isdir(output_directory + "/FilesSpecificResults"):
            os.mkdir(output_directory +"/FilesSpecificResults")
    
    if settings["3D_plot_parameters"] is not None:
        if not os.path.isdir(output_directory + "/Plots/3d_plots/"):
            os.mkdir(output_directory +"/Plots/3d_plots/")        

    if settings["plot_waterfalls"]:
        os.mkdir(output_directory +"/Plots/waterfalls")


    if settings["additional_csv_export"]:
        if not os.path.isdir(output_directory + "/FilesSpecificResults/csv_exports"):
            os.mkdir(output_directory + "/FilesSpecificResults/csv_exports")
        if not os.path.isdir(output_directory + "/FilesSpecificResults/csv_exports/genotype_table"):
            os.mkdir(output_directory + "/FilesSpecificResults/csv_exports/genotype_table")
        if not os.path.isdir(output_directory + "/FilesSpecificResults/csv_exports/counts_table"):
            os.mkdir(output_directory + "/FilesSpecificResults/csv_exports/counts_table")
        if not os.path.isdir(output_directory + "/FilesSpecificResults/csv_exports/unique_counts_table"):
            os.mkdir(output_directory + "/FilesSpecificResults/csv_exports/unique_counts_table")
        if not os.path.isdir(output_directory + "/FilesSpecificResults/csv_exports/flanking_seq_table"):
            os.mkdir(output_directory + "/FilesSpecificResults/csv_exports/flanking_seq_table")
        if settings["match_singletons"]:
            if not os.path.isdir(output_directory + "/FilesSpecificResults/csv_exports/before_matching_singltons"):
                os.mkdir(output_directory + "/FilesSpecificResults/csv_exports/before_matching_singltons")
        


def check_settings_file(settings_file):
    if not os.path.exists(settings_file):
        print("settings file not be found") 
        sys.exit()
    if not settings_file.endswith('.json'):
        print("settings file should be in json format") 
        sys.exit()