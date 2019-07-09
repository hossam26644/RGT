import getopt
import sys, os, glob



def check_input_directory(input_directory):
    if not(os.path.isdir(input_directory)):
        print('invalid input directory')
        sys.exit()
    
    samples = glob.glob(input_directory + "/*.fastq")
    if len(samples)==0:
        print('no fastq files found in the input directory')
        sys.exit()

def check_number_of_threads(number_of_threads):
    try:
        int(number_of_threads)
    except Exception as e:
        print("Please enter a valid number of threads or leave blank") 
        sys.exit()

def check_or_create_output_directory(output_directory):
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
        
    if not (os.path.isdir(output_directory + "/Plots/total-counts-plots/") and
             os.path.isdir(output_directory + "/Plots/specified-units-count-plots/")):
            
            os.mkdir(output_directory +"/Plots/total-counts-plots/")
            os.mkdir(output_directory +"/Plots/specified-units-count-plots/")
            print("Plots directory created inside output directory")

    if not os.path.isdir(output_directory + "/FilesSpecificResults"):
            os.mkdir(output_directory +"/FilesSpecificResults")
            print("FilesSpecificResults created inside output directory")

def check_settings_file(settings_file):
    if not os.path.exists(settings_file):
        print("settings file not be found") 
        sys.exit()
    if not settings_file.endswith('.json'):
        print("settings file should be in json format") 
        sys.exit()