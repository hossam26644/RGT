from . import check_params 
import getopt
import sys, os

def get_user_inputs(argv):
    input_directory = ''
    output_directory = ''
    settings_file = ''
    number_of_threads = None
    try:
        opts, args = getopt.getopt(argv,"hi:o:s:t:", ["help", "input", "output", "sys", "threads"])
    except getopt.GetoptError:
        print ('test.py -i <input directory> -o <output directory> -s <settings json file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h','--help'):
            print('test.py -i <input directory> -o <output directory> -s <settings json file> ')
            sys.exit()
        elif opt in ("-i"):
            input_directory = arg
        elif opt in ("-o"):
            output_directory = arg
        elif opt in ("-s"):
            settings_file = arg
        elif opt in ("-t"):
            check_params.check_number_of_threads(arg)
            number_of_threads = int(arg)

    check_params.check_input_directory(input_directory)
    check_params.check_or_create_output_directory(output_directory)
    check_params.check_settings_file(settings_file)
    

    return(input_directory, output_directory, settings_file, number_of_threads)
