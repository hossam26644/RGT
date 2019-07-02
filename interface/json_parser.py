import json
import sys
import distutils

def extract_parameters(json_file):
    with open(json_file) as json_file:  
        settings = json.load(json_file)
        check_parameters(settings)
        return(settings)

def get_bool_value_from_string(given_string):
    if given_string in [0,"0","False","FALSE","false","no","No","NO","n","N" ]:
        return False
    elif given_string in [1,"1","True","TRUE","true", "yes","Yes","YES","y","Y"]:
        return True
    else:
        raise ValueError("given input can not be determined")

def check_parameters(settings):
    try:
        settings["repeat_units"]
    except Exception as e:
        print("list of repeat units should be provided")
        sys.exit()

    try:
        settings["start_flank"]
        settings["end_flank"]
    except Exception as e:
        print("Warning: no flanking sequence is selected")
        settings["start_flank"] = None
        settings ["end_flank"] = None

    try:
        settings["unique_repeat_units"]
    except Exception as e:
        settings["unique_repeat_units"] = settings["repeat_units"]
   
    try:
        settings["grouping_repeat_units"]
    except Exception as e:
        settings["grouping_repeat_units"] = None
    
    try:
        settings["min_size_repeate"]
    except Exception as e:
        settings["min_size_repeate"] = 5

    try:
        settings["max_interrupt_tract"]
    except Exception as e:
        settings["max_interrupt_tract"] = 5

    try:
        settings["discard_reads_with_no_end_flank"]
        settings["discard_reads_with_no_end_flank"] = get_bool_value_from_string(settings["discard_reads_with_no_end_flank"])
        
    except Exception as e:
        settings["discard_reads_with_no_end_flank"] = False
    try:
        int(settings["discarded_reads_flag_percentage"])
    except Exception as e:
        settings["discarded_reads_flag_percentage"] = 10
    try:
        settings["reverse_strand"]
        settings["reverse_strand"] = get_bool_value_from_string(settings["reverse_strand"])
    except Exception as e:
        settings["reverse_strand"] = False
    try:
    	settings["minimum_no_of_reads"]
    except Exception as e:
    	settings["minimum_no_of_reads"] = 30
    print(settings)