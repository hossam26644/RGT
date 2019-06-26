import json
import sys
import distutils

def extract_parameters(json_file):
    with open(json_file) as json_file:  
        settings = json.load(json_file)
        check_parameters(settings)
        return(settings)

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
        settings["discard_reads_with_no_end_flank"] = distutils.util.strtobool((settings["discard_reads_with_no_end_flank"]))
    except Exception as e:
        settings["discard_reads_with_no_end_flank"] = True
    try:
        int(settings["discarded_reads_flag_percentage"])
    except Exception as e:
        settings["discarded_reads_flag_percentage"] = 10