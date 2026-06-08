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

def check_3d_plot_params(params_3d):
    try:
        params_3d["x_units"]
        params_3d["z_units"]
    except:
        print("3D parameters are not complete, either remove 3D_plot_parameters attripute or completely "
         "insert proper list of x_units and z_units")
        sys.exit()
    keys = ["before_x_seq", "after_x_seq", "before_z_seq", "after_z_seq"]
    for key in keys:
        try:
            params_3d[key]
            if (params_3d[key] == "" or params_3d[key] == [] ):
                raise ValueError("null parameter")
        except:
            params_3d[key] = None

def check_parameters(settings):
    try:
        settings["repeat_units"]
        settings["repeat_units"].sort(key=len,reverse=True)
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
        settings["unique_repeat_units"].sort(key=len,reverse=True)
    except Exception as e:
        settings["unique_repeat_units"] = settings["repeat_units"]
   
    try:
        settings["grouping_repeat_units"]
        settings["grouping_repeat_units"].sort(key=len,reverse=True)
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
        if settings["discard_reads_with_no_end_flank"] != "smart":
            settings["discard_reads_with_no_end_flank"] = get_bool_value_from_string(settings["discard_reads_with_no_end_flank"])
    except Exception as e:
        settings["discard_reads_with_no_end_flank"] = "smart"
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
        settings["PCR_free"]
        settings["PCR_free"] = get_bool_value_from_string(settings["PCR_free"])
    except Exception as e:
        settings["PCR_free"] = False

    try:
    	settings["minimum_no_of_reads"]
    except Exception as e:
        if settings["PCR_free"] == True:
            settings["minimum_no_of_reads"] = 1
        else:
            settings["minimum_no_of_reads"] = 30
    try:
        settings["number_of_allowed_strt_flank_point_mutations"]
    except Exception as e:
        settings["number_of_allowed_strt_flank_point_mutations"] = 1
    try:
        settings["number_of_allowed_end_flank_point_mutations"]
    except Exception as e:
        settings["number_of_allowed_end_flank_point_mutations"] = 1

    try:
        settings["3D_plot_parameters"]
        check_3d_plot_params(settings["3D_plot_parameters"])
    except Exception as e:
        settings["3D_plot_parameters"] = None

    try:
        settings["additional_csv_export"]
        settings["additional_csv_export"] = get_bool_value_from_string(settings["additional_csv_export"])
    except Exception as e:
        settings["additional_csv_export"] = False
    
    try:
        settings["match_singltons"]
    except Exception as e:
        settings["match_singltons"] = 1

    try:
        int(settings["report_consensus_flanking_sequence"])
    except Exception as e:
        settings["report_consensus_flanking_sequence"] = False
    
    print(settings)