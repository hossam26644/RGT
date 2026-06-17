'''docstring'''
from .Repeat import Repeat
from .revComplementry import get_rev_complementry
import copy

class Genotype():
    """docstring for Genotype"""
    def __init__(self, reads, settings):
        self.settings =  copy.deepcopy(settings)
        self.reverse_strand = settings["reverse_strand"]
        if self.reverse_strand:
            self.repeat_units = get_rev_complementry(self.settings["repeat_units"])
            self.unique_repeat_units = get_rev_complementry(self.settings["unique_repeat_units"])
           
            if self.settings["3D_plot_parameters"]!= None:
                for key in self.settings["3D_plot_parameters"]:
                    try:
                        self.settings["3D_plot_parameters"][key] = get_rev_complementry(
                            self.settings["3D_plot_parameters"][key])
                    except Exception as e:
                        pass
                temp = self.settings["3D_plot_parameters"]["before_x_seq"]
                self.settings["3D_plot_parameters"]["before_x_seq"] = self.settings["3D_plot_parameters"]["after_x_seq"]

                self.settings["3D_plot_parameters"]["after_x_seq"] = temp
                
                temp = self.settings["3D_plot_parameters"]["before_z_seq"]
                self.settings["3D_plot_parameters"]["before_z_seq"] = self.settings["3D_plot_parameters"]["after_z_seq"]

                self.settings["3D_plot_parameters"]["after_z_seq"] = temp


        else:
            self.repeat_units = self.settings["repeat_units"]
            self.unique_repeat_units = self.settings["unique_repeat_units"]
            
        self.list_of_repeat_units_lengths = self.get_list_of_repeat_units_lengths()#used to have diff sliding windows lengths
        self.reads = reads
        self.min_size_repeat = self.settings["min_size_repeat"]
        self.max_interrupt_tract = self.settings["max_interrupt_tract"]
        self.grouping_repeat_units = self.settings["grouping_repeat_units"]
        self.plot_3d_settings = self.settings["3D_plot_parameters"]
        self.geno_table = {}
        self.counts_table = {}
        self.unique_counts_table = {}
        self.table_3d = {}
        if self.settings["match_singletons"]: self.before_matching_table = {}
        self.repeats = []
        self.genotype_repeats()

    def clear_tables(self):
        self.geno_table = {}
        self.counts_table = {}
        self.unique_counts_table = {}
        self.table_3d = {}
        if self.settings["match_singletons"]: self.before_matching_table = {}

    def genotype_repeats(self):
        genotyped_repeats = []
        for read in self.reads:
            if read.interflanking_seq is not None:
                genotyped_repeats.extend(self.get_line_repeates(read))
                   

        if self.settings["match_singletons"]:
            filtered_repeats = self.match_singletons(genotyped_repeats)         
            self.repeats = filtered_repeats
        else:
            self.repeats = genotyped_repeats
        
        for repeat in self.repeats: 
            if repeat.get_non_perfect_units_percentage() <= 0.3: #only add repeates with unique percentage > 0.3
                repeat.discarded = False
                self.add_repeat_to_tables(repeat)        

    def match_singletons(self, genotyped_repeats):
        self.clear_tables()
        for repeat in genotyped_repeats:
            self.add_repeat_to_genotable(repeat)
        geno_table = self.get_geno_table()
        
        good_sequences = {k: v for k, v in geno_table.items() if v[0] > 1}
        good_sequences_sorted = [
            v[5] for k, v in sorted(good_sequences.items(), key=lambda item: item[0], reverse=True)
        ]
        # reverse lookup: sequence_str -> key
        seq_to_key = {v[5]: k for k, v in geno_table.items()}

        filtered_matched_repeats = []
        for repeat in genotyped_repeats:
            abuncance = geno_table.get(repeat.get_seq_smart_string(self.reverse_strand))
            if repeat.get_seq() in good_sequences_sorted:
                filtered_matched_repeats.append(repeat)
                continue
            elif abuncance is not None and abuncance[0] > 1:
                filtered_matched_repeats.append(repeat)
                continue
            else:
                seq = next(
                    (seq for n in range(self.settings["match_singletons"])
                        for seq in good_sequences_sorted
                        if Repeat.levenshtein(seq, repeat.get_seq()) == n + 1),
                    None
                )

                aligned_to_repeat = next(
                    (r for r in genotyped_repeats if r.get_seq() == seq),
                    None
                )
                if aligned_to_repeat is None:
                    filtered_matched_repeats.append(repeat)
                else:
                    filtered_matched_repeats.append(aligned_to_repeat)
                
        
        before_matching_table = copy.deepcopy(self.geno_table)     
        self.clear_tables()
        self.before_matching_table = before_matching_table
        
        return filtered_matched_repeats

    def get_line_repeates(self, read):
        genotyped_repeats = []
        i = 0
        repeat = None
        while i <= len(read.interflanking_seq)-min(self.list_of_repeat_units_lengths): #sliding window
            checker = self.window_enters_repeat_sequence(i, self.repeat_units, repeat, read.interflanking_seq)
            if checker[0]:
                '''if window detects a repeat unit, while it is not inside a repeat sequence'''
                repeat = Repeat(read, self.settings["3D_plot_parameters"])
                window = checker[1]
                repeat.start_index = i
                repeat.last_unit_index = i+len(window)
                repeat.repeat_units = self.repeat_units
                repeat.unique_repeat_units_list = self.unique_repeat_units
                repeat.plot_3d_settings = self.plot_3d_settings
                repeat.add_unit(window,i+len(window)) #add a repeat unit count

                i = i+len(window) #Jumb one window
                continue
            
            checker = self.detect_repeat_unit_inside_repeat(i, self.repeat_units, repeat, read.interflanking_seq)
            if checker[0]:
                '''if it detects a repeat while inside the repeat sequence'''
                window = checker[1]
                repeat.add_unit(window,i+len(window)) #add a repeat unit count
                i = i+len(window) #jumb one window
                continue

            elif self.non_matching_unit_within_repeat(i, self.repeat_units, repeat, read.interflanking_seq):
                #print(window,i,repeat.last_unit_index, i-repeat.last_unit_index)
                if i-repeat.last_unit_index < self.max_interrupt_tract:
                    #ignore if length is smaller than max interrupt tract
                    i += 1
                    continue
                #if length is larger than max interrupt tract
                elif repeat.number_of_units >= self.min_size_repeat: #check that number of repeates is larger than the minimum size repeate
                    read.repeats.append(repeat)
                    if read.successfully_extracted:
                        genotyped_repeats.append(repeat)
                repeat = None
            i +=1 

        if repeat != None and repeat.number_of_units >= self.min_size_repeat: #if sequence ends on a repeat
            read.repeats.append(repeat)
            if read.successfully_extracted:
                genotyped_repeats.append(repeat)

        return genotyped_repeats

    def non_matching_unit_within_repeat(self,idx, repeat_units,  repeat_object,sequence):
        window_inside_repeates_flag = repeat_object != None
        if window_inside_repeates_flag:
            if self.do_repeat_unit_exist(idx,repeat_object,sequence)[0]:
                return False
            return True  
        return False

    def detect_repeat_unit_inside_repeat(self,idx, repeat_units, repeat_object, sequence):
        window_inside_repeates_flag = repeat_object != None
        if window_inside_repeates_flag:
            checker = self.do_repeat_unit_exist(idx,repeat_object,sequence)
            if checker[0]:
                return True, checker[1]
        return False, ""

    def window_enters_repeat_sequence(self,idx, repeat_units, repeat_object,sequence): 
        window_inside_repeates_flag = repeat_object != None
        if not window_inside_repeates_flag:
            checker = self.do_repeat_unit_exist(idx,repeat_object,sequence)
            if checker[0]:
                return True, checker[1]
        return False, ""

    def do_repeat_unit_exist(self,idx,repeat_object,sequence):
        similar_seq = ""
        for length in self.list_of_repeat_units_lengths:
            window = sequence[idx:idx+length]
            for repeat_unit in self.repeat_units:
                if window == repeat_unit:
                    return True, window
                elif repeat_object != None and Repeat.hamming_distance(window, repeat_unit)==1:
                    if similar_seq == "":
                        similar_seq = window
        if similar_seq != "":
            return True, similar_seq
        return False, ""

    def add_repeat_to_tables(self, repeat):
        self.add_repeat_to_genotable(repeat)
        self.add_repeat_to_countstable(repeat)
        self.add_repeat_to_unique_countstable(repeat)
        if self.settings["3D_plot_parameters"]!= None:
            self.add_repeat_to_table3d(repeat)

    def add_repeat_to_table3d(self, repeat):
        key = (repeat.x_counts_for_3d, repeat.z_counts_for_3d)
        if key in self.table_3d:
            self.table_3d[key] += 1
        else:
            self.table_3d[key] = 1

    def add_repeat_to_genotable(self, repeat):
        if self.grouping_repeat_units == None:
            repeat_sequence = repeat.get_seq_smart_string(self.reverse_strand)
        else:
            repeat_sequence = repeat.get_grouped_string(self.grouping_repeat_units,self.reverse_strand)
        
        if(repeat_sequence in self.geno_table):
            self.geno_table[repeat_sequence][0] += 1
        else:
            if self.settings["3D_plot_parameters"]!= None:
                self.geno_table[repeat_sequence] = [1,repeat.number_of_units,
                    repeat.x_counts_for_3d, repeat.z_counts_for_3d,
                    repeat.unique_repeat_units_count,repeat.get_seq()]
            else:
                self.geno_table[repeat_sequence] = [1,repeat.number_of_units,
                    "not applicable", "not applicable",
                    repeat.unique_repeat_units_count,repeat.get_seq()]

    def add_repeat_to_countstable(self, repeat):
        number_of_repeat_units = repeat.number_of_units
        
        if(number_of_repeat_units in self.counts_table):
            self.counts_table[number_of_repeat_units] += 1
        else:
            self.counts_table[number_of_repeat_units] = 1

    def add_repeat_to_unique_countstable(self, repeat):
        number_of_unique_repeat_units = repeat.unique_repeat_units_count
        
        if(number_of_unique_repeat_units in self.unique_counts_table):
            self.unique_counts_table[number_of_unique_repeat_units] += 1
        else:
            self.unique_counts_table[number_of_unique_repeat_units] = 1

    def get_geno_table(self):
        return dict(sorted(self.geno_table.items(), key=lambda x: x[1], reverse=True))        

    def get_counts_table(self):
        return dict(sorted(self.counts_table.items(), key=lambda x: x[1], reverse=True))

    def get_unique_counts_table(self):
        return dict(sorted(self.unique_counts_table.items(), key=lambda x: x[1], reverse=True))

    def get_list_of_repeat_units_lengths(self):
        list_of_repeat_units_lengths =[]
        for unit in self.repeat_units:
            list_of_repeat_units_lengths.append(len(unit))
        list_of_repeat_units_lengths = set(list_of_repeat_units_lengths) #set removes repeats
        list_of_repeat_units_lengths = list(list_of_repeat_units_lengths) #making lenghts a list
        list_of_repeat_units_lengths.sort(reverse = True) #sort descendingly
       
        return(list_of_repeat_units_lengths)

    def get_discarded_reads_percentage(self):
        percentage = len(self.reads) - sum([v[0] for k,v in self.geno_table.items()])
        percentage /= len(self.reads)
        percentage *= 100
        return percentage  
