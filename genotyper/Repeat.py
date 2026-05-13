'''docstring'''
from .SmartString import SmartString
from .GroupingString import GroupingString
from .revComplementry import get_rev_complementry

class Read():
    def __init__(self, header=None, raw_read=None, quality_score=None):
        self.header = header
        self.raw_read = raw_read
        self.quality_score = quality_score
        self.interflanking_seq = None
        self.start_flank_found = False
        self.end_flank_found = False
        self.successfully_extracted = False
        self.repeats = []
    
    def state(self):
        return {
            "header": self.header,
            "raw_read": self.raw_read,
            "quality_score": self.quality_score,
            "interflanking_seq": self.interflanking_seq,
            "start_flank_found": self.start_flank_found,
            "end_flank_found": self.end_flank_found,
            "successfully_extracted": self.successfully_extracted,
        }
        
class Repeat(Read):
    """docstring for Genotype"""
    def __init__(self, read):

        for k, v in read.state().items():
            setattr(self, k, v)

        self.repeat_units = None
        self.number_of_units = 0
        self.last_unit_index = None
        self.start_index = None 
        self.single_point_mutation_indexes = []
        self.unconfirmed_units_buffer = 0
        self.non_perfect_units = 0
        self.repeat_sequence = ""
        self.unconfirmed_sequence = "" #sequence of non pure repeates waiting for a confirmed unit to  be added
        self.unique_repeat_units_list = None
        self.unique_repeat_units_count = 0
        self.plot_3d_settings = None
        #for 3D plots
        if self.plot_3d_settings is not None:
            self.x_units = plot_3d_settings["x_units"]
            self.z_units = plot_3d_settings["z_units"]
            self.before_x_seq = plot_3d_settings["before_x_seq"]
            self.after_x_seq = plot_3d_settings["after_x_seq"]
            self.before_z_seq = plot_3d_settings["before_z_seq"]
            self.after_z_seq = plot_3d_settings["after_z_seq"]
            self.x_counts_for_3d = 0
            self.z_counts_for_3d = 0
            self.in_x_region_flag = False
            self.in_z_region_flag = False
            if self.before_x_seq == None:
                self.in_x_region_flag = True

            if self.before_z_seq == None:
                self.in_z_region_flag = True


    def add_unit(self, window, index): #index is the last base index in the sequence
        if window in self.repeat_units: #perfect match            
            self.number_of_units +=1
            self.number_of_units += self.unconfirmed_units_buffer #add the unconfirmed units (points with SNPs)
            self.non_perfect_units += self.unconfirmed_units_buffer
            self.unconfirmed_units_buffer = 0 #zero the points with SNPs
            self.change_last_unit_index(index)
            self.repeat_sequence += self.unconfirmed_sequence
            self.unconfirmed_sequence = ""
            self.repeat_sequence += window

        else:
            self.unconfirmed_units_buffer += 1
            #self.single_point_mutation_indexes.append(self.get_SNP_index(window, index))
            self.unconfirmed_sequence += window

        if window in self.unique_repeat_units_list:
            self.unique_repeat_units_count +=1
        
        if self.plot_3d_settings!= None:
            self.add_x_count(index, window)
            self.add_z_count(index, window)
      
    def add_z_count(self, index, window):
        if self.check_seq_before_repeat(index-len(window), self.before_z_seq):
            self.in_z_region_flag = True
        
        #print(window in self.z_units, window, self.z_units)
        if self.in_z_region_flag and (window in self.z_units):
            self.z_counts_for_3d += 1

        if self.check_seq_after_repeat(index, self.after_z_seq):
            self.in_z_region_flag = False  

    def add_x_count(self, index, window):
        if self.check_seq_before_repeat(index-len(window), self.before_x_seq):
            self.in_x_region_flag = True

        if self.in_x_region_flag and (window in self.x_units):
            self.x_counts_for_3d += 1

        if self.check_seq_after_repeat(index, self.after_x_seq):
            self.interflanking_seq
            self.in_x_region_flag = False

    def check_seq_before_repeat(self, window_start_index, before_seq):
        if before_seq != None:
            if self.interflanking_seq[window_start_index-len(before_seq):window_start_index] == before_seq:
                return True
            else:
                return False 
        else:
            return False
    
    def check_seq_after_repeat(self, index, after_seq):
        if after_seq != None:
            if self.interflanking_seq[index:index+len(after_seq)] == after_seq:
                return True
            else:
                return False
        else:
            return False

    def change_last_unit_index(self, index):
        self.last_unit_index = index

    def get_non_perfect_units_percentage(self):
    	return self.non_perfect_units/self.number_of_units

    def get_seq(self):
    	return self.interflanking_seq[self.start_index:self.last_unit_index]
    
    def get_seq_smart_string(self, rev_strand):
        if rev_strand:
            reversed_seq = get_rev_complementry(self.get_seq())
            smart_string = SmartString.get_smart_string_from_sequence(reversed_seq, 3, self.repeat_units)
        else:
            smart_string = SmartString.get_smart_string_from_sequence(self.get_seq(), 3, self.repeat_units)
        return smart_string

    def get_grouped_string(self, grouping_repeat_units, rev_strand):
        if rev_strand:
            reversed_seq = get_rev_complementry(self.get_seq())
            grouping_string = GroupingString.get_grouped_string_from_sequence(reversed_seq, grouping_repeat_units)
        else:
            grouping_string = GroupingString.get_grouped_string_from_sequence(self.get_seq(), grouping_repeat_units)
        return grouping_string
        
    def __eq__(self, other):
        if isinstance(other, Repeat):
            return self.get_seq() == other.get_seq()
        return False