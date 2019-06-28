'''docstring'''
from .SmartString import SmartString
from .GroupingString import GroupingString
from .revComplementry import get_rev_complementry

class Repeat():
    """docstring for Genotype"""
    def __init__(self,read, last_unit_index,window, number_of_units=0, repeat_units=["CTG"], unique_repeat_units_list=None):
        self.repeat_units = repeat_units
        self.last_unit_index = last_unit_index
        self.number_of_units = number_of_units
        self.start_index = last_unit_index - ((number_of_units+1)*len(repeat_units[0])) #for future use
        self.single_point_mutation_indexes = []
        self.unconfirmed_units_buffer = 0
        self.non_perfect_units = 0
        self.repeat_sequence = ""
        self.unconfirmed_sequence = "" #sequence of non pure repeates waiting for a confirmed unit to  be added
        self.unique_repeat_units_list = unique_repeat_units_list
        self.unique_repeat_units_count = 0
        self.add_unit(window, last_unit_index)
        self.read = read


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
    '''
    def get_SNP_index(self, window, index):
        for i in range(0,len(window)):
            if (window[i] != self.repeat_unit[i]) :
                return index -len(window)+i
	'''
    def change_last_unit_index(self, index):
        self.last_unit_index = index

    def get_non_perfect_units_percentage(self):
    	return self.non_perfect_units/self.number_of_units

    def get_seq(self):
    	return self.read[self.start_index:self.last_unit_index]
    
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
        