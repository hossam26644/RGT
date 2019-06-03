'''docstring'''
from .Repeat import Repeat

class Genotype():
    """docstring for Genotype"""
    def __init__(self, reads, repeat_units=["CTG","CCG"], min_size_repeate=5, max_interrupt_tract=5, unique_repeat_units=None):
        self.repeat_units = repeat_units
        self.reads = reads
        self.min_size_repeate = min_size_repeate
        self.max_interrupt_tract = max_interrupt_tract + len(repeat_units[0])
        self.geno_table = {}
        self.counts_table = {}
        self.get_repeates()
        if unique_repeat_units == None:
            self.unique_repeat_units = self.repeat_units
        else:
            self.unique_repeat_units = unique_repeat_units

    def get_repeates(self):

        for line in self.reads:
            self.get_line_repeates(line)


    def get_line_repeates(self, sequence):
        window_length = len(self.repeat_units[0])

        i = window_length
        repeat = None

        while i < len(sequence): #sliding window
            window = sequence[i-window_length:i]
            if self.window_enters_repeat_sequence(window, self.repeat_units, repeat):
                '''if window detects a repeat unit, while it is not inside a repeat sequence'''
                repeat = Repeat(sequence, i, window,repeat_units=self.repeat_units) #creat a repeat object
                i = i+3 #Jumb one window
                continue

            elif self.detect_repeat_unit_inside_repeat(window, self.repeat_units, repeat):
                '''if it detects a repeat while inside the repeat sequence'''
                repeat.add_unit(window,i) #add a repeat unit count
                i = i+3 #jumb one window
                continue

            elif self.non_matching_unit_within_repeat(window, self.repeat_units, repeat):
                #print(window,i,repeat.last_unit_index, i-repeat.last_unit_index)
                if i-repeat.last_unit_index <= self.max_interrupt_tract:
                    #ignore if length is smaller than max interrupt tract
                    i += 1
                    continue
                #if length is larger than max interrupt tract
                if repeat.number_of_units >= self.min_size_repeate: #check that number of repeates is larger than the minimum size repeate
                    self.add_repeat_to_tables(repeat)

                repeat = None
            i +=1 

        if repeat != None and repeat.number_of_units >= self.min_size_repeate: #if sequence ends on a repeat
            self.add_repeat_to_tables(repeat)




    def non_matching_unit_within_repeat(self, window, repeat_units,  repeat_object):
        window_inside_repeates_flag = repeat_object != None
        if window_inside_repeates_flag:
            for repeat_unit in repeat_units:
                if self.is_window_equals_repeat_unit(window, repeat_unit, repeat_object):
                    return False
            return True 
        return False   

    def detect_repeat_unit_inside_repeat(self, window, repeat_units, repeat_object):
        window_inside_repeates_flag = repeat_object != None
        if window_inside_repeates_flag:
            for repeat_unit in repeat_units:
                if self.is_window_equals_repeat_unit(window, repeat_unit, repeat_object):
                    return True
        return False

    def window_enters_repeat_sequence(self, window, repeat_units, repeat_object): 
        window_inside_repeates_flag = repeat_object != None
        if not window_inside_repeates_flag:
            for repeat_unit in repeat_units:
                if self.is_window_equals_repeat_unit(window, repeat_unit, repeat_object):
                    return True
        return False

    def is_window_equals_repeat_unit(self, window, repeat_unit, repeat_object):
        mismatch = False
        for i in range(0,len(window)):
            if (window[i] != repeat_unit[i]) and (not mismatch) and (repeat_object != None):
                mismatch = True
            elif (window[i] != repeat_unit[i]) and (mismatch or repeat_object == None):
                return False
        return True

    def add_repeat_to_tables(self, repeat):
        self.add_repeat_to_genotable(repeat)
        self.add_repeat_to_countstable(repeat)

    def add_repeat_to_genotable(self, repeat):
        if repeat.get_non_perfect_units_percentage() <= 0.3: #only add repeates with unique percentage > 0.3
            #number_of_repeat_units = repeat.number_of_units
            repeat_sequence = repeat.get_seq()
            
            if(repeat_sequence in self.geno_table):
                self.geno_table[repeat_sequence][0] += 1
            else:
                self.geno_table[repeat_sequence] = [1,repeat.number_of_units]
   
    def add_repeat_to_countstable(self, repeat):
        if repeat.get_non_perfect_units_percentage() <= 0.3: #only add repeates with unique percentage > 0.3
            number_of_repeat_units = repeat.number_of_units
           
            if(number_of_repeat_units in self.counts_table):
                self.counts_table[number_of_repeat_units] += 1
            else:
                self.counts_table[number_of_repeat_units] = 1



    def get_geno_table(self):
        return self.geno_table

    def get_counts_table(self):
        return self.counts_table
