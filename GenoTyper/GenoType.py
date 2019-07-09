'''docstring'''
from .Repeat import Repeat
from .revComplementry import get_rev_complementry

class Genotype():
    """docstring for Genotype"""
    def __init__(self, reads, repeat_units=["CAG","CAA","CCG", "CCA","CCT"], min_size_repeate=5,
                max_interrupt_tract=5, unique_repeat_units=None, grouping_repeat_units=None, reverse_strand=False):
        self.reverse_strand = reverse_strand
        if self.reverse_strand:
            self.repeat_units = get_rev_complementry(repeat_units)
            self.unique_repeat_units = get_rev_complementry(unique_repeat_units)
        else:
            self.repeat_units = repeat_units
        self.list_of_repeat_units_lengths = self.get_list_of_repeat_units_lengths()#used to have diff sliding windows lengths
        self.reads = reads
        self.min_size_repeate = min_size_repeate
        self.max_interrupt_tract = max_interrupt_tract
        self.grouping_repeat_units = grouping_repeat_units
        self.geno_table = {}
        self.counts_table = {}
        self.unique_counts_table = {}

        self.get_repeates()

    def get_repeates(self):

        for line in self.reads:
            self.get_line_repeates(line)


    def get_line_repeates(self, sequence):

        i = 0
        repeat = None

        while i <= len(sequence)-min(self.list_of_repeat_units_lengths): #sliding window
            checker = self.window_enters_repeat_sequence(i, self.repeat_units, repeat, sequence)
            if checker[0]:
                '''if window detects a repeat unit, while it is not inside a repeat sequence'''
                window = checker[1]
                repeat = Repeat(sequence, i, window,repeat_units=self.repeat_units,
                                unique_repeat_units_list=self.unique_repeat_units) #creat a repeat object
                i = i+len(window) #Jumb one window
                continue
            
            checker = self.detect_repeat_unit_inside_repeat(i, self.repeat_units, repeat, sequence)
            if checker[0]:
                '''if it detects a repeat while inside the repeat sequence'''
                window = checker[1]
                repeat.add_unit(window,i+len(window)) #add a repeat unit count
                i = i+len(window) #jumb one window
                continue

            elif self.non_matching_unit_within_repeat(i, self.repeat_units, repeat, sequence):
                #print(window,i,repeat.last_unit_index, i-repeat.last_unit_index)
                if i-repeat.last_unit_index < self.max_interrupt_tract:
                    #ignore if length is smaller than max interrupt tract
                    i += 1
                    continue
                #if length is larger than max interrupt tract
                elif repeat.number_of_units >= self.min_size_repeate: #check that number of repeates is larger than the minimum size repeate
                    self.add_repeat_to_tables(repeat)
                repeat = None
            i +=1 

        if repeat != None and repeat.number_of_units >= self.min_size_repeate: #if sequence ends on a repeat
            self.add_repeat_to_tables(repeat)



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
                elif repeat_object != None and self.hamming_distance(window, repeat_unit)==1:
                    if similar_seq == "":
                        similar_seq = window
        if similar_seq != "":
            return True, similar_seq
        return False, ""

    def hamming_distance(self,s1, s2):
        if len(s1) != len(s2):
            return -1
        return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

    def add_repeat_to_tables(self, repeat):
        self.add_repeat_to_genotable(repeat)
        self.add_repeat_to_countstable(repeat)
        self.add_repeat_to_unique_countstable(repeat)

    def add_repeat_to_genotable(self, repeat):
        if repeat.get_non_perfect_units_percentage() <= 0.3: #only add repeates with unique percentage > 0.3
            if self.grouping_repeat_units == None:
                repeat_sequence = repeat.get_seq_smart_string(self.reverse_strand)
            else:
                repeat_sequence = repeat.get_grouped_string(self.grouping_repeat_units,self.reverse_strand)
            
            if(repeat_sequence in self.geno_table):
                self.geno_table[repeat_sequence][0] += 1
            else:
                self.geno_table[repeat_sequence] = [1,repeat.number_of_units,repeat.unique_repeat_units_count,repeat.get_seq()]
   
    def add_repeat_to_countstable(self, repeat):
        if repeat.get_non_perfect_units_percentage() <= 0.3: #only add repeates with unique percentage > 0.3
            number_of_repeat_units = repeat.number_of_units
           
            if(number_of_repeat_units in self.counts_table):
                self.counts_table[number_of_repeat_units] += 1
            else:
                self.counts_table[number_of_repeat_units] = 1

    def add_repeat_to_unique_countstable(self, repeat):
        if repeat.get_non_perfect_units_percentage() <= 0.3: #only add repeates with unique percentage > 0.3
            number_of_unique_repeat_units = repeat.unique_repeat_units_count
           
            if(number_of_unique_repeat_units in self.unique_counts_table):
                self.unique_counts_table[number_of_unique_repeat_units] += 1
            else:
                self.unique_counts_table[number_of_unique_repeat_units] = 1


    def get_geno_table(self):
        return self.geno_table

    def get_counts_table(self):
        return self.counts_table

    def get_unique_counts_table(self):
        return self.unique_counts_table

    def get_list_of_repeat_units_lengths(self):
        list_of_repeat_units_lengths =[]
        for unit in self.repeat_units:
            list_of_repeat_units_lengths.append(len(unit))
        list_of_repeat_units_lengths = set(list_of_repeat_units_lengths) #set removes repeats
        list_of_repeat_units_lengths = list(list_of_repeat_units_lengths) #making lenghts a list
        list_of_repeat_units_lengths.sort(reverse = True) #sort descendingly
       
        return(list_of_repeat_units_lengths)
