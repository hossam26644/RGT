'''docstring'''

class MatchingSequence():
    """docstring for Genotype"""
    def __init__(self, geno_table_element, idx):
   
        self.sequence_string = geno_table_element[0]
        self.repeat_units_count = geno_table_element[1][1]
        self.abundance = geno_table_element[1][0]
        self.order_in_genotable = idx+1 #idx starts from 1 not zero
        self.unique_units_count = geno_table_element[1][4]
        self.x_count = geno_table_element[1][2] 
        self.z_count = geno_table_element[1][3]
        self.raw_seq = geno_table_element[1][5]
        self.most_common_start_flank = None
        self.most_common_end_flank = None
        self.pre_repeat_consensus = None
        self.post_repeat_consensus = None
   
    def __eq__(self, other): 
        if not isinstance(other, MatchingSequence):
            # don't attempt to compare against unrelated types
            return NotImplemented

        return self.sequence_string == other.sequence_string