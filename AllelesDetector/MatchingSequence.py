'''docstring'''

class MatchingSequence():
    """docstring for Genotype"""
    def __init__(self, sequence_string, repeat_units_count, abundance, order, unique_units_count):
        self.sequence_string = sequence_string
        self.repeat_units_count = repeat_units_count
        self.abundance = abundance
        self.order_in_genotable = order+1 #idx starts from 1 not zero
        self.unique_units_count = unique_units_count
   
    def __eq__(self, other): 
        if not isinstance(other, MatchingSequence):
            # don't attempt to compare against unrelated types
            return NotImplemented

        return self.sequence_string == other.sequence_string