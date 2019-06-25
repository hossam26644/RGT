'''docstring'''
from .PeakIdentifier import PeakIdentifier
from .MatchingSequence import MatchingSequence

class AllelesDetector():
    """docstring for main"""
    def __init__(self, counts_table, geno_table):

        #sorted_geno_list: sorted by abundance, 2D list => list of lists(key index 0, value 1)
        self.sorted_geno_list = (sorted(geno_table.items(), key=lambda x: x[1][0], reverse=True))
        self.color_code = ""
        self.counts_table = counts_table
        self.possible_alleles_list = self.get_possible_alleles_list_from_sorted_geno_list()  
        peak_identifier = PeakIdentifier(counts_table)
        self.peak_repeat_counts = peak_identifier.get_peaks()
        self.plateau = 15 #plateau to ignore small peaks in the region of the detected allele,
                            #in the case of the expanded alleles


    def get_alleles(self):
        matching_sequences = self.get_matches_between_peaks_and_possible_alleles_list()


        #First case, when 2 matches identified
        if len(matching_sequences) == 2:
            message = "Heterozygous"
            self.color_code = "green" 
            #check if an expanded allele exists
            if self.peaks_list_has_peaks_bigger_than_genotyped_alleles(matching_sequences):
                message += " ,Peaks are found after expanded allele, please check manually "
                self.color_code = "yellow"
           
            #check if the two matches have one count, this may mean that the peak is not a true peak
            #the peak could be the overlap of two peaks
            if matching_sequences[0].repeat_units_count == matching_sequences[1].repeat_units_count:
                if (matching_sequences[0].order_in_genotable == 1) and (matching_sequences[1].order_in_genotable == 2):
                    if matching_sequences[1].abundance >= 0.5* matching_sequences[0].abundance:
                        return([matching_sequences[0].sequence_string,matching_sequences[1].sequence_string, message,
                                str(self.peak_repeat_counts)])
                else:
                    neighbour_seq = self.get_neighbour_seq_if_it_is_an_allele(matching_sequences[0])
                    if neighbour_seq != None:
                        self.color_code = "green" 
                        return([matching_sequences[0].sequence_string,neighbour_seq.sequence_string,
                            "Heterozygous", str(self.peak_repeat_counts)])
                    

                message += " ,two matches with one repeat count, peak may not be a true peak, please check manually"
                self.color_code = "red"



            return([matching_sequences[0].sequence_string,matching_sequences[1].sequence_string, message,
                    str(self.peak_repeat_counts)])
        
        #Second case, more than two matches, faulty read
        elif len(matching_sequences) >2:
            self.color_code = "red"
            return([matching_sequences[0].sequence_string, matching_sequences[1].sequence_string,
                    "More than two potential alleles, please check manually", str(self.peak_repeat_counts)])
        
        #Third case, one peak is identified
        elif len(matching_sequences) == 1:

            #check if the next repeat is an allele
            neighbour_seq = self.get_neighbour_seq_if_it_is_an_allele(matching_sequences[0])
            #new_matching_sequences = self.check_if_the_next_repeat_is_an_allele(matching_sequences[0])
            if neighbour_seq != None :
                self.color_code = "green"
                return([matching_sequences[0].sequence_string,neighbour_seq.sequence_string,
                       "Heterozygous", str(self.peak_repeat_counts)])
            
            if self.peaks_list_has_peaks_bigger_than_genotyped_alleles(matching_sequences):
                other_allele = self.get_seq_from_matching_peaks_more_than_counts_of(matching_sequences[0])
                self.color_code = "yellow"
                return([matching_sequences[0].sequence_string, other_allele.sequence_string,
                       "Heterozygous, expanded allele detected with small count, please check",
                       str(self.peak_repeat_counts)])

            message = "Homozygous"
            self.color_code = "green"
            if matching_sequences[0].repeat_units_count != self.possible_alleles_list[0][1][1]:
                message += " ,most abundant repeat sequence not selected, please chek manually"
                self.color_code = "red"
            return([matching_sequences[0].sequence_string, matching_sequences[0].sequence_string,
                message, str(self.peak_repeat_counts)] )
        
        
        #Fourth case, no mathces are found
        elif len(matching_sequences) == 0:
            self.color_code = "red"
            new_matching_sequences = self.explore_if_a_close_allele_exists()
            if len(new_matching_sequences)==2 :
                return([new_matching_sequences[0].sequence_string ,new_matching_sequences[1].sequence_string,
                    "Heterozygous, please check manually", str(self.peak_repeat_counts)])
            
            elif len(new_matching_sequences)==1:
                return([new_matching_sequences[0].sequence_string , new_matching_sequences[0].sequence_string,
                    "Homozygous, please check manually" , str(self.peak_repeat_counts) ] )
           
            return([self.sorted_geno_list[0][1][0],self.sorted_geno_list[1][1][0],
                "No possible alleles found, please check manually", str(self.peak_repeat_counts) ])

    def peaks_list_has_peaks_bigger_than_genotyped_alleles(self, matching_sequences):

        largest_detected_peak = matching_sequences[0].repeat_units_count+self.plateau
       
        for matching_sequence in matching_sequences:
            if matching_sequence.repeat_units_count > largest_detected_peak:
                largest_detected_peak = matching_sequence.repeat_units_count+self.plateau
        for peak in self.peak_repeat_counts:
            if peak > largest_detected_peak:
                return True
        
        return False
  

    def get_neighbour_seq_if_it_is_an_allele(self, matched_sequence):
        candidate_sequence_count = matched_sequence.repeat_units_count+1
        candidate_sequence = self.get_seq_possible_alleles_list_by_repeats_count(candidate_sequence_count,
                                                                    self.possible_alleles_list)
        if (candidate_sequence == None) or (candidate_sequence.abundance < (matched_sequence.abundance*0.5)):
            return None

        sequence_smaller_than_matched_seq_count = matched_sequence.repeat_units_count-1
        sequence_smaller_than_matched_seq = self.get_seq_possible_alleles_list_by_repeats_count(
                                                                sequence_smaller_than_matched_seq_count,
                                                                    self.sorted_geno_list)

        if candidate_sequence.abundance > sequence_smaller_than_matched_seq.abundance*1.1 :
            return candidate_sequence

        return None


    def get_seq_possible_alleles_list_by_repeats_count(self, count, given_list):
        for idx, possibe_allele in enumerate(given_list):
            if (possibe_allele[1][1]) == count:
                matching_sequence = MatchingSequence(possibe_allele[0], possibe_allele[1][1], possibe_allele[1][0], idx)
                return matching_sequence
        return None


    def get_seq_from_matching_peaks_more_than_counts_of(self, matching_sequence):
        detected_peak = matching_sequence.repeat_units_count
        repeat_counts_bigger_than_detected_allele = [x for x in self.peak_repeat_counts if x >= detected_peak+self.plateau]
        
        for idx, possibe_allele in enumerate(self.sorted_geno_list):
            if (possibe_allele[1][1]) in repeat_counts_bigger_than_detected_allele:
                matching_sequence = MatchingSequence(possibe_allele[0], possibe_allele[1][1], possibe_allele[1][0], idx)
                return matching_sequence
        return MatchingSequence("can't find the other allele", 0, 0, 0) #this line should be impossible to reach


    def explore_if_a_close_allele_exists(self):
        new_matching_sequences = []
        new_peak_repeat_counts = list(self.peak_repeat_counts)

        first_key = list(self.counts_table)[0]
        counts_min_threshold = self.counts_table[first_key]* 0.5

        #adding near by points to peak list
        for number_of_repeat_counts in self.counts_table.keys():
            count = self.counts_table[number_of_repeat_counts]
           
            if number_of_repeat_counts in new_peak_repeat_counts:
                continue
            if count > counts_min_threshold:
                new_peak_repeat_counts.append(number_of_repeat_counts)
      
        #checking if a match happens with the near by points, now identified as new peaks
        for idx, possibe_allele in enumerate(self.possible_alleles_list):
            if (possibe_allele[1][1]) in new_peak_repeat_counts:
                matching_sequence = MatchingSequence(possibe_allele[0], possibe_allele[1][1], possibe_allele[1][0], idx)
                new_matching_sequences.append(matching_sequence)
        return new_matching_sequences
        

    def get_matches_between_peaks_and_possible_alleles_list(self):
        matching_sequences = []
        for idx, possibe_allele in enumerate(self.possible_alleles_list):
            if possibe_allele[1][1] in self.peak_repeat_counts:
                matching_sequence = MatchingSequence(possibe_allele[0], possibe_allele[1][1], possibe_allele[1][0], idx)
                matching_sequences.append(matching_sequence)
        return matching_sequences
    
    def get_possible_alleles_list_from_sorted_geno_list(self):
        most_abundant_allele_values = self.sorted_geno_list[0][1]
        
        # minimum threshold for identifing a repeat sequence as possible allele
        min_threshold_abundance = most_abundant_allele_values[0] * 0.2

        for i in range(0,len(self.sorted_geno_list)):
            #first index the repeat, second the values, third the repeat count   
            if self.sorted_geno_list[i][1][0] < min_threshold_abundance: 
                return self.sorted_geno_list[0:i]
        
        return self.sorted_geno_list
