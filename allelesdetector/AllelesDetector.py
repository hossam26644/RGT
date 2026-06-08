'''docstring'''
from .PeakIdentifier import PeakIdentifier
from .MatchingSequence import MatchingSequence
import pandas as pd
from collections import Counter


class AllelesDetector():
    """docstring for main"""
    def __init__(self, genotype, settings):

        self.genotype = genotype        
        self.settings = settings
        self.pcr_free = self.settings["PCR_free"]        
        #sorted_geno_list: sorted by abundance, 2D list => list of lists(key index 0, value 1)
        self.sorted_geno_list = (sorted(self.genotype.get_geno_table().items(), key=lambda x: x[1][0], reverse=True))
        self.color_code = ""
        self.possible_alleles_list = self.get_possible_alleles_list_from_sorted_geno_list()  
        peak_identifier = PeakIdentifier(self.genotype.get_counts_table())
        self.peak_repeat_counts = peak_identifier.get_peaks()
        self.first_allele = None
        self.second_allele = None


        self.result_summery = self.predict_alleles()
        if self.first_allele.abundance < self.settings["minimum_no_of_reads"]:
            self.color_code = "red"
            self.result_summery[2] = "Number of reads is lower than threshold"

        self.longer_allele_repeats_count = max(self.first_allele.repeat_units_count,
            self.second_allele.repeat_units_count)

        if self.settings['discard_reads_with_no_end_flank'] == True:
            self.check_no_expanded_alleles_with_no_end_flank()

        self.flanking_sequence_table = None
        self.allele1_most_common_start_flank = None
        self.allele1_most_common_end_flank = None
        self.allele2_most_common_start_flank = None
        self.allele2_most_common_end_flank = None
        
        if settings["report_consensus_flanking_sequence"]:
            self.report_consensus_flanking_sequence()
        else:
            self.report_allele_flanking_sequences()
                   

    def predict_alleles(self):
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
                        self.first_allele = matching_sequences[0]
                        self.second_allele = matching_sequences[1]
                        return([matching_sequences[0].sequence_string,matching_sequences[1].sequence_string, message,
                                str(self.peak_repeat_counts)])
                else:
                    try:
                        neighbour_seq = self.get_neighbour_seq_if_it_is_an_allele(matching_sequences[0])
                        if neighbour_seq != None:
                            self.color_code = "green" 
                            self.first_allele = matching_sequences[0]
                            self.second_allele = neighbour_seq
                            return([matching_sequences[0].sequence_string,neighbour_seq.sequence_string,
                                "Heterozygous", str(self.peak_repeat_counts)])
                    except ValueError as e:
                        self.color_code = "red"
                        self.first_allele = matching_sequences[0]
                        self.second_allele = matching_sequences[1]
                        return([matching_sequences[0].sequence_string,matching_sequences[1].sequence_string, str(e),
                                str(self.peak_repeat_counts)])
                    

                message += " ,two matches with one repeat count, peak may not be a true peak, please check manually"
                self.color_code = "red"

            if self.pcr_free:
                chosen_alleles = set([matching_sequences[0].sequence_string, matching_sequences[1].sequence_string])
                top_alleles = set([self.sorted_geno_list[0][0], self.sorted_geno_list[1][0]])
                if chosen_alleles != top_alleles:
                    message += " ,most abundant repeat sequences not selected, please check manually"
                    self.color_code = "yellow"

            self.first_allele = matching_sequences[0]
            self.second_allele = matching_sequences[1]
            return([matching_sequences[0].sequence_string,matching_sequences[1].sequence_string, message,
                    str(self.peak_repeat_counts)])
        
        #Second case, more than two matches, faulty read
        elif len(matching_sequences) >2:
            self.color_code = "red"
            self.first_allele = matching_sequences[0]
            self.second_allele = matching_sequences[1] 
            return([matching_sequences[0].sequence_string, matching_sequences[1].sequence_string,
                    "More than two potential alleles, please check manually", str(self.peak_repeat_counts)])
        
        #Third case, one peak is identified
        elif len(matching_sequences) == 1:

            #check if the next repeat is an allele
            try:
                neighbour_seq = self.get_neighbour_seq_if_it_is_an_allele(matching_sequences[0])
            except ValueError as e:
                self.color_code = "red"
                self.first_allele = matching_sequences[0]
                self.second_allele = matching_sequences[1]
                return([matching_sequences[0].sequence_string,matching_sequences[1].sequence_string, str(e),
                        str(self.peak_repeat_counts)])
            
            #new_matching_sequences = self.check_if_the_next_repeat_is_an_allele(matching_sequences[0])
            if neighbour_seq != None :
                self.color_code = "green"
                self.first_allele = matching_sequences[0]
                self.second_allele = neighbour_seq                
                return([matching_sequences[0].sequence_string,neighbour_seq.sequence_string,
                       "Heterozygous", str(self.peak_repeat_counts)])
            
            if self.peaks_list_has_peaks_bigger_than_genotyped_alleles(matching_sequences):
                other_allele = self.get_seq_from_matching_peaks_more_than_counts_of(matching_sequences[0])
                self.color_code = "yellow"
                self.first_allele = matching_sequences[0]
                self.second_allele = other_allele              
                return([matching_sequences[0].sequence_string, other_allele.sequence_string,
                       "Heterozygous, expanded allele detected with small count, please check",
                       str(self.peak_repeat_counts)])

            message = "Homozygous"
            self.color_code = "green"
            if matching_sequences[0].repeat_units_count != self.possible_alleles_list[0][1][1]:
                message += " ,most abundant repeat sequence not selected, please chek manually"
                self.color_code = "red"
            self.first_allele = matching_sequences[0]
            self.second_allele = matching_sequences[0]                
            return([matching_sequences[0].sequence_string, matching_sequences[0].sequence_string,
                message, str(self.peak_repeat_counts)] )
        
        
        #Fourth case, no mathces are found
        elif len(matching_sequences) == 0:
            self.color_code = "red"
            new_matching_sequences = self.explore_if_a_close_allele_exists()
            if len(new_matching_sequences)==2 :
                self.first_allele = new_matching_sequences[0]
                self.second_allele = new_matching_sequences[1]                
                return([new_matching_sequences[0].sequence_string ,new_matching_sequences[1].sequence_string,
                    "Heterozygous, please check manually", str(self.peak_repeat_counts)])
            
            elif len(new_matching_sequences)==1:
                self.first_allele = new_matching_sequences[0]
                self.second_allele = new_matching_sequences[0]                
                return([new_matching_sequences[0].sequence_string , new_matching_sequences[0].sequence_string,
                    "Homozygous, please check manually" , str(self.peak_repeat_counts) ] )
           
            return([self.sorted_geno_list[0][1][0],self.sorted_geno_list[1][1][0],
                "No possible alleles found, please check manually", str(self.peak_repeat_counts) ])


    def peaks_list_has_peaks_bigger_than_genotyped_alleles(self, matching_sequences):

        largest_detected_peak = matching_sequences[0].repeat_units_count
       
        for matching_sequence in matching_sequences:
            if matching_sequence.repeat_units_count > largest_detected_peak:
                largest_detected_peak = matching_sequence.repeat_units_count
        for peak in self.peak_repeat_counts:
            if peak > largest_detected_peak:
                return True
        
        return False
  

    def get_neighbour_seq_if_it_is_an_allele(self, matched_sequence):

        if self.pcr_free: 
            plus_sequence_count = matched_sequence.repeat_units_count+1
            minus_sequence_count = matched_sequence.repeat_units_count-1

            candidate_plus = self.get_seq_possible_alleles_list_by_repeats_count(plus_sequence_count,
                                                                    self.possible_alleles_list)
            candidate_minus = self.get_seq_possible_alleles_list_by_repeats_count(minus_sequence_count,
                                                                    self.possible_alleles_list)
            if (candidate_minus is None or candidate_minus.abundance == 0) and (
                candidate_plus is None or candidate_plus.abundance == 0):
                return None
            elif candidate_plus and candidate_plus.abundance >= 1 and (
                    candidate_minus is None or (candidate_plus.abundance > candidate_minus.abundance)):
                return candidate_plus
            elif candidate_minus and candidate_minus.abundance >= 1 and (
                    candidate_plus is None or (candidate_minus.abundance > candidate_plus.abundance)):
                return candidate_minus
            else:
                raise ValueError("Unexpected case, both candidate plus and candidate minus \
                    have abundance higher than 50% of the matched sequence, please check manually")

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
                matching_sequence = MatchingSequence(possibe_allele, idx)
                return matching_sequence
        return None


    def get_seq_from_matching_peaks_more_than_counts_of(self, matching_sequence):
        detected_peak = matching_sequence.repeat_units_count
        repeat_counts_bigger_than_detected_allele = [x for x in self.peak_repeat_counts if x > detected_peak]
        
        for idx, possibe_allele in enumerate(self.sorted_geno_list):
            if (possibe_allele[1][1]) in repeat_counts_bigger_than_detected_allele:
                matching_sequence = MatchingSequence(possibe_allele, idx)
                return matching_sequence
        return MatchingSequence("can't find the other allele", 0, 0, 0, 0) #this line should be impossible to reach


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
                matching_sequence = MatchingSequence(possibe_allele, idx)

                new_matching_sequences.append(matching_sequence)
        return new_matching_sequences
        

    def get_matches_between_peaks_and_possible_alleles_list(self):
        matching_sequences = []
        for idx, possibe_allele in enumerate(self.possible_alleles_list):
            if possibe_allele[1][1] in self.peak_repeat_counts:
                matching_sequence = MatchingSequence(possibe_allele, idx)
                matching_sequences.append(matching_sequence)
        return matching_sequences
    

    def get_possible_alleles_list_from_sorted_geno_list(self):
        most_abundant_allele_values = self.sorted_geno_list[0][1]
        
        # minimum threshold for identifing a repeat sequence as possible allele
        if self.pcr_free:
            min_threshold_abundance = 1 #one read
        else:
            min_threshold_abundance = most_abundant_allele_values[0] * self.settings["min_peak_percentage_threshold"]

        for i in range(0,len(self.sorted_geno_list)):
            #first index the repeat, second the values, third the repeat count   
            if self.sorted_geno_list[i][1][0] < min_threshold_abundance: 
                return self.sorted_geno_list[0:i]
        
        return self.sorted_geno_list


    def check_no_expanded_alleles_with_no_end_flank(self):
        longer_allele_repeats_count = max(self.first_allele.repeat_units_count,
            self.second_allele.repeat_units_count)
        max_detected_expanded = 0
        for read in self.genotype.reads:
            # only check reads with no end flank
            if read.start_flank_found and not read.end_flank_found:
                for repeat in read.repeats:
                    if repeat.number_of_units and repeat.number_of_units > longer_allele_repeats_count:
                        max_detected_expanded = max(max_detected_expanded, repeat.number_of_units)
                        self.color_code = "red"
        
        if max_detected_expanded > longer_allele_repeats_count:
            self.result_summery[2] += f" ,Expanded allele detected with no end flank"
            self.result_summery[2] += f" ,max detected expanded allele has {max_detected_expanded} repeat units, please check manually"

    def report_consensus_flanking_sequence(self):
        """ for repeats making allele1 and allele 2 get every
            allele.pre_repeat_structure and allele.post_repeat_structure
            then align the pre_repeat structure of each allele by the right most
            base (as there will be shorter fragments), and align the post_repeat 
            structure of each allele by the left most base (as there will be shorter fragments)
            Then, at each position, identify the most common base and construct a consensus 
            sequence for the pre_repeat_structure and post_repeat_structure of each allele.
        """
        def consensus(sequences, align_right=False):
            if not sequences:
                return ""
            max_len = max(len(s) for s in sequences)
            padded = [s.rjust(max_len) if align_right else s.ljust(max_len) for s in sequences]
            result = []
            for i in range(max_len):
                col = Counter(s[i] for s in padded if s[i] != ' ')
                if not col:
                    result.append(' ')
                elif len(col) > 1 and col.most_common(2)[0][1] == col.most_common(2)[1][1]:
                    result.append('*')
                else:
                    result.append(col.most_common(1)[0][0])
            return ''.join(result)
        
        def collect(allele):
            pre, post = [], []
            for r in self.genotype.repeats:
                if r.repeat_sequence != allele.raw_seq:
                    continue
                if r.start_flank_found:
                    pre.append(r.pre_repeat_structure)
                if r.end_flank_found:
                    post.append(r.post_repeat_structure)
            return {"pre": Counter(pre), "post": Counter(post)}

        frames = []
        for allele in (self.first_allele, self.second_allele):
            counts = collect(allele)
            seqs_pre = list(counts["pre"].elements())
            seqs_post = list(counts["post"].elements())
            allele.pre_repeat_consensus = consensus(seqs_pre, align_right=True)
            allele.post_repeat_consensus = consensus(seqs_post, align_right=False)

            for position, counter in counts.items():
                df = pd.DataFrame(
                    counter.items(),
                    columns=[f"{allele.sequence_string} {position} flank structures",
                            f"{allele.sequence_string} {position} flank abundance"]
                )
                frames.append(df)

        self.flanking_sequence_table = pd.concat(frames, axis=1)

        self.result_summery.extend([
            self.first_allele.pre_repeat_consensus,
            self.first_allele.post_repeat_consensus,
            self.second_allele.pre_repeat_consensus,
            self.second_allele.post_repeat_consensus,
        ])
                
                        
    def report_allele_flanking_sequences(self):
        for allele in (self.first_allele, self.second_allele):
            counts = self._count_flanks(allele, self.genotype.reads)
            allele.most_common_start_flank = counts["start"].most_common(1)[0][0] if counts["start"] else None
            allele.most_common_end_flank = counts["end"].most_common(1)[0][0] if counts["end"] else None
            allele._flank_counts = counts

        self.flanking_sequence_table = pd.concat([
            self._counts_to_df(a, a._flank_counts)
            for a in (self.first_allele, self.second_allele)
        ], axis=1)

        self.result_summery.extend([
            self.first_allele.most_common_start_flank,
            self.first_allele.most_common_end_flank,
            self.second_allele.most_common_start_flank,
            self.second_allele.most_common_end_flank,
        ])
    
    def _counts_to_df(self, allele: MatchingSequence, counts: dict) -> pd.DataFrame:
        frames = []
        for position in ("start", "end"):
            df = pd.DataFrame(
                counts[position].items(),
                columns=[f"{allele.sequence_string} {position} flank structures",
                        f"{allele.sequence_string} {position} flank abundance"]
            )
            frames.append(df)
        return pd.concat(frames, axis=1)

    def _count_flanks(self, allele: MatchingSequence, reads) -> dict:

        """Returns {position: Counter} for 'start' and 'end' flanks of given allele."""
        counts = {"start": Counter(), "end": Counter()}
        for repeat in self.genotype.repeats:
            if repeat.repeat_sequence != allele.raw_seq:
                continue
            if repeat.start_flank_found:
                counts["start"][repeat.start_flank_seq] += 1
            if repeat.end_flank_found:
                counts["end"][repeat.end_flank_seq] += 1
        return counts
