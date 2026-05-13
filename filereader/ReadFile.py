'''docstring'''
from genotyper.Repeat import Read

class ReadFile():
    """docstring for main"""
   
    def __init__(self, filename, settings):
       
        self.start_flank = settings["start_flank"]
        self.end_flank = settings["end_flank"]
        self.number_of_allowed_frwrd_flank_point_mutations = settings["number_of_allowed_frwrd_flank_point_mutations"]
        self.number_of_allowed_bcwrd_flank_point_mutations = settings["number_of_allowed_bcwrd_flank_point_mutations"]
        self.discard_reads_with_no_end_flank = settings["discard_reads_with_no_end_flank"]
        self.number_of_discarded_reads = 0

        self.reads = self.get_reads(filename)

    def get_reads(self, file_name):
        reads = []

        with open(file_name, "r") as f:
            while True:
                header = f.readline()
                if not header:
                    break

                sequence = f.readline().strip()
                plus = f.readline()
                quality = f.readline()
                read = Read(header, sequence, quality)
                if self.start_flank == None and self.end_flank == None:
                    read.interflanking_seq = sequence
                else:
                    self.extract_reads_between_flanks(read)
                reads.append(read)
        return reads

    def extract_reads_between_flanks(self, read):
        start_flank_length = len(self.start_flank)
        end_flank_length = len(self.end_flank)

        seq_start_index = -1
        seq_end_index = -1
        
        for i in range(start_flank_length, len(read.raw_read)-end_flank_length): #loop through the read searching for the start flank
            current_window = read.raw_read[i-start_flank_length:i]
            if self.flank_equal(current_window, self.start_flank, True): #start flank found
                seq_start_index = i 
                read.start_flank_found = True
                break
        if seq_start_index == -1: #no start flank found
            self.number_of_discarded_reads += 1
            return

        for i in range(seq_start_index,len(read.raw_read)-end_flank_length+1): #loop through the remaining of the read searching for the end flank
            current_window = read.raw_read[i:i+end_flank_length]
            if self.flank_equal(current_window,self.end_flank, False):
                seq_end_index = i
                read.end_flank_found = True
                break
        read.interflanking_seq = read.raw_read[seq_start_index:seq_end_index]
        read.quality_score = read.quality_score[seq_start_index:seq_end_index]

        if read.end_flank_found or not(self.discard_reads_with_no_end_flank):
            read.successfully_extracted = True
        else:
            self.number_of_discarded_reads += 1

    
    def flank_equal(self, window, flank, frwrd=True):
        allowed_mismatches = self.number_of_allowed_frwrd_flank_point_mutations if frwrd \
            else self.number_of_allowed_bcwrd_flank_point_mutations
        
        current_mismatches = 0
        for i in range(len(window)):
            if window[i] != flank[i]:
                current_mismatches += 1
                if current_mismatches > allowed_mismatches:
                    return False             
        return True

    def get_discarded_reads_percentage(self):
        percentage = len([self.reads[i] for i in range(len(self.reads)) if not self.reads[i].successfully_extracted])
        percentage /= len(self.reads)
        percentage *= 100
        return percentage        

