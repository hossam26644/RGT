'''docstring'''

class ReadFile():
    """docstring for main"""
   
    def __init__(self, file_name, start_flank=None, end_flank=None,
                number_of_allowed_flank_point_mutations=1, discard_reads_with_no_end_flank=True):
        self.get_raw_reads(file_name)
        self.number_of_discarded_reads = 0
        if start_flank == None and end_flank == None:
            self.reads = self.raw_reads
        else:
            self.start_flank = start_flank
            self.end_flank = end_flank
            self.number_of_allowed_flank_point_mutations = number_of_allowed_flank_point_mutations
            self.discard_reads_with_no_end_flank = discard_reads_with_no_end_flank
            self.reads = self.extract_reads_between_flanks()


    def get_raw_reads(self,file_name):
        text_file = open(file_name, "r")
        self.raw_reads = []
        self.all_lines = text_file.read().splitlines()
        for i in range(0,len(self.all_lines)):
            if (i%4==1):
                self.raw_reads.append(self.all_lines[i]) 

    def extract_reads_between_flanks(self):
        reads = []
        start_flank_length = len(self.start_flank)
        end_flank_length = len(self.end_flank)

        for line in self.raw_reads: #loop through all raw reads line by line
            seq_start_index = -1
            seq_end_index = -1
            end_flank_found = False
            for i in range(start_flank_length, len(line)-end_flank_length): #loop through the read searching for the start flank
                current_window = line[i-start_flank_length:i]
                if self.flank_equal(current_window, self.start_flank): #start flank found
                    seq_start_index = i 
                    break
            if seq_start_index == -1: #no start flank found
                self.number_of_discarded_reads += 1
                continue #skip the checking for the end flank

            for i in range(seq_start_index,len(line)-end_flank_length+1): #loop through the remaining of the read searching for the end flank
                current_window = line[i:i+end_flank_length]
                if self.flank_equal(current_window,self.end_flank):
                    seq_end_index = i
                    end_flank_found = True
                    break

            if end_flank_found or not(self.discard_reads_with_no_end_flank):
                reads.append(line[seq_start_index:seq_end_index])
            else:
                self.number_of_discarded_reads += 1


            
        return reads
    
    def flank_equal(self, window, flank):
        current_mismatches = 0
        for i in range(len(window)):
            if window[i] != flank[i]:
                current_mismatches += 1
                if current_mismatches > self.number_of_allowed_flank_point_mutations:
                    return False             
        return True

    def get_discarded_reads_percentage(self):
        percentage = self.number_of_discarded_reads / len(self.raw_reads)
        percentage *= 100
        return percentage
