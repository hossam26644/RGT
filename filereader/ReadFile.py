'''docstring'''
from genotyper.Repeat import Read

class ReadFile():
    """docstring for main"""
   
    def __init__(self, filename, settings):
       
        self.settings = settings
        self.start_flank = settings["start_flank"]
        self.end_flank = settings["end_flank"]
        self.number_of_allowed_strt_flank_point_mutations = settings["number_of_allowed_strt_flank_point_mutations"]
        self.number_of_allowed_end_flank_point_mutations = settings["number_of_allowed_end_flank_point_mutations"]
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
        sf_len, ef_len = len(self.start_flank), len(self.end_flank)
        max_sd, max_ed = self.number_of_allowed_strt_flank_point_mutations, self.number_of_allowed_end_flank_point_mutations
        rr = read.raw_read

        best_sd, seq_start_index, start_flank_start = max_sd + 1, -1, -1
        for i in range(len(rr)):
            for w in sorted(range(sf_len - max_sd, sf_len + max_sd + 1), key=lambda x: abs(x - sf_len)):
                j = i - w
                if j < 0: continue
                if abs(w - sf_len) >= best_sd: continue
                d = Read.levenshtein(rr[j:i], self.start_flank)
                if d < best_sd:
                    best_sd, seq_start_index, start_flank_start = d, i, j
                    read.start_flank_found = True
                    if d == 0: break
            if best_sd == 0: break

        if seq_start_index == -1:
            self.number_of_discarded_reads += 1
            return

        best_ed, seq_end_index, end_flank_end = max_ed + 1, -1, -1
        for i in range(seq_start_index, len(rr)):
            for w in sorted(range(ef_len - max_ed, ef_len + max_ed + 1), key=lambda x: abs(x - ef_len)):
                if i + w > len(rr): continue
                if abs(w - ef_len) >= best_ed: continue
                d = Read.levenshtein(rr[i:i + w], self.end_flank)
                if d < best_ed:
                    best_ed, seq_end_index, end_flank_end = d, i, i + w
                    read.end_flank_found = True
                    if d == 0: break
            if best_ed == 0: break

        read.start_flank_seq = rr[start_flank_start:seq_start_index]
        read.end_flank_seq = rr[seq_end_index:end_flank_end]
        read.interflanking_seq = rr[seq_start_index:seq_end_index]
        read.quality_score = read.quality_score[seq_start_index:seq_end_index]

        if self.settings["report_consensus_flanking_sequence"]:
            start_len = seq_start_index - self.settings["report_consensus_flanking_sequence"]
            start_len = max(0, start_len)
            read.pre_repeat_structure = rr[start_len:seq_start_index]

            end_len = seq_end_index + self.settings["report_consensus_flanking_sequence"]
            end_len = min(len(rr), end_len)
            read.post_repeat_structure = rr[seq_end_index:end_len]

        if read.end_flank_found or not self.discard_reads_with_no_end_flank:
            read.successfully_extracted = True
        else:
            self.number_of_discarded_reads += 1

    def get_discarded_reads_percentage(self):
        percentage = len([self.reads[i] for i in range(len(self.reads)) if not self.reads[i].successfully_extracted])
        percentage /= len(self.reads)
        percentage *= 100
        return percentage        

