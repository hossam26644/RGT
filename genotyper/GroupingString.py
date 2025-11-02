'''docstring'''

class GroupingString():
    """docstring for Genotype"""

    @staticmethod
    def get_grouped_string_from_sequence(sequence, grouping_units):
        grouped_seq = sequence
        for unit in grouping_units:
            window_ln = len(unit)
            window_strt_idx = 0
            buffer_str = ""
            count = 0

            while window_strt_idx < len(grouped_seq):
                window = grouped_seq[window_strt_idx:window_strt_idx+window_ln]
                
                if window == unit:
                    if count == 0:
                        buffer_str += "[" + unit + "]"
                        count += 1
                        window_strt_idx += window_ln
                    else:
                        count += 1
                        window_strt_idx += window_ln
                
                else:
                    

                    if count != 0:
                        buffer_str += str(count) + window
                        count = 0
                        window_strt_idx += window_ln

                    else:
                        buffer_str += window[0]
                        window_strt_idx += 1

                    if "[" in window: #to skip already grouped repeats
                        new_idx = GroupingString.get_idx_after_the_closing_sq_bracket(window_strt_idx, grouped_seq)
                        buffer_str += grouped_seq[window_strt_idx:new_idx]
                        window_strt_idx = new_idx
                        
            if count != 0:
                buffer_str += str(count)
            grouped_seq = buffer_str
            buffer_str = ""
        return grouped_seq

   
    @staticmethod
    def get_idx_after_the_closing_sq_bracket(idx, sequence):
        for i in range(idx, len(sequence)):
            if sequence[i] == "]":
                return i 