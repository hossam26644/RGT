'''
Rewrites strings by detecting repeated adjacent
substrings and converting them into a compact representation.

Example:
    "abcabcabc" -> "[abc]3"
    "abababab"  -> "[ab]4"

The module implements a sliding-window pattern matching algorithm
using a "slow runner / fast runner" approach:

- The slow runner defines a candidate repeat unit.
- The fast runner checks whether the same substring repeats
  immediately after it.
- Consecutive matches are counted and encoded as:
      [repeat_unit]<count>
'''

class SmartString():
    """docstring for Genotype"""
    @staticmethod
    def get_buffer_smart_string(string_buffer, window_size, preferabel_repeats):
        '''
        This method is applied to raw seqeunces.
        If the buffer is large enough, it attempts another compression pass
        using a slightly larger window size.

        Args:
            string_buffer (str):
                The substring buffer to optimize.

            window_size (int):
                Current substring comparison window size.

            preferabel_repeats (list[str]):
                Optional list of preferred repeat patterns that should
                be favored during alignment fine-tuning.

        Returns:
            str:
                The compressed or original buffer string.
        '''        
        if len(string_buffer) > 2*window_size:
            string_buffer = SmartString.get_smart_string_from_sequence(string_buffer, window_size+1, preferabel_repeats)
        return (string_buffer)


    @staticmethod
    def get_smart_string_from_sequence(sequence, window_size, preferabel_repeats=[]):
        ''' fast runner slow runner algorithm 
            both runners are windows with length equals repeat size
            slow runner defines a repeat unit, fast runner (initialy after the 
            slow runner with no gap) checks if both are identical
            if they are, fast runner moves by one window size and adds to the repeat counter


       The algorithm compares two equal-sized sliding windows:
            - slow window
            - fast window

        If both windows contain identical substrings, the substring is
        treated as a repeat unit and consecutive repetitions are counted.

        Example:
            sequence = "abcabcabc"
            window_size = 3

            returns:
                "[abc]3"

        Compression only occurs for directly adjacent repeated patterns.

        Args:
            sequence (str):
                Input string to analyze.

            window_size (int):
                Size of the comparison windows.

            preferabel_repeats (list[str], optional):
                Preferred repeat patterns used during alignment adjustment.
                Defaults to an empty list.

        Returns:
            str:
                Compressed smart-string representation of the sequence.
                    
        '''
        slow_index = 0
        fast_index = window_size
        smart_string = ""
        repeat_unit = ""
        string_buffer = ""
        number_of_repeat_units = 0

        while fast_index < len(sequence) + window_size:
            slow_window = sequence[slow_index:slow_index+window_size]
            fast_window = sequence[fast_index:fast_index+window_size]


            if (slow_window == fast_window) and (repeat_unit == ""):
                smart_string += SmartString.get_buffer_smart_string(string_buffer, window_size, preferabel_repeats)
                string_buffer = ""

                slow_index,fast_index, smart_string = SmartString.fine_tune(sequence, slow_index,
                                                              fast_index, preferabel_repeats,
                                                              window_size, smart_string)
                slow_window = sequence[slow_index:slow_index+window_size]
                fast_window = sequence[fast_index:fast_index+window_size]

                repeat_unit += "["+slow_window + "]"
                fast_index += window_size
                number_of_repeat_units =2

            elif (slow_window == fast_window) and (repeat_unit != ""):
                number_of_repeat_units += 1
                fast_index += window_size

            elif (slow_window != fast_window) and (repeat_unit != ""):
                smart_string += SmartString.get_buffer_smart_string(string_buffer, window_size, preferabel_repeats)
                string_buffer = ""

                smart_string += repeat_unit + str(number_of_repeat_units)
                repeat_unit = ""
                number_of_repeat_units = 0
                slow_index = fast_index
                fast_index += window_size

            elif (slow_window != fast_window) and (repeat_unit == ""):
                string_buffer += slow_window[0]
                #smart_string += slow_window
                fast_index += 1
                slow_index += 1
        
        smart_string += SmartString.get_buffer_smart_string(string_buffer, window_size, preferabel_repeats)
        string_buffer = ""
        return(smart_string)

    def fine_tune(sequence, slow_index, fast_index, preferabel_repeats, window_size, smart_string):

        '''
        Adjust repeat alignment to favor preferred repeat patterns.

        This method shifts the slow/fast windows within the current
        range to check whether a more desirable repeat unit exists.

        If a preferred repeat pattern is found, the indices are updated
        so compression starts at the improved alignment.

        Example:
            Instead of:
                "a[bca]3"

            Fine tuning may produce:
                "[abc]3a"

        Args:
            sequence (str):
                Full input sequence.

            slow_index (int):
                Current slow window start index.

            fast_index (int):
                Current fast window start index.

            preferabel_repeats (list[str]):
                List of preferred repeat units.

            window_size (int):
                Current comparison window size.

            smart_string (str):
                Current compressed output string.

        Returns:
            tuple:
                (
                    updated_slow_index,
                    updated_fast_index,
                    updated_smart_string
                )
        '''        
        for i in range (1, window_size):
            new_slow_index = slow_index +i
            new_slow_window = sequence[new_slow_index:new_slow_index+window_size]
            new_fast_index = fast_index+i
            new_fast_window = sequence[new_fast_index:new_fast_index+window_size]
           
            if (new_slow_window in preferabel_repeats) and new_slow_window == new_fast_window:
                smart_string =  smart_string + sequence [slow_index:new_slow_index] 
                slow_index = new_slow_index
                fast_index = new_fast_index
                break
        return slow_index, fast_index, smart_string
