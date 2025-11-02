'''docstring'''

class SmartString():
    """docstring for Genotype"""
    @staticmethod
    def get_buffer_smart_string(string_buffer, window_size, preferabel_repeats):
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
