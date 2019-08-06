
def inverse(character):
    if character == 'A':
        output = 'T'
    elif character == 'C':
        output = 'G'
    elif character == 'G':
        output = 'C'
    elif character == 'T':
        output = 'A'	
    else:
        output = character
    return output


def calc_complementry(text):
    length = len(text)
    array = []
    for i in range (0,length):
        array.append(inverse(text[i]))
    return ''.join(reversed(list(array)))

def get_rev_complementry(_input):
    if(type(_input)==list):
        output = []
        for item in _input:
            output.append(calc_complementry(item))
        return output
    else:
        return calc_complementry(_input)