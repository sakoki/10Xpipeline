import re 

def format_to_tuple(colname):
    parts = colname.split(',')
    parts = [re.sub('[\(\'\)]', '', x.strip()) for x in parts]
    return tuple(parts)