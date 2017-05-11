import json
import sys
import copy
import os

# Create new dictionary from json 

def new_dict(filename):
    mydict = {}
    data = json.load(open(filename, "r"))
    data.update(mydict)
 # json.dump(data, open('recent.json', "w"))
    return data

# Alter true/falses so only one is present for each output file
def get_bool_keys(raw_dict):
    """
    takes a dictionary and returns a list of keys whose values are boolean
    """
    key_list = []
    for key in raw_dict:
        value = raw_dict[key]
        if isinstance(value, bool):
            key_list.append(key)

    return key_list


def alter_me(key_lists, raw_dict):
    """
    alters values of keys in dictionary
    """    
    dict_list = []
    for key_list in key_lists:
        new_dict = copy.deepcopy(raw_dict)
        for key in key_list:
            new_dict[key] = False
        dict_list.append(new_dict)

    return dict_list

def combiner(key_list):
    """
    recurses through all combinations of present/absent parameters within key_list
    """
    if len(key_list) == 0:
        return [[]] 
    
    working_list = copy.deepcopy(key_list)
    cur_key = working_list.pop(-1)
    start_list = combiner(working_list)
    return_list = []
    for sub_list in start_list:
        copy_list = copy.deepcopy(sub_list)
        copy_list.append(cur_key)
        return_list.append(sub_list)
        return_list.append(copy_list)

    return return_list    

# Save output files 
def output_to_json(dict_list, filename):
    filebase = os.path.basename(filename)
    basename = filebase.strip(".json")
    for i, out_dict in enumerate(dict_list):
        outname = "".join([basename,"_", str(i), ".json"])
        print outname
        with open(outname, "w") as outfile:
            json.dump(out_dict, outfile)

filename = sys.argv[1] 
my_dict = new_dict(filename)
key_list = get_bool_keys(my_dict)
list_o_lists = combiner(key_list)
dict_list = alter_me(list_o_lists, my_dict)
#for item in dict_list:
#    print item
output_to_json(dict_list, filename)
