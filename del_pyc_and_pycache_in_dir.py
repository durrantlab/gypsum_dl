import __future__

import glob
import sys
import copy
import os
def get_folders():

    script_dir = str(os.path.dirname(os.path.realpath(__file__)))
                

    folders_all_list = [x for x in glob.glob(script_dir+"/*/") if os.path.isdir(x)==True]
    bottom_dir = False
    new_folders = copy.deepcopy(folders_all_list)

    while bottom_dir ==False:
        temp = []
        for folder in new_folders:
            sub_list = [x for x in glob.glob(folder+"/*/") if os.path.isdir(x)==True]
            temp.extend(sub_list)
        
        folders_all_list.extend(temp)
        folders_all_list=list(set(folders_all_list))
        if len(temp) == 0:
            bottom_dir = True
            break
        
        else:
            new_folders = copy.deepcopy(temp)


    folders_all_list=list(set(folders_all_list))
    return folders_all_list


folders_all_list = get_folders()

final_dir_list = []
to_del_folder = []
for folder in folders_all_list:
    if "__pycache__" in folder:
        to_del_folder.append(folder)
        os.system("rm -rf {}".format(folder))

    else:
        final_dir_list.append(folder)




files_to_del = []
for folder in folders_all_list:
    temp_file_list = [x for x in glob.glob(folder+"/*.pyc")]
    for file_del in temp_file_list:
            
        if ".pyc" in file_del:
            to_del_folder.append(file_del)
            os.system("rm {}".format(file_del))


print("")
print("folders to start : ", len(folders_all_list))
print("folder del : ", len(files_to_del))
print("files del : ", len(to_del_folder))

print("")
print("folders at end : ", len(get_folders()))
print("")


files_to_del = []
for folder in folders_all_list:
    temp_file_list = [x for x in glob.glob(folder+"/#*#")]
    for file_del in temp_file_list:
        if "#" in file_del:
            print(file_del)

        files_to_del.extend(temp_file_list)

print("len # files to delete: ", len(files_to_del))
for file_del in files_to_del:
    os.system("rm {}".format(file_del))
