#first inport the AGORA model file names as a list into your python workspace, replace path with your own
with open('/Users/theogoodchild-michelman/Desktop/research/MGH/final_healthy_abun/IBD_final_tables/AGORA2_fullnames.txt','r') as input:
    
    AGORA_names = []
    for i in input:
        AGORA_names.append(i.strip())


#first pip install difflib on your terminal (not on python)
#species_list is a list that has all of the original species names from the CDGEM excel
import difflib
AGORA_name = []
spec_list = []
for name in species_list: #species_list is a list that has all of the original species names from the CDGEM excel
    match = difflib.get_close_matches(name, AGORA_names, cutoff=0.6)   
    if match: 
        AGORA_name.append(match[0])
        print(name)
        print(match[0])
        print("---------------------------------")
    else:
        spec_list.append(name)

#the proper names are now saved in the spec_list list
#you can inspect the printed output to see if difflib matched the names well

#SAVING THE correct NAMES AS A CSV FILE, you can rename myfile for whatever name you want

import csv

with open(..., 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(mylist)

     
