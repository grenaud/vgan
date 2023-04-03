from ete3 import NCBITaxa, Tree
import csv
import argparse
ncbi = NCBITaxa()

parser = argparse.ArgumentParser()
parser.add_argument('outFile', help= "euka abundance output file")
args = parser.parse_args()

clades = []
res_list = []
index = {}
tsv_file = open(args.outFile)
read_tsv = csv.reader(tsv_file, delimiter="\t")
for row in read_tsv:
    clades.append(row[0])
    if (row[1] == "yes"):

        res_list.append(row[0])
        index[f"{row[0]}"] = row[2:6]


del clades[0]
taxid_list = []
for i in clades:
    try:
        taxid = ncbi.get_name_translator([f'{i}'])
        res = taxid[f"{i}"]
        taxid_list.append(res)
    except:
        pass

flat_list = [item for sublist in taxid_list for item in sublist] 

not_detected = [item for item in clades if item not in res_list]   

tree = ncbi.get_topology(flat_list)


for node in tree.traverse():

    if (node.common_name == ""):
        node.name = node.sci_name
    else:
        node.name = node.sci_name + ' - ' + node.common_name

    if (node.sci_name in not_detected):
        print(node.sci_name)
        node.name = node.name + ' - ' + "not detected"

    elif (node.sci_name in res_list):
        print(node.sci_name)
        node.name = node.name  + ' - ' + "detected" + ' - ' + index[f"{node.sci_name}"][0] + ' - ' + index[f"{node.sci_name}"][1] + ' - ' + index[f"{node.sci_name}"][2] + ' - ' + index[f"{node.sci_name}"][3]
    

final = tree.get_ascii(attributes=["name"])

print(final)

tree.write(format = 1,outfile="tree.nw")
