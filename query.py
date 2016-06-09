queryDrugs = open("./query_drug_drugbank.txt")

drugs = open("./drugs-marked.txt")

d = []

for l in drugs:
    lineList = l.strip().split("\t")
    d.append([lineList[0].strip(), int(lineList[1])])
    
    
    
drugs.close()

for j in queryDrugs:
    for l in d:
        if l[0] == j.strip():
            l[1] = 1
            
queryDrugs.close()
            
output = open("./query_drug_processed.txt", "w")

for i in d:
    output.write("%s\t%d\n" % (i[0], i[1]))
    
    
output.close()