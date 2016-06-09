from __future__ import division

drugs = []
proteins = []
pp = []
dd = []
newdd = []

    
def parseDrugTarget():
    
    drugFile = open("./drugs.txt") # File containing DB ID and name of drugs
    
    for line in drugFile:
        lineList = line.split("\t")
        drugs.append((lineList[0].strip(), lineList[1].strip(), []))
        
    drugFile.close()
    
    proFile = open("./proteins.txt") # File containing UniProt ID and names of proteins
    
    for line in proFile:
        proteins.append((line.strip(), []))
    drugNameIDFile = open("./drug_target_uniprot_links.txt") 
    
    drugNameIDFile.readline()
    
    for line in drugNameIDFile:
        lineList = line.split("\t")
        newName = lineList[1].strip()
        if newName[0] == '\"':
            newName = newName[1:-1]
        newID = lineList[0].strip()
        newProtein = lineList[3].strip()
        newProteinIndex = -1
        newDrug = 0
        for j in range(len(proteins)):
            if proteins[j][0] == newProtein:
                newProteinIndex = j
                break
            
        for i in drugs:
            if i[0] == newID:
                i[2].append(j) # i[2] has indices of proteins
                newDrug = drugs.index(i)
    
        for j in proteins:
            if (newProtein == j[0]):
                j[1].append(newDrug) #j[1] has indices of drugs
            
        
    drugNameIDFile.close()
        
    print len(drugs), len(proteins)
        
        
        
def initPP():
    
    global dd
    global pp
    
    drugLength = len(drugs)
    proteinLength = len(proteins)
    dd = map(lambda x: map(lambda x: 0.0, range(drugLength)), range(drugLength))
    pp = map(lambda x: map(lambda x: 0.0, range(proteinLength)), range(proteinLength))
    
    #print "0 filled"
    
    '''for i in range(proteinLength):
        pp.append([0.0] * proteinLength)
        
    for i in range(drugLength):
        pp.append([0.0] * drugLength)   ''' 
        
    for i in range(proteinLength):
        pp[i][i] = 1.0
        for j in range(i):
            if set(proteins[i][1]).intersection(proteins[j][1]) != set():
                pp[i][j] = 0.5
                pp[j][i] = 0.5
            #print "%d, %d processed" % (i, j)

def constructDD():
    global dd
    global pp
    drugLength = len(drugs)
    proteinLength = len(proteins)

    for d1 in range(drugLength):
        #print d1
        dd[d1][d1] = 1.0
        for d2 in range(d1):
            temp = 0.0
            d1targets = drugs[d1][2]
            d2targets = drugs[d2][2]
            d1length = len(d1targets)
            d2length = len(d2targets)                
            for p1 in d1targets:
                for p2 in d2targets:
                    temp += pp[p1][p2]
                    
            dd[d1][d2] = temp / (d1length * d2length)
            dd[d2][d1] = dd[d1][d2]
            
    dindex = drugs.index            
    for p1 in range(proteinLength):
        #print p1
        pp[p1][p1] = 1.0
        for p2 in range(p1):
            temp = 0.0
            p1drugs = proteins[p1][1]
            p2drugs = proteins[p2][1]
            p1length = len(p1drugs)
            p2length = len(p2drugs)
            
            for d1 in p1drugs:
                for d2 in p2drugs:
                    temp += dd[d1][d2]
                    
            pp[p1][p2] = temp / (p1length * p2length)
            pp[p2][p1] = pp[p1][p2]                


def writeDD():
    
    outputFile1 = open("./result2.txt", "w")
    outputFile1.write("Drug1\t Drug2\t Correlation\n")
    outputFile2 = open("./result3.txt", "w")
    outputFile2.write("Drug1\t Drug2\t Correlation\n")  
    outputFile3 = open("./result4.txt", "w")
    outputFile3.write("Drug1\t Drug2\t Correlation\n")    
    outputFile4 = open("./result35.txt", "w")
    outputFile4.write("Drug1\t Drug2\t Correlation\n")    
    drugLength = len(drugs)
    for i in range(drugLength):
        for j in range(i):
            if (newdd[i][j] > 0.2):
                outputFile1.write("%s\t%s\t%f\n" % (drugs[i][1], drugs[j][1], dd[i][j]))
            if (newdd[i][j] > 0.3):
                outputFile2.write("%s\t%s\t%f\n" % (drugs[i][1], drugs[j][1], dd[i][j]))                
            if (newdd[i][j] > 0.4):
                outputFile3.write("%s\t%s\t%f\n" % (drugs[i][1], drugs[j][1], dd[i][j]))
            if (newdd[i][j] > 0.35):
                outputFile4.write("%s\t%s\t%f\n" % (drugs[i][1], drugs[j][1], dd[i][j]))                
            
    outputFile1.close()
    outputFile2.close()
    outputFile3.close()
    outputFile4.close()
    
    
def parseDEG():
    
    global dd, newdd
    DEGFile = open("./Table_S2.csv") # Name of Drug-DEG correlation database
    
    DEGFile.readline()
    drugLength = len(drugs)
    corr = map(lambda x: map(lambda x: 0.0, range(drugLength)), range(drugLength))
    
    for line in DEGFile:
        lineList = line.split(",")
        if lineList[8].strip() != "agent":
            continue
        
        agentA = lineList[10].strip()
        agentB = lineList[13].strip()
        
        drugA = -1
        drugB = -1
        
        for k in range(len(drugs)):
            if drugs[k][1] == agentA:
                drugA = k
            elif drugs[k][2] == agentB:
                drugB = k
                
            if (drugA != -1) and (drugB != -1):
                score = float(lineList[2].strip())
                if score > 0:
                    corr[drugA][drugB] = float(lineList[2].strip())
                    corr[drugB][drugA] = corr[drugA][drugB]
                break
            
    newdd = map(lambda x: map(lambda x: 0.0, range(drugLength)), range(drugLength))
    
    for i in range(drugLength):
        for j in range(drugLength):
            newdd[i][j] = 0.5 * corr[i][j] + dd[i][j]
                
            
     

parseDrugTarget()
print "parsing finished"
initPP()
print "init finished"
print len(drugs), len(dd)
for i in range(20):
    constructDD()
    print "%d cycle finished" % i
print "constructing DD finished"

parseDEG()
writeDD()