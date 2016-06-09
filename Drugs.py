drugs = []
proteins = []
pp = []
dd = []

class Drug(object):
    DBid = ""
    name = ""
    DEGs = []
    targets = []
    
    def __init__(self, DBid, name):
        self.DBid = DBid
        self.name = name
        
    def __str__(self):
        return "Drug name: %s, DB ID: %s" % (self.name, self.DBid)
    
    def __eq__(self, rhs):
        return (self.name == rhs.name) and (self.DBid == rhs.DBid)
    
    def getID(self):
        return self.DBid
    
    def getName(self):
        return self.name
    
    def addTarget(self, target):
        self.targets.append(target)      
        
    def addDEG(self, DEG):
        self.DEGs.append(DEG)
        
    def getTargets(self):
        return self.targets
    
    def getDEG(self):
        return self.DEGs
    
    def getTargetCorrelation(self, drug2):
        count = 0
        for i in self.targets:
            if i in drug2.targets:
                count += 1
                
        return count
    
class Protein(object):
    
    ID = ""
    drugs = []
    
    def __init__(self, ID):
        self.ID = ID
        
    def __eq__(self, rhs):
        return self.ID == rhs.ID
    
    def __str__(self):
        return "protein " + self.ID
    
def parseDrugTarget():
    
    drugNameIDFile = open("./test.csv") #open("./drug_target_uniprot_links.csv")
    
    drugNameIDFile.readline()
    
    for line in drugNameIDFile:
        add = True
        addP = True
        lineList = line.split(",")
        newName = lineList[0].strip()
        newID = lineList[1].strip()
        newDrug = Drug(lineList[0].strip(), lineList[1].strip())
        newProtein = Protein(lineList[3].strip())
        for i in drugs:
            if (newDrug == i):
                add = False
                newDrug = i
                
        if add:
            drugs.append(newDrug)
        
           
        for j in proteins:
            if (newProtein == j):
                addP = False
                newProtein = j
        if addP:
            proteins.append(newProtein)
            
        newDrug.targets.append(newProtein)
        newProtein.drugs.append(newDrug)
        
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
            if set(proteins[i].drugs).intersection(proteins[j].drugs) != set():
                pp[i][j] = 0.5
                pp[j][i] = 0.5
            print "%d, %d processed" % (i, j)

def constructDD():
    global dd
    global pp
    drugLength = len(drugs)
    proteinLength = len(proteins)
    pindex = proteins.index
    for d1 in range(len(dd)):
        print d1
        dd[d1][d1] = 1.0
        for d2 in range(d1):
            temp = 0.0
            d1targets = drugs[d1].targets
            d2targets = drugs[d2].targets
            d1length = len(d1targets)
            d2length = len(d2targets)                
            for p1 in d1targets:
                for p2 in d2targets:
                    p1index = pindex(p1)
                    p2index = pindex(p2)
                    temp += pp[p1index][p2index]
                    
            dd[d1][d2] = temp / (d1length * d2length)
            dd[d2][d1] = dd[d1][d2]
    dindex = drugs.index            
    for p1 in range(proteinLength):
        print p1
        pp[p1][p1] = 1.0
        for p2 in range(p1):
            temp = 0.0
            p1drugs = proteins[p1].drugs
            p2drugs = proteins[p2].drugs
            p1length = len(p1drugs)
            p2length = len(p2drugs)
            
            for d1 in p1drugs:
                for d2 in p2drugs:
                    d1index = dindex(d1)
                    d2index = dindex(d2)
                    temp += dd[d1index][d2index]
                    
            pp[p1][p2] = temp / (p1length * p2length)
            pp[p2][p1] = pp[p1][p2]                


def writeDD():
    
    outputFile = open("./result.txt", "w")
    outputFile.write("Drug1\t Drug2\t Correlation")
    drugLength = len(drugs)
    for i in range(drugLength):
        for j in range(i):
            if (dd[i][j] != 0):
                outputFile.write("%s\t %s\t %f" % (drugs[i].name, drugs[j].name, dd[i][j]))
                
    outputFile.close()
    
def parseDEG():
    
    DEGFile = open("./Table_S2.csv")
    
    DEGFile.readline()
    
    for line in DEGFile:
        lineList = line.split(",")
        if lineList[8].strip() != "agent":
            continue
        
        agentA = lineList[10].strip()
        agentB = lineList[13].strip()
        
        drugA = None
        drugB = None
        
        for k in drugs:
            if k.name == agentA:
                drugA = k
            elif k.name == agentB:
                drugB = k
                
            if (drugA != None) and (drugB != None):
                break
            
        

parseDrugTarget()
print "parsing finished"
initPP()
print "init finished"
print len(drugs), len(dd)
for i in range(5):
    constructDD()
    print "%d cycle finished" % i
print "constructing DD finished"
writeDD()