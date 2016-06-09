def parseDrugTarget():
    
    drugNameIDFile = open("./drug_target_uniprot_links.txt")
    
    drugNameIDFile.readline()
    drugs = []
    proteins = []
    for line in drugNameIDFile:
        add = True
        addP = True
        lineList = line.split("\t")
        newName = lineList[1].strip()
        if newName[0] == '\"':
            newName = newName[1:-1]
        newID = lineList[0].strip()
        newDrug = (newID, newName)
        newProtein = lineList[3].strip()
        for i in drugs:
            if (newDrug == i):
                add = False
                
        if add:
            drugs.append(newDrug)
        
           
        for j in proteins:
            if (newProtein == j):
                addP = False
        if addP:
            proteins.append(newProtein)
            
    drugNameIDFile.close()
        
    outputDrugFile = open("./drugs.txt", "w")
    for i in drugs:
        outputDrugFile.write("%s\t%s\n" % i)
        
    outputDrugFile.close()
    
    outputProteinFile = open("./proteins.txt", "w")
    for j in proteins:
        outputProteinFile.write("%s\n" % j)
        
        
parseDrugTarget()