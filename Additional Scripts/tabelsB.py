import pickle

intToChar = {1:"A", 2:"B", 3:"C", 4:"D", 5:"E", 6:"F", 7:"G", 8:"H"}

digDict = {1: "foregut", 2:"hindgut", 3:"monogastric"}

DigesDictNums = {1: 1, 3: 1, 12: 1, 4: 1, 6: 1, 9: 1, 13: 1, 14: 1,15:1, 16: 1, 19: 1, 11: 2,
                31: 2, 32: 2, 33: 2, 34:2,35:2, 2: 3, 29: 3, 8: 3, 20: 3, 21: 3, 27: 3, 23: 3,
                24: 3, 17: 3, 18: 3, 26: 3, 36:3, 'MM': 3, 22: 3, 25: 3, 5: 3, 7: 3}

phyloDictNums = {1: 2, 3: 2, 12: 2, 4: 3, 6: 3, 9: 3, 13: 3, 14: 3, 16: 3, 19: 3, 11: 1,
                31: 1, 32: 1, 33: 1, 34: 1, 2: 4, 29: 4, 8: 7, 20: 7, 21: 7, 27: 7, 23: 7,
                24: 7, 17: 6, 18: 6, 26: 6, 36: 6, 'MM': 6, 22: 5, 25: 5, 5: 8, 7: 8}

##values: amount of cogs
parameters = {
              "cog30": [23190+1],
              "operon30": [38522+1],

              "cog40": [28965+1],
              "operon40": [29092+1],

              "cog50": [33057+1],
              "operon50": [23190+1],

              "cog60": [36780+1],
              "operon60": [18558+1],

              "cog80": [42971+1],
              "operon80": [11540+1],

              "cog90": [45424+1],
              "operon90": [9053+1],
              } #maximum operon


def makeTables(class1, class2, taskID, param):
    if (taskID == 1):
        return dig_cog(class1,class2, param)
    if (taskID == 2):
        return phylo_cog(class1, class2, param)
    if (taskID == 3):
        return dig_operon(class1, class2, param)
    if (taskID == 4):
        return phylo_operon(class1, class2, param)

def dig_cog(class1,class2, param):
    array1 = makeAnimalCogMatrix(class1, param)
    array2 = makeAnimalCogMatrix(class2, param)
    name1 = digDict[class1]
    name2 = digDict[class2]
    cogsAmount = parameters["cog"+param][0]
    return array1, array2, name1, name2, cogsAmount


def dig_operon(class1, class2, param):
    array1 = makeAnimalOperonMatrix(class1, param)
    array2 = makeAnimalOperonMatrix(class2, param)
    name1 = digDict[class1]
    name2 = digDict[class2]
    cogsAmount = parameters["operon"+param][0]
    return array1, array2, name1, name2, cogsAmount


def phylo_cog(class1, class2, param):
    array1 = makePhyloCogMatrix(class1, param)
    array2 = makePhyloCogMatrix(class2, param)
    name1 = intToChar[class1]
    name2 = intToChar[class2]
    cogsAmount = parameters["cog"+param][0]
    return array1, array2, name1, name2, cogsAmount


def phylo_operon(class1, class2, param):
    array1 = makePhyloOperonMatrix(class1, param)
    array2 = makePhyloOperonMatrix(class2, param)
    name1 = intToChar[class1]
    name2 = intToChar[class2]
    cogsAmount = parameters["operon"+param][0]
    return array1, array2, name1, name2, cogsAmount


def makePhyloOperonMatrix(classI, param):
    newArray = []
    indexesDict = {}
    file = open("data/OGB_catalog_plasmidome_safari_"+param+"_ins10_q1_0_q2_5_l2_instances.fasta", "r")
    nextidx = 0
    for line in file:
        if line[0] == ">":
            line1 = line.strip(">")
            line1=line1.split()
            operon = line1[0]
        else:
            line1 = line.split()
            line1 = line1[1].split("_")
            animal = line1[0]
            if animal != "MM":
                if int(animal) in phyloDictNums.keys():
                    if (phyloDictNums[int(animal)] == classI):
                        if not animal in indexesDict.keys():
                            indexesDict[animal] = nextidx
                            newArray.append([0 for i in range(parameters["operon"+param][0])])
                            nextidx += 1
                        newArray[indexesDict[animal]][int(operon)] = 1
              #          print("marking operon ",operon,"in animal ",animal," from class ",classI)
            else: ##its 'MM' line
                if (phyloDictNums[animal] == classI):
                    if not animal in indexesDict.keys():
                        indexesDict[animal] = nextidx
                        newArray.append([0 for i in range(parameters["operon" + param][0])])
                        nextidx += 1
                    newArray[indexesDict[animal]][int(operon)] = 1
    return newArray



##linesType: line is animal
def makePhyloCogMatrix(classI, param):
    newArray = []
    indexesDict = {}
    file = open("data/geneteams_"+param+".txt", "r")
    nextidx = 0
    for line in file:
        line1 = line.strip().split()
        plasmid=line1[0].strip().split("_")
        if plasmid[0] != 'MM':
            if int(plasmid[0]) in phyloDictNums.keys():
                if (phyloDictNums[int(plasmid[0])] == classI):
                    for cog in line1[1:]:
                        if cog!='X':
                            if not plasmid[0] in indexesDict.keys():
                                indexesDict[plasmid[0]] = nextidx
                                newArray.append([0 for i in range(parameters["cog"+param][0])])
                                nextidx += 1
                            newArray[indexesDict[plasmid[0]]][int(cog)] = 1
        else: ##its 'MM' line
            if (phyloDictNums[plasmid[0]] == classI):
                for cog in line1[1:]:
                    if cog != 'X': ##MM is tha last idx
                        if not plasmid[0] in indexesDict.keys():
                            indexesDict[plasmid[0]] = nextidx
                            newArray.append([0 for i in range(parameters["cog" + param][0])])
                            nextidx += 1
                        newArray[indexesDict[plasmid[0]]][int(cog)] = 1
    return newArray



def makeAnimalOperonMatrix(classI, param):
        newArray = []
        indexesDict = {}
        file = open("data/OGB_catalog_plasmidome_safari_" + param + "_ins10_q1_0_q2_5_l2_instances.fasta", "r")
        nextidx = 0
        for line in file:
            if line[0] == ">":
                line1 = line.strip(">")
                line1 = line1.split()
                operon = line1[0]
            else:
                line1 = line.split()
                line1 = line1[1].split("_")
                animal = line1[0]
                if animal != "MM":
                    if int(animal) in DigesDictNums.keys():
                        if (DigesDictNums[int(animal)] == classI):
                            if not animal in indexesDict.keys():
                                indexesDict[animal] = nextidx
                                newArray.append([0 for i in range(parameters["operon" + param][0])])
                                nextidx += 1
                            newArray[indexesDict[animal]][int(operon)] = 1
                            #          print("marking operon ",operon,"in animal ",animal," from class ",classI)
                else:  ##its 'MM' line
                    if (DigesDictNums[animal] == classI):
                        if not animal in indexesDict.keys():
                            indexesDict[animal] = nextidx
                            newArray.append([0 for i in range(parameters["operon" + param][0])])
                            nextidx += 1
                        newArray[indexesDict[animal]][int(operon)] = 1
        return newArray

##linesType: line is animal
def makeAnimalCogMatrix(classI, param):
    newArray = []
    indexesDict = {}
    file = open("data/geneteams_"+param+".txt", "r")
    nextidx = 0
    for line in file:
        line1 = line.strip().split()
        plasmid=line1[0].strip().split("_")
        if plasmid[0] != 'MM':
            if int(plasmid[0]) in DigesDictNums.keys():
                if (DigesDictNums[int(plasmid[0])] == classI):
                    for cog in line1[1:]:
                        if cog!='X':
                            if not plasmid[0] in indexesDict.keys():
                                indexesDict[plasmid[0]] = nextidx
                                newArray.append([0 for i in range(parameters["cog"+param][0])])
                                nextidx += 1
                    #        print("animal: " + plasmid[0] + " cog: " + cog)
                            newArray[indexesDict[plasmid[0]]][int(cog)] = 1
        else: ##its 'MM' line
            if (DigesDictNums[plasmid[0]] == classI):
                for cog in line1[1:]:
                    if cog != 'X': ##MM is tha last idx
                        if not plasmid[0] in indexesDict.keys():
                            indexesDict[plasmid[0]] = nextidx
                            newArray.append([0 for i in range(parameters["cog" + param][0])])
                            nextidx += 1
                    #    print("animal: " + plasmid[0] + " cog: " + int(cog))
                        newArray[indexesDict[plasmid[0]]][int(cog)] = 1
    return newArray



def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)




