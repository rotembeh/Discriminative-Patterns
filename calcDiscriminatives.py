from __future__ import division
from preprocessing import preProc
from SearchTree import Tree
import copy
import os
import pickle

##finds discriminative patterns between class1Matrix and class2Matrix.
def calcDisc(class1matrix, class2matrix, taskName, cogsPrecent, class1, class2, cogsAmount):
    print("Starting CalcDisc:", class1, class2, taskName, cogsPrecent)
    if(taskName == "digestion-Cogs" or taskName == "digestion-Operons"): ##parameteres for digestion tasks
        pivot = 3 ##will seperate the results to two files, big patterns- with more than pivot cogs, and small patterns file.  
        bound = 0.8 ##maxSup threshold
    else: ## parameters for phylogenetic tasks
        pivot = 4
        bound = 0.8
    taskName = taskName + "-" + cogsPrecent
    arr1 = class1matrix
    arr2 = class2matrix
    print "Running Preprocessing..."
    if (len(arr2)>len(arr1)): ##choosing the smaller class as class1 (not mandatory)
        fullData, Class1SeqsAmount, Class2SeqsAmount, originalCogsIdxList, initArray = preProc(arr1, arr2, cogsAmount)
    else:
        temp = class1 ##swiching names
        class1 = class2
        class2 = temp
        fullData, Class1SeqsAmount, Class2SeqsAmount, originalCogsIdxList, initArray=preProc(arr2, arr1, cogsAmount)
    print "Preprocessing Done"
    print "Building Tree with BOUND=",bound,"..."
    myTree=Tree(bound, initArray, fullData, Class1SeqsAmount)
    print "Buiding Tree Done"
    print "Active Cells: ", myTree.activeCellsAmount
    print "************************"
    print "Adding levels..."
    while(myTree.activeCellsAmount>0):
        print "ADDING NEW LEVEL..."
        myTree.addNewLevel()
        print "LEVEL ", myTree.level, " ADDED"
        print "ActiveCellsAmount: ", myTree.activeCellsAmount
    print "BUILDING LEVELS DONE, END LEVEL: ", myTree.level
    if not os.path.exists("Results"):
        os.makedirs("Results")
    if not os.path.exists("Results/"+taskName):
        os.makedirs("Results/"+taskName)
    if not os.path.exists("Results/"+taskName+"/"+class1+"-vs-"+class2):
        os.makedirs("Results/"+taskName+"/"+class1+"-vs-"+class2)
    file=open("Results/"+taskName+"/"+class1+"-vs-"+class2+"/Discriminatives-"+str(bound)+".txt","w")
    longpatts=open("Results/"+taskName+"/"+class1+"-vs-"+class2+"/Long(at_least_"+str(pivot+1)+")DiscriminativePatterns-"+str(bound)+".txt", "w")
    shortpatts=open("Results/"+taskName+"/"+class1+"-vs-"+class2+"/Short(at_most_"+str(pivot)+")DiscriminativePatterns-"+str(bound)+".txt", "w")
    print "Processing results..."
    outputNodes=myTree.getLeafsAlphas()
    results= convertCogs(outputNodes,originalCogsIdxList)
    print "Converting Done. Writing into text files... "
    file.write(str(results))
    longPatterns=[]
    shortPatterns=[]
    separateResluts(results,longPatterns,shortPatterns,pivot)
    longpatts.write(str(longPatterns))
    shortpatts.write(str(shortPatterns))
    print "DONE."
    print "************************"
    print "BOUND=",bound," Class1Size=", Class1SeqsAmount, "Class2Size=", Class2SeqsAmount, ""
    print "Results:"
    print "TOTAL Discriminative Patterns: ",len(results)
    print "Long Patterns (>",pivot,"): ", len(longPatterns)
    print "Short Patterns(<=",pivot,"): ", len(shortPatterns)
    print "-------------------------"
    print "-------------------------"


#@making arrays from input data files of one class
def makeInitialArrayFromFile(file):
    array=file.read()
    array=array.strip('[')
    array=array.strip(']')
    array=array.split(']')
    array[0]=array[0].split(", ")
    array[0] = list(map(int, array[0]))
    for i in range(1,len(array),1):
        array[i] = array[i].strip(", [")
        array[i] = array[i].split(", ")
        array[i] = list(map(int, array[i]))
    return array

def printMatrix(a):
    for row in a:
        print(' '.join([str(elem) for elem in row]))
    print ""

#@converts calculataed cogs to the originals, using originalCogsIdxList
def convertCogs(discUnconverted, originalCogsIdxList):
    discList=copy.deepcopy(discUnconverted)
    for i in range(len(discUnconverted)):
        for j in range(len(discUnconverted[i])):
            discList[i][j]=originalCogsIdxList[discUnconverted[i][j]]
    return discList

#@separate to long and short discriminative patterns
def separateResluts(total, longlist, shortlist, pivot):
    for i in range(len(total)):
        if len(total[i])>pivot:
            longlist.append(total[i])
        else:
            shortlist.append(total[i])

def openPickle(path):
    with open(path, 'rb') as input:
        matrix = pickle.load(input)
    return matrix