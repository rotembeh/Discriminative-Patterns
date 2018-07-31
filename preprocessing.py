from __future__ import division
import sys
from itertools import chain, combinations

def cogsShrink(arr,numOfCogs):
    # INPUT: arr is two dimentions array with sequences cogs data, of two classes.
    # cogsShrink makes list (activeCogsList) of cogs in use in our two classes.
    # newArray is array without unused cogs.
    #saves the list (activeCogsList) into CogsFile.txt file
    activeCogsList=[]
    for cog in range(0,numOfCogs,1):
        for seq in range(0,len(arr),1):
            if (arr[seq][cog]!=0):
                activeCogsList.append(cog)
                break;
    activeCogsAmount=len(activeCogsList)
    newArray = [[0 for x in range(activeCogsAmount)] for y in range(len(arr))]
    for seq in range(0, len(arr), 1):
        for cog in range(0,activeCogsAmount,1):
            newArray[seq][cog]=arr[seq][activeCogsList[cog]]
    return activeCogsAmount,activeCogsList,newArray

#preProc calculates maxSup2(alpha) table (MaxSupBasic) for alpha in size 2. (for class 2).
def preProc(arr1,arr2, numOfCogs):
    Class1SeqsAmount=(len(arr1));
    Class2SeqsAmount=(len(arr2));
    for seq in arr2:
        arr1.append(seq);
    activeCogsAmount, activeCogsList, A = cogsShrink(arr1,numOfCogs) #TOTAL COGS OF TWO CLASSES
    MaxSupBasic = [[0 for x in range(activeCogsAmount)] for y in range(activeCogsAmount)]
    pairIterator=combinations(range(activeCogsAmount), 2)
    PairsAmount=(activeCogsAmount*(activeCogsAmount-1))/2;
    k=0
    while (k<PairsAmount):
        pair=pairIterator.next()
        MaxSupBasic[pair[0]][pair[1]]=calcMaxSupBasic(pair[0],pair[1],A, Class1SeqsAmount, Class2SeqsAmount)
        k+=1
    return A, Class1SeqsAmount, Class2SeqsAmount, activeCogsList, MaxSupBasic

def calcMaxSupBasic(cog1,cog2,A, Class2IDX, Class2SeqsAmount):
    pairSupport=0
    for i in range(Class2IDX, len(A),1):
        if (A[i][cog1]==1 & A[i][cog2]==1):
            pairSupport+=1
    if (pairSupport==0):
        return 0
    ans=pairSupport/Class2SeqsAmount
    return ans

def SaveActiveCogsList(list):
    F=open("activeCogsList.txt","w")
    idx=0;
    for cog in list:
        F.write(str(idx));
        F.write(":");
        F.write(str(cog));
        F.write("\n");
        idx+=1

def printMatrix(a):
    for row in a:
        print(' '.join([str(elem) for elem in row]))
    print ""

