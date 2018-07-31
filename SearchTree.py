import copy
roots=[]
fullDataPointer=0
Class1Size=0
cogsAmount=0
initArrayPointer=0
BOUND=0
activeCellsTemp=0

#@Tree is array of roots. each root represents alpha with two members(cogs)
class Tree:
    def __init__(self, bound, initArray, fullData, Class1SeqsAmount):
        global roots
        roots = []
        global fullDataPointer
        global initArrayPointer
        global BOUND
        global Class1Size
        global cogsAmount
        Class1Size=Class1SeqsAmount
        BOUND=bound
        self.level=0
        self.BOUND=bound
        self.activeCellsAmount=0
        fullDataPointer=fullData
        initArrayPointer=initArray
        cogsAmount=len(initArrayPointer)
        for i in range(0,len(initArrayPointer),1):
            for j in range(i+1, len(initArrayPointer), 1):
                alpha=[i,j]
                relSup1Alpha=calcRelsup1(alpha)
                maxSup2Alpha=initArrayPointer[i][j]
                if ((relSup1Alpha-maxSup2Alpha)>=self.BOUND): #prone condotion
                    coupleNode=Node(alpha,maxSup2Alpha,0,0)
                    roots.append(coupleNode)
                    self.activeCellsAmount+=1

    def addNewLevel(self):
        global roots
        global activeCellsTemp
        activeCellsTemp=0
        for i in range(len(roots)):
            roots[i].expandRoot(self.level)
        self.level+=1
        self.activeCellsAmount=activeCellsTemp

    def getLeafsAlphas(self):
        leafsAlphas=[]
        for i in range(len(roots)):
            rootIleafs=roots[i].findLeafs([])
            for j in range(len(rootIleafs)):
                leafsAlphas.append(rootIleafs[j])
        return leafsAlphas

#@FIELDS: lastActiveInPath: indicaties that all the children of the Node isn't active
#@ active: indicates that this cell needs to be expanded
class Node:
    def __init__(self, alpha, maxSup, depth, father):
        self.depth=depth
        self.father=father
        self.maxSup = maxSup
        self.alpha= alpha
        self.children = []
        self.active=True
        self.lastActiveInPath=True

    def expandRoot(self, level):
        global activeCellsTemp
        if (self.depth==level and self.active==True):
            for i in range(self.alpha[-1]+1,cogsAmount,1):
                newAlpha=copy.deepcopy(self.alpha)
                newAlpha.append(i)
                relsup1Alpha=calcRelsup1(newAlpha)
                maxSup2Alpha=self.calcMaxsup2(newAlpha)
                newNode = Node(newAlpha, maxSup2Alpha, self.depth + 1, self)
                self.children.append(newNode)
                if ((relsup1Alpha-maxSup2Alpha)>=BOUND):
                    newNode.active=True
                    newNode.father.lastActiveInPath=False
                    activeCellsTemp+=1
                else:
                    newNode.active=False #the unic place of Proning (except when building the tree)
                    newNode.lastActiveInPath=False #optional
        else:
            if (len(self.children)>0):
                for i in range(len(self.children)):
                    return self.children[i].expandRoot(level)
            #else: leaf with death!level--> proned leaf. or not active leaf

    #@calcMaxsup2 using the formula of three members(a1,a2,a3) to find calcMaxsup2
    #alpha=(i1...in), n>=3. alpha1=(i1...in-1). alpha2=(i1...in-2,in), alpha3=(in-1,in)
    #self is the father, this functions calculates value MaxSup2 to his optional new child alpha's('newAlpha')
    def calcMaxsup2(self, newAlpha):
        alpha1val=self.maxSup
        alpha3val=initArrayPointer[newAlpha[-2]][newAlpha[-1]]
        if (len(newAlpha) == 3):
            alpha2val=initArrayPointer[newAlpha[0]][newAlpha[2]]
        elif (len(newAlpha)>3):
            grandfather=self.father
            alpha2Node=binarySearch(grandfather.children, newAlpha[-1])
            if alpha2Node!=0:
                alpha2val=alpha2Node.maxSup
            else:
                alpha2val=0
        return max(alpha1val, alpha2val, alpha3val)

    #@findLeafs for output calculation: finds last active in path nodes and active leafs.
    def findLeafs(self, leafList):
        if len(self.children)==0 and self.active==True:
            leafList.append(self.alpha)
            return leafList
        if len(self.children)>0 and self.lastActiveInPath==True:
            leafList.append(self.alpha)
            return leafList
        elif len(self.children)>0 and self.lastActiveInPath==False:
            for i in range(len(self.children)):
                if self.children[i].active==True:
                    self.children[i].findLeafs(leafList)
        return leafList

def calcRelsup1(alpha):
     counter=0
     for i in range(Class1Size):
         flag=1
         for j in range(len(alpha)):
             if (fullDataPointer[i][alpha[j]]!=1):
                 j=len(alpha)
                 flag=0
         if (flag==1):
             counter+=1
     return (counter/Class1Size)

def binarySearch(childs, cog):
    if len(childs) == 0:
        return 0
    else:
        midpoint = len(childs)//2
        if childs[midpoint].alpha[-1]==cog:
            return childs[midpoint]
        else:
          if cog<childs[midpoint].alpha[-1]:
            return binarySearch(childs[:midpoint],cog)
          else:
            return binarySearch(childs[midpoint+1:],cog)

