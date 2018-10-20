from calcDiscriminatives import calcDisc
from tabelsB import makeTables

taskID={"digestion-Cogs":1, "phylogenetic-Cogs":2, "digestion-Operons":3, "phylogenetic-Operons":4}
cogsPrecent = ["30", "40", "50", "60", "80", "90"]



def main():
    print("digestion-Cogs")
    sendParameters("digestion-Cogs", 3)
    print("phylogenetic-Cogs")
    sendParameters("phylogenetic-Cogs", 8)
    print("digestion-Operons")
    sendParameters("digestion-Operons", 3)
    print("phylogenetic-Operons")
    sendParameters("phylogenetic-Operons", 8)

def sendParameters(taskName, classesAmount):
    classesAmount += 1
    for param in cogsPrecent:
        for i in range(1, classesAmount, 1):
            for j in range (i+1, classesAmount, 1):
                    class1matrix, class2matrix, class1name, class2name, cogsAmount = makeTables(i, j, taskID[taskName], param)
                    if (len(class1matrix) > 0 and len(class2matrix) > 0):
                        calcDisc(class1matrix, class2matrix, taskName, param, class1name, class2name, cogsAmount)
                    else:
                        print("skipping classes ", class1name, class2name)

if __name__ == "__main__":
    main();
