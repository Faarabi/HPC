# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 11:18:59 2020

@author: fatemeh
"""


import numpy as np
import itertools
from sympy import Matrix, S, linsolve, symbols, var
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
import timeit

start = timeit.default_timer()


# Class of Atom, each atom has a name which is an int and a rate which is Lambda or arrival rate for that atom
class Atom:
    def __init__(mysillyobject, name, arrivalrate, degradationRate, preferenceList, location):
        mysillyobject.name = name
        mysillyobject.arrivalrate = arrivalrate
        mysillyobject.degradationRate=degradationRate
        mysillyobject.preferenceList=preferenceList
        mysillyobject.location=location
    def ownList(self):
        x=int(self.name)
        return(self.preferenceList[x-1])
        
class Servers:
    def __init__(self, name, preferenceList, currentStatus, servicerate, numberofServers , location):
        self.name=name
        self.preferenceList=preferenceList
        self.currentStatus=currentStatus
        self.servicerate=servicerate
        self.numberofServers=numberofServers
        self.location=location
    # newList: It is a list contains of all the possible status that a server can get
    #newDict: It is a dictionary that assign a number to each possible status
    #e.g. one status is [1,0] it means that server 1 is in the first step of serving 
    #atom 0 and it is the second priority for atom 0 if it was third priority it was like[2,0]
    def serverStatus(self):
        newList=[0,1]
        newDict={0:0,1:1}
        m=2
        for i in range(0, len(self.preferenceList)):
            for j in range(1, self.numberofServers):
                if self.preferenceList[i][j]==self.name:
                    for cal in range (1, j+1):
                        newDict[m]=[(j,cal),i]
                        #newDict={m+cal : [(j,cal),cal] for cal in range(1, j+1) }
                        newList.append([(j,cal),i])
                        m=m+1
        numberOfPossibleStatus=len(newList)
        return(newList, numberOfPossibleStatus, newDict)

    #It gives  a dictionary of rates for each status we will go
    def RatesOut(self):  
        rateDict={}    
        for key in self.serverStatus()[2]:
            if key==0:
                A={}
                for i in range(2, len(self.serverStatus()[2])):
                    if self.serverStatus()[2][i][0][1]==1:
                        A[i]=atom[self.serverStatus()[2][i][1]].arrivalrate
                A[1]= sum(atom[i].arrivalrate for i in range(len(self.preferenceList)) if self.preferenceList[i][0]==self.name)
                rateDict[key]=A
            elif key==1:
                rateDict[key]={0:self.servicerate}
            else:
                
                first=self.serverStatus()[2][key][0][0]
                second=self.serverStatus()[2][key][0][1]
                difference = first-second
                if (second< first):
                    dictWithoutZeroandOne=dict(self.serverStatus()[2])
                    del dictWithoutZeroandOne[0]
                    del dictWithoutZeroandOne[1]
                    a=[ key for key, value in dictWithoutZeroandOne.items()  if value == [(first, second+1), dictWithoutZeroandOne[key][1] ] ]
                    rateDict[key]={a[0]: self.servicerate+self.servicerate/((self.numberofServers-difference)*atom[self.serverStatus()[2][key][0][1]].degradationRate)}
                
                    
                if (second == first):
                    rateDict[key]={1:self.servicerate+self.servicerate/(self.numberofServers*atom[self.serverStatus()[2][key][0][1]].degradationRate)}
                              
        return(rateDict)
    
    


def makeMatrix(obj1, obj2):
    matrixSize=obj1.serverStatus()[1]*obj2.serverStatus()[1]
    myMatrix= np.zeros(shape=(matrixSize, matrixSize))
    myList=list(itertools.product(list(obj1.serverStatus()[2].keys()), list(obj2.serverStatus()[2].keys())))

    for i in range(0,matrixSize):
        From=myList[i]
        To=list(From)

        server1=Servers(obj1.name, obj1.preferenceList, From[0], obj1.servicerate, obj1.numberofServers, obj1.location ) 
        server2=Servers(obj2.name, obj1.preferenceList, From[1], obj2.servicerate, obj2.numberofServers, obj2.location)
        
        servers=[server1, server2]
        for k in range(0, obj1.numberofServers):
            #We want to ignor status 0 and 1 which do not have a list, for example status 2 has a list of [(2, 1), 1] so 
            #we consider the ones with len of list >1
            #It is for status 0 that can go to more than 1 statuses
            if len(servers[k].RatesOut()[servers[k].currentStatus])>1: #we come to this condition just if currentstatus = 0 
                
                for j in servers[k].RatesOut()[servers[k].currentStatus].keys() :
                    To[k]= j
                    #To[k]=i,i is the name of status, the tatus is [(-,-),x] we are looking to find atom name which is x
                    if To[k]==1:
                        myMatrix[i,myList.index(tuple(To))]=servers[k].RatesOut()[From[k]][To[k]]
                        To=list(From)
                    else:   
                        whichAtomToGo = servers[k].serverStatus()[2][To[k]][1] 
                        
                        #If all the servers preferer to server k for the atom be busy then system can go to state To[k]
                        ##if atom[whichAtomToGo].ownList().index(k+1)!=0:
                        if all(From[mm-1] for mm in (atom[whichAtomToGo].ownList()[0:atom[whichAtomToGo].ownList().index(k+1)] ))==True :
                            myMatrix[i, myList.index(tuple(To))]=servers[k].RatesOut()[From[k]][To[k]]
                            To=list(From)
                                
            else: #we come to this condition just if currentstatus be any thing but 0
                
                To[k]= [ii for ii in (servers[k].RatesOut()[servers[k].currentStatus].keys())][0]
                myMatrix[i,myList.index(tuple(To))]=servers[k].RatesOut()[From[k]][To[k]]
                To=list(From)
        #Subtract the sum of all items in row from the diagonal item (The sum of each row should be zero)---> -5 , 2, 3         
    for index in range(matrixSize):
       myMatrix[index,index]=-np.sum(myMatrix, axis=1)[index]  
          
    return(myMatrix)


def Size(obj1,obj2):
    matrixSize=obj1.serverStatus()[1]*obj2.serverStatus()[1]
    return(matrixSize)

###############################################################################

#--------------------------#2 servers----------------------------
pList=[[1,2],[1,2],[1,2],[1,2]]#,[1,3,2],[3,2,1]]
numberOfServers=2
currentStatusList=[0,0]
ServiceRate=[2.5,3.5]
ServerLocations=[(20,10), (20,30)]
#------------------Atoms-----------------------------------------
numberOfAtoms=4
ArrivalRate=[1, 1, 1, 1]#, 1, 1]
DegredationRate=[100, 100, 100, 10]# 0.3, 0.3]#1,1,1]
#AtomLocations=[(1,1),(2,1),(1,2),(2,2),(1,3),(2,3)]
AtomLocations=[(10,10), (20,10), (30,10), (10,20)]#, (20,30), (30,30)]

server=[]
for i in range(0, numberOfServers):
    server.append(Servers( i+1, pList, currentStatusList[i], ServiceRate[i], numberOfServers,ServerLocations[i] ))

atom=[]
for j in range(0, numberOfAtoms):
    atom.append(Atom(j+1,  ArrivalRate[j], DegredationRate[j], pList, AtomLocations[j]))
   
###############################################################################

###############################################################################

np.set_printoptions(threshold=np.inf)
np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})


def Dist(obj,server):
    return(np.linalg.norm(np.array(obj.location)-np.array(server.location)))

############################ compute E_nj and MeanTravelTime ##################
def PP(Input, obj1,obj2,obj3):
    Main=makeMatrix(server[0],server[1]).transpose()
    Main=Main.tolist()
    ll=[]
    bb=[]
    for kk in range(0,12):
        ll.append(1)
        bb.append(0)   
    bb.pop(1)    
    Main.append(ll)
    bb.append(1)
    Main.pop(2)
    A = Matrix(Main)
    bbb = Matrix(bb)
    aa = np.array(Main)
    bbb = np.array(bb)
    A= list(np.linalg.solve(aa,bbb))
    return(A)


def Mean(Input, obj1,obj2):
    for j in range(0,len(atom)):
        atom[j].degradationRate=Input[j]
    myList=list(itertools.product(list(server[0].serverStatus()[2].keys()), list(server[1].serverStatus()[2].keys())))
    E= np.zeros(shape=(len(server),len(atom)))
    for n in range(0, len(server)):
        for j in range(0,len(atom)): 
            #index Of Server N in Atom J
            if pList[j].index(n+1)!=0:
                for i in range (0, len(myList)):
                    if all(myList[i][mm-1] for mm in (pList[j][0:pList[j].index(n+1)] ))==True and myList[i][pList[j].index(n+1)] ==0:
                        E[n][j]+=PP(Input, obj1,obj2)[i]
            else:
                E[n][j]=sum( PP(Input, obj1,obj2)[t]  for t in range(len( myList)) if myList[t][n-1]==0)
    fI=[]      
    for n in range(0, len(server)):        
        fI.append(sum( PP(Input, obj1,obj2)[t]  for t in range(len( myList)) if myList[t][n]!=0))
  
    P_Q=1-sum(PP(Input, obj1, obj2))
    
    P_allserverBusy=sum( PP(Input, obj1, obj2)[i] for i in range(len(myList)) if all(myList[i][x]==1 for x in range(len(server))) == True)
    
    SumOfLambdas=sum((atom[m].arrivalrate) for m in range(len(atom)))
    T_Q=0
    for i in range(len(atom)):
        for j in range(len(atom)):
            T_Q+=(atom[i].arrivalrate/SumOfLambdas) * (atom[j].arrivalrate / SumOfLambdas) *Dist(atom[i], atom[j])
    
    
    for n in range(0, len(server)):
        for j in range(0,len(atom)):  
             temp=E[n][j]* (atom[j].arrivalrate / SumOfLambdas )*Dist(atom[j], server[n]) 
     
    MeanTravelTime=temp+ T_Q  * (P_Q + P_allserverBusy)
    return(MeanTravelTime, fI)







def main():
    Inputs=[100, 70, 10, 1]
    print(PP(Inputs, server[0], server[1]))
    print(Mean(Inputs, server[0], server[1])[0], Mean(Inputs, server[0], server[1])[1])
    print(Size(server[0],server[1]))
    stop = timeit.default_timer()
    print('Time: ', stop - start) 

    
    
    
if __name__=="__main__":
    main()