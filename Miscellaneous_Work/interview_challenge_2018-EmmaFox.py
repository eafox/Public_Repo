############
# APPROACH #
############

#The data given behaves as an absorbing Markov chain.
 
#To find the stable distribution of an absorbing Markov chain (in other
# words, the most likely final distribution of the ore between the terminal states),
# you need to first write up the transition matrix (making sure it is in standard form).
#Then, you calculate the limiting distribution. 
#This can be plugged back into a matrix of the same dimensions as the 
# transition matrix and multiplied by a matrix of the initial state distributions.
#This will then yield the stable distribution of the Markov chain.
#https://en.wikipedia.org/wiki/Absorbing_Markov_chain
#https://www.youtube.com/watch?v=u89Sd514EDI

#An important distinction mentioned in the instructions is that a path
# from the current state to a terminal state will always exist.
#This is a necessary condition for the limiting distribution of the 
# absorbing chain to exist. 

#I am working in the standard form below:
# I | 0
#-------
# R | Q

#Based off the following video: https://www.youtube.com/watch?v=qhnFHnLkrfA 

########
# NOTE #
########

#If running this script from the commandline, the input list of lists
# must be put in quotation marks

###########
# IMPORTS #
###########

import sys, ast
import numpy as np
from numpy.linalg import inv

##########################
# LOAD AND FORMAT MATRIX #
##########################

def Load_matx(m):
    """Take matrix input from commandline and format for analysis"""
    #Convert counts to conversion rates
    m=[[num/sum(row) for num in row] if sum(row)>0 else row for row in m]
    #Convert list of list to matrix
    M=np.matrix(m,dtype=float)
    
    #Find index of unreachable states (will have column sum equal to 0)
    unreach=M.sum(axis=0)
    elimIdx=[x for x in np.where(unreach==0)[1].tolist() if x!=0]
    #Remove unreachable states for the time being
    #(For the analysis to work, all rows must sum to 1 and having an unreachable
    # state will disrupt this)
    M=np.delete(M, elimIdx, axis=1)
    M=np.delete(M, elimIdx, axis=0)
    
    return([M,elimIdx])

def Stand_form(M):
    """Convert matrix of transition probabilities to standard form of
    absorbant Markov chain"""
    #Make template for standard form matrix
    Mstand=np.matrix(np.zeros(M.shape))
    
    #Fill in the I (identity) section of the standard form matrix by 
    # identifying absorbant/terminal states
    absorb=M.sum(axis=1)
    identIdx=[x for x in np.where(absorb==0)[0].tolist() if x!=0]
    for pos,val in enumerate(identIdx):
        Mstand[pos,pos]=1 
    
    #Fill in the R and Q sections of the standard form matrix
    saveIdx=np.where(M.any(axis=1))[0].tolist()
    saveRows=M[saveIdx]
    #New order of axes
    ord_new=identIdx+saveIdx
    #Put new section in with data reflecting new axes order
    Mstand[len(identIdx):]=saveRows[:,ord_new]
    return([Mstand,identIdx])
    #I did consider just returning the R and Q components at this stage
    # instead of the full standard matrix but thought it might be useful
    # to have the standard form accessible at some point.

###################################
# CALCULATE LIMITING DISTRIBUTION #
###################################

def Lim_dist(Mstand,identIdx):
    """Take a standard form of an absorbant Markov chain matrix and
    return the limting distribution matrix"""
    #Get R and Q sections 
    R=Mstand[len(identIdx):,:len(identIdx)]
    Q=Mstand[len(identIdx):,len(identIdx):]
    
    #Calculate the fundamental matrix (F=(I-Q)^-1)
    F=inv(np.identity(len(Q))-Q)
    
    #Calculate new section that will replace R
    newSec=F.dot(R)
    
    #Assemble fundamental matrix
    LDmatx=np.matrix(np.zeros(Mstand.shape))
    #Add in identity section
    for pos,val in enumerate(identIdx):
        LDmatx[pos,pos]=1 
    #Add in product of F and R
    LDmatx[len(identIdx):,:len(identIdx)]=newSec
    
    return(LDmatx)

#######################################
# CALCULATE STABLE/FINAL DISTRIBUTION #
#######################################

def Fin_dist(LDmatx,elimIdx,identIdx):
    """Create an initial state matrix (reflecting all ore in state 0
    to calculate the final distribution between terminal states.
    Format answer to reflect question requirements"""
    #Create matrix for initial state with 100% of the ore being in 
    # state 0 (wherever that is located on the new axes)
    initState=np.matrix(np.zeros(len(LDmatx)))
    initState[0,len(identIdx)]=1
    
    #Solve dot product of initial and limiting distribution matrices
    endDist=initState.dot(LDmatx).tolist()[0]
    #Put states back in order
    endDist_ord=endDist[len(identIdx):]+endDist[:len(identIdx)]
    #Add in any unreachable states
    if elimIdx:
        for i in elimIdx:
            endDist_ord.insert(i,0.)
    
    #Get rid of non-terminal states to give answer
    ans=[]
    for pos,val in enumerate(endDist_ord):
        if val!=0 or pos in elimIdx:
            ans.append(val)
            
    
    return(ans)

##################    
# FULL EXECUTION #
##################

def main(listOfLists):
    """Runs all functions to return the final distributions from the
    input list of lists"""
    M1, eIdx=Load_matx(listOfLists)
    M2, iIdx=Stand_form(M1)
    M3=Lim_dist(M2,iIdx)
    M4=Fin_dist(M3,eIdx,iIdx)

    print("The final predicted proportion of ore in each terminal state is:")
    print(M4)
    
    return(M4)
    

###############################
# TAKE INPUT FROM COMMANDLINE #
###############################

if __name__ == "__main__":
    #Run the input if given
    if len(sys.argv)>1:
        main(ast.literal_eval(sys.argv[1]))
        
    #If no input given, run the test cases
    else:
        print("Distribution for test case 1:")
        main([[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]])

        print("Distribution for test case 2:")
        main([[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])


################
# RUN EXAMPLES #
################

#This section contains the two test cases and code to run through them
# step by step, examining the matrix created after each step.
#The final 2 lines will run the full pipeline on each sample.

#~ testr=[[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]] 
#~ testr=[[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]  

#~ M1, eIdx=Load_matx(testr)
#~ M1
#~ M2, iIdx=Stand_form(M1)
#~ M2
#~ M3=Lim_dist(M2,iIdx)
#~ M3
#~ M4=Fin_dist(M3,eIdx,iIdx)
#~ M4

