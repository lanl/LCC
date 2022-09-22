#!/usr/bin/env python3
#Example: ./getOtrthoCongruent.py -n0 1 1 1 -lattice "latticeFile.dat"
#
import numpy as np
import argparse
from orthoMods import * 

#Parse input files 
parser = argparse.ArgumentParser(description="""Construct two orthogonal Vectors with a periodicity congruent with the lattice""")
parser.add_argument("-n0", help="Passing given vector", type=float, nargs="+")
parser.add_argument("-lattice", help="Lattice filename", type=str )
parser.add_argument("-tol", help="Tolerance", type=float )
options = parser.parse_args()
if(len(options.n0) != 3):
    print("ERROR: Three floats are needed as input")
    exit(0)

#Read the lattice
myFile = open(options.lattice,"r")
i = -1
matM = np.zeros((3,3))
for lines in myFile:
        lines_split = lines.split()
        if(len(lines_split) >= 1):
            i = i + 1
            matM[0,i] = float(lines_split[0])
            matM[1,i] = float(lines_split[1])
            matM[2,i] = float(lines_split[2])

print("My Lattice: \n",matM)

n0 = np.zeros((3))
n0[:] = options.n0[:]
tol = options.tol

print("My Given vector: \n",n0)

#Search for the given congruent equivalent vector.
#First compute M^-1:
invM = np.linalg.inv(matM)
print("My Lattice inverse: \n",invM)

#We now multiply M^-1 with our given vector:
n0TimesInvM = np.dot(invM,n0)
print("Pseudo integers : \n",n0TimesInvM)

#Do a loop to find the "l" so that l*n0 can be written as 
#an integer (a,b,c) linear combination of (a1,a2,a3).
l, ln0, error0 = scan_error(n0TimesInvM,0.001,100000,tol,"n0.err",True)
print("Resulting vector n0=",ln0)
print("Resulting scaling factor l=",l)
print("Resulting error =",error0)

#Scan for orthogonal vectors. We do not care about the 
#orientation as log as they are all orthogonal to each other
myFile = open("n1n2.err","w")
bestError = 3.0
n = 100
for i in range(n):
    epsilon = 2*float(i)/float(n) #Epsilon is a parameter that will change the orientation
    n1, n2 = get_two_orthogonal_to(n0,epsilon)
    
    n1TimesInvM = np.dot(invM,n1)
    l1, ln1, error1 = scan_error(n1TimesInvM,0.001,100000,tol,"n1.err",False)

    n2TimesInvM = np.dot(invM,n2)
    l2, ln2, error2 = scan_error(n2TimesInvM,0.001,100000,tol,"n2.err",False)

    totError = error1 + error2
    print(" ")
    print("ITER,Epsilon,Total Error,l1,l2",i,epsilon,totError,l1,l2)
    string = str(epsilon)+" "+str(totError)+" "+str(l1)+" "+str(l2)+"\n"
    myFile.write(string)
    if(totError < 2*tol):
        break
    else:
        if(totError < bestError):
            #Best vals so far ...
            bestError = totError 
            bestError1 = error1
            bestError2 = error2
            bestl1  = l1
            bestl2  = l2
            bestN1 = n1
            bestN2 = n2
    if(i == n-1):
        totError = bestError 
        error1 = bestError1
        error2 = bestError2
        l1 = bestl1
        l2 = bestl2
        n1 = bestN1
        n2 = bestN2

print(" ")
print("=========================================")
print("             Final results               ")
print("=========================================\n")
print("Given lattice vectors:")
print("a1 =",matM[0,:])
print("a2 =",matM[1,:])
print("a3 =",matM[2,:],"\n")
print("Given vector with best elongation")
print("n0 =",n0[0],n0[1],n0[2],", l0 =",l,"\n")
print("Best two orthogonal vectors find by the algorithm with best elongation")
print("n1 =",n1[0],n1[1],n1[2],", l1 =",l1)
print("n2 =",n2[0],n2[1],n2[2],", l2 =",l2,"\n")
print("Translation errors")
print("Err0 =",error0)
print("Err1 =",error1)
print("Err2 =",error2)




