#!/usr/bin/env python3
import numpy as np
import argparse

#Get two orthogonal vectors based on a free parameter epsilon
# vN: Given vector
# epsilon: Parameter to change the orientation
# vN1: First output vector
# vN2: Second output vector
def get_two_orthogonal_to(vN,epsilon):
    #Given and resulting vectors allocation
    vN1 = np.zeros((3))
    vN2 = np.zeros((3))
   
    k1 = 1.0
    l2 = k1
    l1 = epsilon*k1

    h = vN[0]
    k = vN[1]
    l = vN[2]

    #Loop to track the indexing needed to 
    #permute the entries of the given vector 
    #such that the first entry is always non-zero
    indexTrack = np.zeros((3),dtype=int)
    indexBackTrack = np.zeros((3),dtype=int)
    if(abs(vN[0]) > 0.0):
        indexTrack[0] = 0
        indexTrack[1] = 1
        indexTrack[2] = 2
        indexBackTrack[0] = 0
        indexBackTrack[1] = 1
        indexBackTrack[2] = 2
    else:
        if(abs(vN[1]) > 0.0):
            indexTrack[1] = 0
            indexTrack[0] = 1
            indexTrack[2] = 2
            indexBackTrack[0] = 1
            indexBackTrack[1] = 0
            indexBackTrack[2] = 2
        else:
            if(abs(vN[2]) > 0.0):
                indexTrack[2] = 0
                indexTrack[0] = 2
                indexTrack[1] = 1
                indexBackTrack[0] = 2
                indexBackTrack[1] = 1
                indexBackTrack[2] = 0

    #Getting permutted h k l values
    h = vN[indexTrack[0]]
    k = vN[indexTrack[1]]
    l = vN[indexTrack[2]]

    h1 = (-k*k1 -l*l1)/h
    k2 = (-l1*l2 - k*k1*l*l2/h**2 - l*l*l1*l2/h**2)/(k*k*k1/h**2 + l*l1*k/h**2 + k1)
    h2 = (-k*k2 -l*l2)/h

    n1 = np.zeros((3))
    n2 = np.zeros((3))
    #Normalize
    n1[0] = h1/np.sqrt(h1**2 + k1**2 + l1**2)
    n1[1] = k1/np.sqrt(h1**2 + k1**2 + l1**2)
    n1[2] = l1/np.sqrt(h1**2 + k1**2 + l1**2)
    n2[0] = h2/np.sqrt(h2**2 + k2**2 + l2**2)
    n2[1] = k2/np.sqrt(h2**2 + k2**2 + l2**2)
    n2[2] = l2/np.sqrt(h2**2 + k2**2 + l2**2)

    #Do the inverse of the permutations
    vN1[indexBackTrack[0]] = n1[0]
    vN1[indexBackTrack[1]] = n1[1]
    vN1[indexBackTrack[2]] = n1[2]
    vN2[indexBackTrack[0]] = n2[0]
    vN2[indexBackTrack[1]] = n2[1]
    vN2[indexBackTrack[2]] = n2[2]

    print("\nGiven initial vector:")
    print("   ",vN[:])
    print("\nOrthogonal choice:")
    print("   ",vN1[:])
    print("   ",vN2[:])

    sumDot = np.dot(vN,vN1) + np.dot(vN1,vN2) + np.dot(vN2,vN)
    print("\nSum of the two-by-two dot prods=",sumDot)
    if(abs(sumDot) > 1E-10):
        print("ERROR: Vectors are not orthogonal")

    return vN1,vN2

#Scan for the smallest elongation factor l
#vect: Vector to be enlarged
#dl: Step to search for l
#n: Number of steps
#tol: Tolerance in the error
#outFileName: Name of the output file to store (l,error)
def scan_error(vect,dl,n,tol,outFileName,writel):


    if(writel == True):
        myFile = open(outFileName,"w")
        myFile.write("#Error for n0\n")

    minError = 1.5
    for i in range(n):
        l = 0.5 + i*dl #A hundredth of and angstrom
        errVect = l*vect
        aux1 = min(abs(errVect[0] - np.floor(errVect[0])),\
            abs(errVect[0] - np.ceil(errVect[0])))
        aux2 = min(abs(errVect[1] - np.floor(errVect[1])),\
            abs(errVect[1] - np.ceil(errVect[1])))
        aux3 = min(abs(errVect[2] - np.floor(errVect[2])),\
            abs(errVect[2] - np.ceil(errVect[2])))
        error = aux1 + aux2 + aux3
        string = str(l)+" "+str(error)+"\n"
        if(writel == True):
            myFile.write(string)
        if(error < tol):
            break
        else:
            if(error < minError):
                minError = error
                lMin = l
        if(i == n-1):
            print("Maxiter reached in scan loop. I will use the l that minimizes the error")
            l = lMin
            error = minError
            
    lTimesVect = l*vect
    
    return l,lTimesVect,error






