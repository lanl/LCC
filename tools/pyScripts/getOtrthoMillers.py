#!/usr/bin/env python3
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="""Construct two orthogonal Miller indices""")

parser.add_argument("-hkl", help="Passing hkl indices", type=int, nargs="+")

options = parser.parse_args()
if(len(options.hkl) != 3):
    print("ERROR: Three integers are needed as input")
    exit(0)

#Given and resulting vectors allocation
vN = np.zeros((3))
vN1 = np.zeros((3))
vN2 = np.zeros((3))

vN[0] = options.hkl[0]
vN[1] = options.hkl[1]
vN[2] = options.hkl[2]

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

ave = (abs(vN[0]) + abs(vN[1]) + abs(vN[2]))/3
intAve = np.ceil(ave)

#Getting permutted h k l values
h = vN[indexTrack[0]]
k = vN[indexTrack[1]]
l = vN[indexTrack[2]]

#Pick the values for the free parameters 
k1 = 1
l1 = 0
l2 = 1

#Compute the rest of the values to fullfil orthogonality condition
#h1 = (-k*k1 -l*l1)/h
#k2 = (-l1*l2 - k*k1*l*l2/h**2 - l*l*l1*l2/h**2)/(k*k*k1/h**2 + l*l1*k/h**2 + k1)
#h2 = (-k*k2 -l*l2)/h

h1 = -k/h
k2 = -(k*l/h**2)/(k**2/h**2 + 1)
h2 = (l*k**2/h**3)/(k**2/h**2 + 1) - l/h 

#Get the min non-zero of all the entries
myMin = 1000
if(abs(h1) > 0):
    myMin = min(myMin,abs(h1))
if(abs(k1) > 0):
    myMin = min(myMin,abs(k1))
if(abs(l1) > 0):
    myMin = min(myMin,abs(l1))
if(abs(h2) > 0):
    myMin = min(myMin,abs(h2))
if(abs(k2) > 0):
    myMin = min(myMin,abs(k2))
if(abs(l2) > 0):
    myMin = min(myMin,abs(l2))

#Do the inverse of the permutations
#and normalize
vN1[indexBackTrack[0]] = h1/myMin 
vN1[indexBackTrack[1]] = k1/myMin
vN1[indexBackTrack[2]] = l1/myMin
vN2[indexBackTrack[0]] = h2/myMin
vN2[indexBackTrack[1]] = k2/myMin 
vN2[indexBackTrack[2]] = l2/myMin

print("\nGiven initial Miller vector:")
print("   ",vN[:])
print("\nOrthogonal choice:")
print("   ",vN1[:])
print("   ",vN2[:])

sumDot = np.dot(vN,vN1) + np.dot(vN1,vN2) + np.dot(vN2,vN)
print("\nSum of the two-by-two dot prods=",sumDot)
