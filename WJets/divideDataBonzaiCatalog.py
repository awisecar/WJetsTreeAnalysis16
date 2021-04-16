#! /usr/bin/env python2
import os
import sys
import datetime
import csv
from array import array

print("\nBegin!")

####################################

# Select number of sublists, year
numSubLists = 10
year = 2016
# year = 2017
# year = 2018

####################################

directory = "DataW_txt_" + str(year) + "/"

# Open the data file and read all of the lines
fileInName = directory + "SMu_13TeV_Data_dR_5311_List.txt"
print('Opening: '+fileInName)
fileIn = open('./'+fileInName, "rt")
textAll = fileIn.read()
fileIn.close()

textAll = (textAll.split('\n')) #split full text into a list where each element is separated by newline
del textAll[-1] # delete last element of list (empty line from file read in)

numFiles = len(textAll)
print("\n#files total: "+str(numFiles))
print("#separate lists specified: "+str(numSubLists))

# see how many lines per file we will make
# if we divide by int then the result rounds down (in python2 it always rounds down to the "floor")
# so adding a +1 is basically rounding up
numFilesPerSubList = (numFiles/numSubLists) + 1
print("#files per separate list: "+str(numFilesPerSubList))

# Divide main list into sublists (#sublists is #output files)
contentList = []
for i in range(0, numFiles, numFilesPerSubList):
    contentList.append(textAll[i:i + numFilesPerSubList])
print("\nDivided all filenames into "+str(len(contentList))+" sublists!")

# Write contents of each sublist into separate output files
# Iterate through sublists, then iterate through elements of each sublist
for fileIdx in range(0, numSubLists):
    fileOutName = "SMu_13TeV_Data_dR_5311_List_"+str(fileIdx)+".txt"
    print('\nMaking list: '+fileOutName)
    print('#files to write: '+str(len(contentList[fileIdx])))
    outFile = open(directory + fileOutName, "w") # open output file
    for iLine, line in enumerate(contentList[fileIdx]):
        # print(line)
        outFile.write(line+"\n")
    outFile.close() # close output file

print("\nFinished!")