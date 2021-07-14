from random import randrange
import csv
import sys
import os
from datetime import date
import datetime
import time
from pathlib import Path
import traceback 
import sys
 
BLANK_CHAR = "-"
GAP_OPENING_PENALTY = -3
GAP_EXTENTION_PENALTY = -2
LOOP_COUNT = 30000
INPUT_FOLDER = "input/"
REFERENCE_FOLDER = "reference/"
OUTPUT_FOLDER = "output"
RESULT_FILE_NAME = "result.csv"

distanceMatrixDict = {}
def loadDistanceMatrix():
    with open('distanceMatrix.csv') as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            #if line_count == 0:
            #    print(f'Column names are {", ".join(row)}')
            #    line_count += 1
            #else:
                #print(f'Row data {", ".join(row)}')
            distanceMatrixDict [row[""]] = row
            line_count += 1
        print(f'Processed {line_count} lines.')
        for key1 in distanceMatrixDict.keys():
            print (distanceMatrixDict[key1])
            for key2 in distanceMatrixDict [key1].keys():
                if key2 != "":
                    #print ("ley1 " + key1 + " key2 "+ key2)
                    distanceMatrixDict[key1][key2]= distanceMatrixDict[key2][key1]
        #print (distanceMatrixDict)
        #print ("Testing distance between A and Q " + distanceMatrixDict["A"]["Q"] + " " + distanceMatrixDict["Q"]["A"] )
        #print ("Testing distance between K and L " + distanceMatrixDict["K"]["L"] + " " + distanceMatrixDict["L"]["K"] )

def getTfaDictionary(file):
    myDict={}
    val=""
    key=""
    for x in file:
        if x[0]== ">":
            if val != "":
                myDict [key] = val
                val=""
            key = x.replace("\n", "").replace ("_oo", "");
        else: 
            val=val+x.replace("\n", "").replace (".", BLANK_CHAR)
    if val != "":
        myDict [key] = val  
    return myDict
    

#  Function to generate MSf from the input file. This is the core algorithm. 
#  Need to be modified
def generateMsf (tfaDictionary):
    generatedDictionary ={}
    for key, val in tfaDictionary.items():
        #Padding with dots to make it 96 char long
        paddedVal = val + BLANK_CHAR*(96 - len(val))
        generatedDictionary [key] = paddedVal
    return generatedDictionary
####################
### Generating different sequence and then calculating distance
####################
def generateSequenceRandom (tfaDictionary):
    generatedDictionary ={}
    for key, val in tfaDictionary.items():
        #inserting random space
        randLoc = randrange(96)
        generatedDictionary [key] = val[:randLoc] + BLANK_CHAR + val[randLoc:]
    return generatedDictionary
 #################
 #Function to generate different sequences and then calculating distances
 ##################
 
 
def calculateDistance(inputGeneratedDictionary):
    generatedDictionary = inputGeneratedDictionary.copy()
    keyList = list(generatedDictionary.keys())
    distance = 0
    for key in keyList:
        val1 = generatedDictionary.pop(key)
        for val2 in generatedDictionary.values():
            distance = distance + calculateDistanceForPairUsingMatrix (val1, val2)
    return distance    

def calculateDistanceForPair(str1, str2):
    retVal=0
    for index, val in enumerate(str1):
        if index >= len(str2):
             retVal = retVal + 0 #Ignore of second string is smaller than first.
        elif val == BLANK_CHAR and str2[index]:
            retVal = retVal + 0
        elif val == str2[index]:
            retVal = retVal + 1
        else:
            retVal = retVal -1 
    return retVal
    
def calculateDistanceForPairUsingMatrix(str1, str2):
    retVal=0
    gapOpened1 = False
    gapOpened2 = False
    for index, val in enumerate(str1):
        if index >= len(str2):
            retVal = retVal + 0 #Ignore of second string is smaller than first.
        elif val == BLANK_CHAR and str2[index] == BLANK_CHAR:
            retVal = retVal + 0 #Ignore if both the chars are blank.
        elif val == BLANK_CHAR:
            if  gapOpened1:
                retVal = retVal + GAP_EXTENTION_PENALTY
            else:
                retVal = retVal + GAP_OPENING_PENALTY
                gapOpened1 = True
                
        elif str2[index] == BLANK_CHAR:
            if  gapOpened2:
                retVal = retVal + GAP_EXTENTION_PENALTY
            else:
                retVal = retVal + GAP_OPENING_PENALTY
                gapOpened2 = True
        else:
            gapOpened1 = False
            gapOpened2 = False
            retVal = retVal + int(distanceMatrixDict[val][str2[index]])
    return retVal

def columnScore (referenceMsa, inputMsa):
    columnScoreArray = []
    i = 0
    msaLength = len(list(referenceMsa.values())[0])
    inputMsaLength = len(list(inputMsa.values())[0]) 
    for i in range (msaLength):
        columnScoreArray.append (0)
    minRange = 0
    maxRange = 0
    for i in range (msaLength):
        if minRange < i-5:
            minRange = i-5
        maxRange = i +10
        if maxRange > inputMsaLength:
            maxRange = inputMsaLength
        for newI in range (minRange, maxRange):
            matched = True
            for key, str1 in referenceMsa.items():
                try:
                    val1 =str1[i]
                except:
                    val1 = ""
                    print ("ERROR the reference MSA has less number of columns for key=" + key + " MSALength= "+ str(msaLength)+ " key Index = " + str(i))
                    break
                try:   
                    val2 = str(inputMsa.get(key))[newI]
                except:
                    val2 = ""
                    print ("ERROR the input MSA has less number of columns for key=" + key + " MSALength= "+ str(msaLength)+ " key Index = " + str(newI))
                    break
                if  val1 != val2:
                    matched = False
                    break
            if matched:
                columnScoreArray[i] = 1
                minRange = newI
                break
    print ( " Number of columns matched = " + str (sum(columnScoreArray)) + " total columns = " + str(len(columnScoreArray)))
    return (sum(columnScoreArray)/len(columnScoreArray))
            
    
def printDictionary (dict):
    for key, val in dict.items():
        print (key + ": " + val)

def saveDictionary (dict, outputPath, fileName):
    oFile = open(outputPath + "/" + fileName, "w")
    for key, val in dict.items():
       oFile.write (key + "\n")
       oFile.write (val+ "\n")
    oFile.close()
def updateResult (resultFileName, fileName, initialDistance,refDistance, inputColScore, leastDistance,  outputColScore):
    rFile = open(resultFileName, "a")
    rFile.write (fileName + "," + str(initialDistance) + "," + str (refDistance) + "," + str (inputColScore) + "," + str (leastDistance) + "," + str (outputColScore) + "\n")
    rFile.close()
    
#####################################################
# Start of the program                          #####
# Change the file name to process another file. #####
# File should be in same directory where the    #####
# program is running.                           #####
#####################################################
loadDistanceMatrix()


######################
#Below code to test the compare function.
#print (calculateDistanceForPair ("abc.xz", "abc.y."))
#exit()
########################
#Code to test the calculate distance function
#dict1 = {"k1":"abc", "k2":"abc", "k3": "abc", "k4": "xbz"}
#print (calculateDistance(dict1))
#exit()
#Code to test the matrix based distance calculation
#print ("Testing the distance function with test pairs")
print ("AAAAA, AAAAA " + str(calculateDistanceForPairUsingMatrix ("AAAAA", "AAAAA")) )  
print ("AA-AAA, AA-AAA " + str(calculateDistanceForPairUsingMatrix ("AA-AAA", "AA-AAA")) )  
#print ("AAAAA , ARNAA "  +str((calculateDistanceForPairUsingMatrix ("AAAAA", "ARNAA")) ))
#print ("A--AA , ARNAA "  +str((calculateDistanceForPairUsingMatrix ("A--AA", "ARNAA") )))
#print ("AAAAA, A--AA "  +str((calculateDistanceForPairUsingMatrix ("AAAAA", "A--AA")) ))
#print ("AAAAARRRRR, A--AA--RRR "  +str((calculateDistanceForPairUsingMatrix ("AAAAARRRRR", "A--AA--RRR")) ))
#exit()



#f = open("BB11001.tfa", "r")
outDir = OUTPUT_FOLDER + "-" + str(datetime.datetime.now()).replace (":", "-")
print ("creating output directory " + outDir )
p = Path(outDir)
p.mkdir(exist_ok=True)
rFile = open(outDir + "/" + RESULT_FILE_NAME, "w")
rFile.write ("File Name,Input File Distance,Reference File Distance,Input Column Score,Output File Distance, Output Column Score\n" )
rFile.close()
fileList = os.listdir (INPUT_FOLDER)
for fileName in fileList:
    try:
        print ("==================================Processing File " + fileName + " =================================")
        inputFile = open(INPUT_FOLDER + fileName, "r")

        tfaDictionary = getTfaDictionary(inputFile)
        inputFile.close()

        refFile = open(REFERENCE_FOLDER + fileName, "r")

        refTfaDictionary = getTfaDictionary(refFile)
        refFile.close()
        inputColScore = columnScore (refTfaDictionary, tfaDictionary)
        print ("Column score for input file: " + str (inputColScore))

        #print ("-------------Input from the file -----------")
        #printDictionary(tfaDictionary)
        initialDistance = calculateDistance(tfaDictionary)
        refDistance = calculateDistance(refTfaDictionary)
        print ("Sum of pair score for reference file: " + str (refDistance))
        print ("Sum of pair score for input file: " + str (initialDistance))

        #msfDictionary = generateMsf(tfaDictionary)
        leastDistance = initialDistance
        bestMsfDictionary = tfaDictionary.copy()
        #exit()
        for x in range(LOOP_COUNT):
            if x % 500 == 0:
                sys.stdout.write(".") 
                sys.stdout.flush()
            #print ("------------Generated MSA----------")
            msfDictionary = generateSequenceRandom (tfaDictionary)
            #printDictionary (msfDictionary)
            #print ("------------Distance/Match for generated MSA -------")
            newDist = calculateDistance(msfDictionary)
            #print (newDist)
            if newDist > leastDistance:
                leastDistance = newDist
                bestMsfDictionary = msfDictionary.copy()
                #print ("This is the best match")
                #break
        print ("")
        print ("loop ended of number of iterations=" + str(x + 1) + " maxMatch =" + str(leastDistance) + " intitialMatch = " + str(initialDistance))
        print ("Best alginment")
        #printDictionary (bestMsfDictionary)
        saveDictionary (bestMsfDictionary, outDir, fileName)
        outputColScore = columnScore (refTfaDictionary, bestMsfDictionary)
        updateResult (outDir + "/" + RESULT_FILE_NAME, fileName, initialDistance,refDistance, inputColScore, leastDistance,  outputColScore)
        print ("Column score for output file: " + str (outputColScore))
        print ("Sum of pair score for output file: " + str (leastDistance))
    except Exception as ex:
        print ("!!!!!! Error processing the file " + fileName )
        traceback.print_exception(*sys.exc_info())
       



