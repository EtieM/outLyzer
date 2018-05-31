#! /usr/bin/env python
#! -*- coding: utf-8 -*-


import sys
import os
import re
import subprocess
import numpy as np
from scipy.stats import t
import math
import argparse
import time
from multiprocessing import Process
from argparse import RawTextHelpFormatter
from shutil import copyfile
from collections import defaultdict


### MultiProcessing

def launchCommandLine(commandLine):
    os.system(commandLine)
    
def ifFolderNotExistCreate(path):
    if (os.path.exists(path)) == False:
        os.makedirs(path)

def manageCommandListToLaunch(commandList,maxJob):
    jobs = []
    for commandLine in commandList:
        if len(jobs) == maxJob:
            while len(jobs) >= maxJob:
                for job in jobs:
                    if not job.is_alive():
                        jobs.remove(job)
                time.sleep(3)
        process = Process(target=launchCommandLine, args=(commandLine, ))
        process.start()
        jobs.append(process)
    while len(jobs) > 0 :
        for job in jobs:
            if not job.is_alive():
                jobs.remove(job)
        time.sleep(3)

### outLyzer

def dataListCreation(datafile,sep):
    """This function reads csv file containing variant identification
    and create a list of lists to analyse it.
    Argument 1: datafile = csv file to analyse"""
    liste = datafile.readlines()
    liste2 =[]
    for i in liste:
        j =i.split(sep)
        liste2.append(j)
    return liste2

def createBedList(bedFilePath):
    bedFile = open(bedFilePath,'r')
    bedList = dataListCreation(bedFile, '\t')
    bedFile.close()
    for line in bedList:
        line[2] = line[2].replace('\n','')
    return bedList

def cutBedFile(bedFilePath,cutNb,outputPath):
    ifFolderNotExistCreate(outputPath)
    bedFile = open(bedFilePath,'r')
    bedList = dataListCreation(bedFile, '\t')
    bedFile.close()
    lineNb = int(math.ceil(len(bedList)/float(cutNb)))
    BedIntervals = range(0,len(bedList),lineNb)
    fileNb = 1
    for interval in BedIntervals[1:]:
        fileToWrite = open('%sbedFile_%i.bed'%(outputPath,fileNb),'w')
        for line in bedList[interval-lineNb:interval]:
            lineToWrite = '\t'.join(line)
            fileToWrite.write(lineToWrite)
        fileToWrite.close()
        fileNb += 1
    fileToWrite = open('%sbedFile_%i.bed'%(outputPath,fileNb),'w')
    for line in bedList[BedIntervals[-1]:len(bedList)]:
        lineToWrite = '\t'.join(line)
        fileToWrite.write(lineToWrite)
    fileToWrite.close()
    bedList = os.listdir(outputPath)
    finalBedList = []
    for bedFile in bedList:
        finalBedList.append('%s%s'%(outputPath,bedFile))
    return finalBedList


def createOutLyzerCommandList(output,pythonPath,samtoolsPath,progPath,cut,bedList,bamFile,refFile,studentValue,balance,Qscore,SDQ,WinSize,AS,HSM,multiplicativeFactor,verbose,FRcor,WSmin):
    outLyzerCommandList = []
    for bedFile in bedList:
        commandLine = '%s %s subprocess -samtools %s -bam %s -cut %i -output %s -ref %s -bed %s -t %f -bal %f -Q %i -SDQ %i -WS %i -x %f -verbose %i -WSmin %i'%(pythonPath,progPath,samtoolsPath,bamFile,cut,output,refFile,bedFile,studentValue,balance,Qscore,SDQ,WinSize,multiplicativeFactor,verbose,WSmin)
        if AS:
            commandLine = commandLine+' -AS'
        if HSM:
            commandLine = commandLine+' -HSM %s'%HSM
        if FRcor:
            commandLine = commandLine+' -FRcor'
        outLyzerCommandList.append(commandLine)
    return outLyzerCommandList


def launchPileupCommandForBedFile(bamFile,samtoolsPath,referenceFile,bedFile,output):
    samtoolsOutput = subprocess.check_output('%s mpileup -d 100000 -Q 0 -A -R -B -f %s -x -l %s %s 2>> %ssamtoolsError.txt'%(samtoolsPath,referenceFile,bedFile,bamFile,output),shell=True)
    pileupList = samtoolsOutput.split('\n')
    pileupList.pop()
    return pileupList

def launchPileupCommandForPosition(bamFile,samtoolsPath,referenceFile,genomicPosition,windowSize):
    name = bamFile.split('.')[0].split('_')[0]
    name = name.split('/')[-1]
    chromosome = genomicPosition.split(':')[0]
    positiontoTest = int(genomicPosition.split(':')[-1])
    samtoolsPosition = subprocess.check_output('%s mpileup -d 10000 -Q 0 -A -R -B -f %s -x -r %s:%i-%i %s'%(samtoolsPath,referenceFile,chromosome,positiontoTest,positiontoTest,bamFile),shell=True)
    if samtoolsPosition == '':
        print 'No reads mapped on position'
    else:
        genomicInterval = '%s:%i-%i'%(chromosome,round(positiontoTest-(windowSize/2)),round(positiontoTest+(windowSize/2)))
        samtoolsOutput = subprocess.check_output('%s mpileup -d 10000 -Q 0 -A -R -B -f %s -x -r %s %s'%(samtoolsPath,referenceFile,genomicInterval,bamFile),shell=True)
        pileupList = samtoolsOutput.split('\n')
        pileupList.pop()
        return pileupList

def calcThompsonTau(n,studentValue):
    tStu = t.ppf(1-(studentValue/2),n-2)
    tau = (tStu*(n-1))/(math.sqrt(n)*math.sqrt(n-2+tStu**2))
    return tau

def thompsonTest(listToTest,studentValue):
    condition = 'continue'
    mean = np.mean(listToTest)
    stdDev = np.std(listToTest)
    deltaMax = abs(listToTest[-1]-mean)
    tau = calcThompsonTau(len(listToTest),studentValue)
    testValue = tau*stdDev
    if deltaMax > testValue:
        listToTest.pop()
    else:
        condition = 'stop'
    return listToTest,condition
    
    
def reject_outliersReads(altReadList,studentValue):
    filteredAltReadList = []
    for number in altReadList:
        num = int(number)
        if num != 0:
            filteredAltReadList.append(num)
    filteredAltReadList.sort()
    outliersReject = 'continue'
    if not not filteredAltReadList :
        while outliersReject == 'continue':
            newAltReadList,condition = thompsonTest(filteredAltReadList,studentValue)
            if condition=='stop':
                break
    else:
        newAltReadList = []
    return newAltReadList
        

def calcConfidenceIntervalFromReads(altReadList,studentValue,DP):
    if altReadList.count('0')> 0.90*len(altReadList):
        maxOutliers = math.ceil(DP*0.025)
    else:
        altReadsOutliers = reject_outliersReads(altReadList,studentValue)
        if not len(altReadsOutliers) == 0:
            maxOutliers = max(altReadsOutliers)
        else:
            maxOutliers = 0
    return maxOutliers

def forwardReverseOrigin(Indel):
    if str.isupper(Indel):
        return 'Forward'
    else:
        return 'Reverse'

def countMaxIndel(indelDic):
    maxIndel = ''
    maxCount = 0
    for element in indelDic:
        totalReads = int(indelDic[element]['Forward'])+int(indelDic[element]['Reverse'])
        if totalReads > maxCount:
            maxCount = totalReads
            maxIndel = element
    return maxIndel

def countDelInPileup(baseCount):
    delDic = {}
    while bool(re.search('[-]',baseCount)):
        indelSearch = re.search('[-]([0-9]+)+([ATCGNatcgn]+)',baseCount)
        indelReal = indelSearch.group(0)[0:int(indelSearch.group(1))+len(indelSearch.group(1))+1]
        indelContent = indelSearch.group(2)[0:int(indelSearch.group(1))]
        orientation = forwardReverseOrigin(indelContent)
        if not delDic.has_key(indelReal.upper()):
            delDic[indelReal.upper()] = {'Forward':0,'Reverse':0}
            delDic[indelReal.upper()][orientation] = baseCount.count(indelReal)
        else:
            delDic[indelReal.upper()][orientation] = baseCount.count(indelReal)
        baseCount = baseCount.replace(indelReal,'')
    if bool(delDic):
        mainDel = countMaxIndel(delDic)
        return mainDel,delDic[mainDel]
    else:
        return 'None',delDic
            
def countInsInPileup(baseCount):
    delDic = {}
    while bool(re.search('[+]',baseCount)):
        insSearch = re.search('[+]([0-9]+)+([ATCGNatcgn]+)',baseCount)
        insReal = insSearch.group(0)[0:int(insSearch.group(1))+len(insSearch.group(1))+1]
        insContent = insSearch.group(2)[0:int(insSearch.group(1))]
        orientation = forwardReverseOrigin(insContent)
        if not delDic.has_key(insReal.upper()):
            delDic[insReal.upper()] = {'Forward':0,'Reverse':0}
            delDic[insReal.upper()][orientation] = baseCount.count(insReal)
        else:
            delDic[insReal.upper()][orientation] = baseCount.count(insReal)
        baseCount = baseCount.replace(insReal,'')
    if bool(delDic):
        mainDel = countMaxIndel(delDic)
        return mainDel,delDic[mainDel]
    else:
        return 'None',delDic
            

def countPointMutations(baseCount,pileupInfosDic):
    baseCountNoIndel = eraseIndelfromPileup(baseCount)
    baseCountDic={'For':baseCountNoIndel.count('.'),'Rev':baseCountNoIndel.count(','),'a':baseCountNoIndel.count('a'),'A':baseCountNoIndel.count('A'),'c':baseCountNoIndel.count('c'),'C':baseCountNoIndel.count('C'),'g':baseCountNoIndel.count('g'),'G':baseCountNoIndel.count('G'),'t':baseCountNoIndel.count('t'),'T':baseCountNoIndel.count('T')}
    qualDic={'A':[],'C':[],'G':[],'T':[]}
    baseList = ['a','A','c','C','g','G','t','T']
    for position,base in enumerate(baseCountNoIndel):
                    if base in baseList:
                        phredQual = ord(pileupInfosDic['baseQual'][position])
                        if base == 'a' or base =='A':
                            qualDic['A'].append(phredQual)
                        elif base == 'c' or base =='C':
                            qualDic['C'].append(phredQual)
                        elif base == 'g' or base =='G':
                            qualDic['G'].append(phredQual)
                        elif base == 't' or base =='T':
                            qualDic['T'].append(phredQual)
    return baseCountDic,qualDic

def eraseIndelfromPileup(baseCount):
    while bool(re.search('[+-]',baseCount)):
        indelSearch = re.search('[+-]([0-9]+)+[ATCGNatcgn]+',baseCount)
        indelReal = indelSearch.group(0)[0:int(indelSearch.group(1))+len(indelSearch.group(1))+1]
        baseCount = baseCount.replace(indelReal,'')
    return baseCount
 
def mainPileupLineInfosInDic(pileupLine):
    pileupInfosDic = {}
    pileupInfosDic['refBase'] = pileupLine[2]
    pileupInfosDic['depth'] = pileupLine[3]
    pileupInfosDic['baseCount'] = re.sub(r"(\^.)",r"",pileupLine[4].replace('$',''))
    pileupInfosDic['baseQual'] = pileupLine[5]
    return pileupInfosDic

def defineMainPointMutation(pileupLine, pileupInfosDic):
    baseCountDic,qualDic = countPointMutations(pileupLine, pileupInfosDic)
    nucleotideList=['A','C','G','T']
    finalAltBase = ''
    finalAltBaseFreq = 0
    finalBaseQual = 0
    finalQualStd = 0
    altReads = 0
    for base in nucleotideList:
        if not base == pileupInfosDic['refBase']:
            totalBaseCount = baseCountDic[base]+baseCountDic[base.lower()]
            if totalBaseCount == 0:
                phredQual = 'NA'
                phredQualStd = 'NA'
            else:
                phredQual = round(float(sum(qualDic[base]))/(baseCountDic[base.lower()]+baseCountDic[base]),2)-33
                phredQualStd = round(np.std(qualDic[base]),2)
            altFreqBase = round(((baseCountDic[base]+baseCountDic[base.lower()])/float(pileupInfosDic['depth']))*100,3)
            if altFreqBase > finalAltBaseFreq:
                finalAltBase = base
                finalAltBaseFreq = altFreqBase
                finalBaseQual = phredQual
                finalQualStd = phredQualStd
                altReads = baseCountDic[base]+baseCountDic[base.lower()]
    return finalAltBase,finalAltBaseFreq,finalBaseQual,finalQualStd,altReads

def defineLeadingMutationType(delReads,insReads,baseCountDic,finalAltBase):
    leadingMutation = 'WT'
    leadingMutationReads = 0
    if bool(delReads):
        totalDelReads = int(delReads['Forward'])+int(delReads['Reverse'])
        if totalDelReads > leadingMutationReads:
            leadingMutationReads = totalDelReads
            leadingMutation = 'Deletion'
    else:
        totalDelReads = 0
    if bool(insReads):
        totalInsReads = int(insReads['Forward'])+int(insReads['Reverse'])
        if totalInsReads > leadingMutationReads:
            leadingMutationReads = totalInsReads
            leadingMutation = 'Insertion'
    else:
        totalInsReads = 0
    if finalAltBase != '':
        pointAltReads = int(baseCountDic[finalAltBase])+int(baseCountDic[finalAltBase.lower()])
        if pointAltReads > leadingMutationReads:
            leadingMutationReads = pointAltReads
            leadingMutation = 'PointMutation'
    else:
        pointAltReads = 0
    return leadingMutation


def countAltReads(baseCount):
    baseCountNoIndel = eraseIndelfromPileup(baseCount)
    indelCount = baseCount.count('-')+baseCount.count('+')
    altReads = baseCountNoIndel.count('n')+baseCountNoIndel.count('N')+baseCountNoIndel.count('a')+baseCountNoIndel.count('A')+baseCountNoIndel.count('c')+baseCountNoIndel.count('C')+baseCountNoIndel.count('g')+baseCountNoIndel.count('G')+baseCountNoIndel.count('t')+baseCountNoIndel.count('T')
    totalAltReads = indelCount + altReads
    return totalAltReads

def extractpileupLineInfos(pileupLine):
    mutationStats = {}
    mutationStats['Chromosome']= pileupLine[0]
    mutationStats['GenomicPosition'] = pileupLine[1]
    mutationStats['refBase'] = pileupLine[2]
    mutationStats['DP'] = pileupLine[3] 
    baseCount = re.sub(r"(\^.)",r"",pileupLine[4].replace('$',''))
    pileupInfosDic = mainPileupLineInfosInDic(pileupLine)
    delName,delReads = countDelInPileup(baseCount)
    insName,insReads = countInsInPileup(baseCount)
    baseCountDic, qualDic = countPointMutations(baseCount, pileupInfosDic)
    finalAltBase,finalAltBaseFreq,finalBaseQual,finalQualStd,altReads = defineMainPointMutation(baseCount, pileupInfosDic)
    leadingMutation = defineLeadingMutationType(delReads, insReads, baseCountDic, finalAltBase)
    totalAlt = countAltReads(baseCount)
    mutationStats['totalAlt'] = int(totalAlt)
    mutationStats['WT_Bal'] = '%i/%i'%(int(baseCountDic['For']),int(baseCountDic['Rev']))
    if leadingMutation == 'Insertion':
        mutationStats['altReadsMutation'] = insReads['Forward']+insReads['Reverse']
        mutationStats['AF'] = round((float(insReads['Forward'])+float(insReads['Reverse']))*100 / float(pileupInfosDic['depth']),3)
        mutationStats['alt'] = insName
        mutationStats['Qual'] = '.'
        mutationStats['Qual_StdDev'] = '.'
        mutationStats['balance'] = '%s/%s'%(insReads['Forward'],insReads['Reverse'])
    elif leadingMutation == 'Deletion':
        mutationStats['altReadsMutation'] =  delReads['Forward']+delReads['Reverse']
        mutationStats['AF'] = round((float(delReads['Forward'])+float(delReads['Reverse']))*100 / float(pileupInfosDic['depth']),3)
        mutationStats['alt'] = delName
        mutationStats['Qual'] = '.'
        mutationStats['Qual_StdDev'] = '.'
        mutationStats['balance'] = '%s/%s'%(delReads['Forward'],delReads['Reverse'])
    elif leadingMutation =='PointMutation':
        mutationStats['altReadsMutation'] = altReads
        mutationStats['AF'] = finalAltBaseFreq
        mutationStats['alt'] = finalAltBase
        mutationStats['Qual'] = finalBaseQual
        mutationStats['Qual_StdDev'] = finalQualStd
        mutationStats['balance'] = '%s/%s'%(baseCountDic[finalAltBase],baseCountDic[finalAltBase.lower()])
    else:
        mutationStats['altReadsMutation'] = '0'
        mutationStats['AF'] = '.'
        mutationStats['alt'] = 'WT'
        mutationStats['Qual'] = '.'
        mutationStats['Qual_StdDev'] = '.'
        mutationStats['balance'] = '.'
    return mutationStats

def writeVcfHeader(commandLine,vcfFile):
    openVcfFile = open(vcfFile,'w')
    patientID = vcfFile.split('/')[-1].split('.')[0].split('_')[0]
    header="""##fileformat=VCFv4.2
##fileDate=%s
##source=%s 
##reference=/data/genome/hg19/hg19_order.fa
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=SDQ,Number=1,Type=Float,Description="Standard Deviation of average Phred Score">
##INFO=<ID=OUTLIER,Number=1,Type=Integer,Description="Background Noise Threshold (Number of Reads)">
##FILTER=<ID=noise_background,Description="Mutation's Rate is Higher than Surrounding Background Noise">
##FORMAT=<ID=BAL,Number=2,Type=String,Description="Count of alternative Forward / Reverse Reads">
##FORMAT=<ID=WTbal,Number=1,Type=String,Description="Forward / Reverse WT reads at the position ">
##FORMAT=<ID=MSS,Number=1,Type=String,Description="Proximity analysis of Stretch and repetition motifs around mutation ">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n"""%(time.strftime("%Y%m%d"),commandLine,patientID)
    openVcfFile.write(header)
    openVcfFile.close()

def defineWindowToScan(altReadList,positionInExon,windowSize):
    if len(altReadList)< windowSize:
        newAltReadList = altReadList
    else:
        if positionInExon-round(windowSize/2) < 0:
            newAltReadList = altReadList[0:windowSize]
        elif positionInExon+round(windowSize/2)>len(altReadList):
            newAltReadList = altReadList[-windowSize:]
        else:
            newAltReadList = altReadList[int(positionInExon-round(windowSize/2)):int(positionInExon+round(windowSize/2))]
    return newAltReadList

def getsubs(loc, s):
    substr = s[loc:]
    i = -1
    while(substr):
        yield substr
        substr = s[loc:i]
        i -= 1

def testMotifRec(seq):
    check=False
    for i in xrange(2,int(len(seq)/2),1):
        first = seq[0:i]
        second = seq[i:2*i]
        if first == second:
            if len(first)!=first.count(first[0]):
                check=True
                break
    return check

def stretchValid(seq):
    if len(seq)==seq.count(seq[0]):
        return True
    else:
        return False

def motifDistanceFromMutation(motif,seq, mutPos):
    subStart = seq.find(motif)
    dist1 = subStart-mutPos
    dist2 = (subStart+len(motif))-mutPos
    distList = [abs(dist1),abs(dist2)]
    finalDist = 0
    if dist1<=0 and dist2>0:
        finalDist = 0
    else:
        finalDist = min(distList)
    return finalDist

def findMotif(seq,minStretch,mutPos):
    occ = defaultdict(int)
    for i in range(len(seq)):
        for sub in getsubs(i,seq):
            occ[sub] += 1
    maxStretch = ''
    maxMotif = ''
    for motif in occ:
        if occ[motif]>=2 and len(motif)>6 and len(motif)>len(maxMotif) and not stretchValid(motif):
            if testMotifRec(motif):
                maxMotif = motif
        if len(motif)>=minStretch and stretchValid(motif) and len(motif)>len(maxStretch):
            if len(motif)!=0:
                maxStretch=motif
    motDic = {}
    if len(maxStretch)>0:
        motDic['stretch']=[maxStretch,motifDistanceFromMutation(maxStretch, seq, mutPos)]
    if len(maxMotif)>0:
        motDic['motif']=[maxMotif,motifDistanceFromMutation(maxMotif, seq, mutPos)]
    return motDic

def defineMotifStatus(seq,mutPos):
    motifStatus = ''
    motifDic = findMotif(seq, 5, mutPos)
    if motifDic.has_key('stretch'):
        motifStatus += 'S%id%i'%(len(motifDic['stretch'][0]),motifDic['stretch'][1])
    if motifDic.has_key('motif'):
        if len(motifStatus)!=0:
            motifStatus += '|'
        motifStatus += 'M%id%i'%(len(motifDic['motif'][0]),motifDic['motif'][1])
    if len(motifStatus) == 0:
        motifStatus = 'None'
    return motifStatus
    
        

def likelihoodMutationState(mutationStats,balance,quality,SDquality,FRcor):
    warning = 0
    if mutationStats[4] != 'WT':
        if FRcor:
            WT_F = float(mutationStats[12].split('/')[0])
            WT_R = float(mutationStats[12].split('/')[1])
            Alt_F = float(mutationStats[11].split('/')[0])
            Alt_R = float(mutationStats[11].split('/')[1])
            if len(mutationStats[3])==1 and len(mutationStats[4])==1:
                diffF = round((WT_F+Alt_F) / (float(mutationStats[5])),3)-0.5
            else:
                diffF = round((WT_F) / (float(mutationStats[5])),3)-0.5
            if diffF > 0:
                balanceProp = round(((Alt_F-(Alt_F*diffF))/(Alt_F+Alt_R)),3)
            elif diffF < 0:
                balanceProp = round(((Alt_R+(Alt_R*diffF))/(Alt_F+Alt_R)),3)
            else:
                balanceProp = round((Alt_F/(Alt_F+Alt_R)),3)
        else:
            balanceProp = float(mutationStats[11].split('/')[0]) /( float(mutationStats[11].split('/')[0])+float(mutationStats[11].split('/')[1]))
    else:
        balanceProp = 'NA'
    if mutationStats[4].startswith('+') or mutationStats[4].startswith('-'):
        if balanceProp >= (balance) and balanceProp <= (1-balance):
            likelihoodStatus=1
        else:
            likelihoodStatus=0
    elif mutationStats[4] != 'WT':
        if float(mutationStats[9]) >= quality and balanceProp >= balance and balanceProp <= 1-balance and float(mutationStats[10])<SDquality:
            likelihoodStatus=1
        else:
            likelihoodStatus=0
    else:
        likelihoodStatus=0
    return likelihoodStatus

def calcExonBackGroundNoise(altReadList,depthList,studentValue):
    meanDepth = sum(depthList)/len(depthList)
    maxOutlier = calcConfidenceIntervalFromReads(altReadList, studentValue,meanDepth)
    sensitivityLimit = float(maxOutlier)/meanDepth
    return sensitivityLimit

def defineMutationPosInSeq(altReadList,positionInExon,windowSize):
    if len(altReadList)< windowSize:
        mutationPosition = positionInExon
    else:
        if positionInExon-round(windowSize/2) < 0:
            mutationPosition = positionInExon
        elif positionInExon+round(windowSize/2)>len(altReadList):
            mutationPosition = round(windowSize/2)+(positionInExon+round(windowSize/2)-len(altReadList))
        else:
            mutationPosition = round(windowSize/2)
    return mutationPosition

def detectMutations(exonMatrix,studentValue,balance,quality,SDquality,windowSize,multiplicativeFactor,fileToWrite,FRcor):
    altReadList = exonMatrix[:,7].tolist()
    refSeqList = exonMatrix[:,3].tolist()
    depthList = [int(i) for i in exonMatrix[:,5].tolist()]
    exonSensitivityLimit = calcExonBackGroundNoise(altReadList, depthList, studentValue)
    for rowNum in range(0,len(exonMatrix)):
        windowList = defineWindowToScan(altReadList, int(rowNum), windowSize)
        maxOutliers = calcConfidenceIntervalFromReads(windowList, studentValue,int(exonMatrix[rowNum,5]))
        if int(exonMatrix[rowNum,6])>(maxOutliers*multiplicativeFactor):
            isReal = likelihoodMutationState(exonMatrix[rowNum,:],balance,quality,SDquality,FRcor)
            if isReal:
                refSeqToScan = ''.join(defineWindowToScan(refSeqList, int(rowNum), 60))
                mutPos = defineMutationPosInSeq(refSeqList, int(rowNum), 60)
                motifStatus = defineMotifStatus(refSeqToScan,mutPos)
                if '+' in exonMatrix[rowNum,4] or '-' in exonMatrix[rowNum,4]:
                    ref,alt = transformIndelAnnotation(exonMatrix[rowNum,3], exonMatrix[rowNum,4])
                    lineToWrite = '%s\t%s\t.\t%s\t%s\t%s\tPASS\tDP=%s;AF=%s;SDQ=%s;OUTLIER=%s\tBAL:WTbal:MSS\t%s:%s:%s\n'%(exonMatrix[rowNum,1],exonMatrix[rowNum,2],ref,alt,exonMatrix[rowNum,9],exonMatrix[rowNum,5],exonMatrix[rowNum,8],exonMatrix[rowNum,10],maxOutliers,exonMatrix[rowNum,11],exonMatrix[rowNum,12],motifStatus)
                    fileToWrite.write(lineToWrite)
                else: 
                    lineToWrite = '%s\t%s\t.\t%s\t%s\t%s\tPASS\tDP=%s;AF=%s;SDQ=%s;OUTLIER=%s\tBAL:WTbal:MSS\t%s:%s:%s\n'%(exonMatrix[rowNum,1],exonMatrix[rowNum,2],exonMatrix[rowNum,3].upper(),exonMatrix[rowNum,4].upper(),exonMatrix[rowNum,9],exonMatrix[rowNum,5],exonMatrix[rowNum,8],exonMatrix[rowNum,10],maxOutliers,exonMatrix[rowNum,11],exonMatrix[rowNum,12],motifStatus)
                    fileToWrite.write(lineToWrite)
    return maxOutliers,exonSensitivityLimit

def transformIndelAnnotation(ref,alt):
    regex = re.compile('[\+\-][0-9]+')
    if '+' in alt:
        newAlt = regex.sub(ref,alt).upper()
        newRef = ref.upper()
    else:
        newRef = regex.sub(ref,alt).upper()
        newAlt = ref.upper()
    return newRef,newAlt

def hotSpotStartPos(hotSpotBed):
    hotSpotBedFile = open(hotSpotBed,'r')
    hotSpotList = dataListCreation(hotSpotBedFile, '\t')
    hotSpotBedFile.close()
    hotSpotDic = {}
    for line in hotSpotList:
        chrom = line[0]
        start = line[1]
        infos = line[2].replace('\n','').replace('\r','')
        hotSpotDic['%s_%s'%(chrom,start)] = infos
    return hotSpotDic


def readPileup(pileup,bedNum,studentValue,balance,quality,SDquality,windowSize,vcfPath,multiplicativeFactor,HSMfile,FRcor,WSmin):
    previousPos = ''
    firstLine = 1
    countPosInExon = 0
    unCoveredRegions = []
    lastLine = '%s_%s'%(pileup[-1].split('\t')[0],pileup[-1].split('\t')[1])
    fileToWrite = open('%s/tempCalling_%i.vcf'%(vcfPath,int(bedNum)),'a')
    if HSMfile != None:
        hotSpotDic = hotSpotStartPos(HSMfile)
        metricsFile = open('%s/tempHSM_%s.txt'%(vcfPath,bedNum),'w')
    hotSpotMark = {}
    exonSensitivity = {}
    for pileupLine in pileup:
        pileupLine = pileupLine.split('\t')
        Chromosome = pileupLine[0]
        genomicPos = int(pileupLine[1])
        mutationStats = extractpileupLineInfos(pileupLine)
        rowToAdd = [countPosInExon,mutationStats['Chromosome'],mutationStats['GenomicPosition'],mutationStats['refBase'],mutationStats['alt'],mutationStats['DP'],mutationStats['altReadsMutation'],int(mutationStats['totalAlt']),mutationStats['AF'],mutationStats['Qual'],mutationStats['Qual_StdDev'],mutationStats['balance'],mutationStats['WT_Bal']]
        if HSMfile != None:
            if hotSpotDic.has_key('%s_%i'%(Chromosome,genomicPos)):
                hotSpotMark['%s_%i'%(Chromosome,genomicPos)] = int(mutationStats['DP'])
        if firstLine == 1:
            exonMatrix = [rowToAdd]
            countPosInExon += 1
            firstLine = 0
            previousPos = genomicPos
        else:
            if previousPos+1 == genomicPos:
                exonMatrix.append(rowToAdd)
                previousPos = genomicPos
                countPosInExon += 1
                if '%s_%i'%(Chromosome,genomicPos) == lastLine:
                    exonMatrix = np.asarray(exonMatrix)
                    if len(exonMatrix)>WSmin:
                        maxOutlier,exonSensitivityLimit = detectMutations(exonMatrix,studentValue,balance,quality,SDquality,windowSize,multiplicativeFactor,fileToWrite,FRcor)
                        exonSensitivity['%s_%s_%s'%(exonMatrix[0,1],exonMatrix[0,2],exonMatrix[-2,2])]=exonSensitivityLimit
                        if HSMfile != None:
                            if not not hotSpotMark:
                                for position in hotSpotMark:
                                    metricsFile.write('%s\t%f\n'%(hotSpotDic[position],round(float(maxOutlier)/int(hotSpotMark[position])*100*multiplicativeFactor,2)))
                    else:
                        uncovRegion = '%s_%s_%s'%(exonMatrix[0,1],exonMatrix[0,2],exonMatrix[-1,2])
                        unCoveredRegions.append(uncovRegion)
            else:
                exonMatrix = np.asarray(exonMatrix)
                if len(exonMatrix)>WSmin:
                    maxOutlier,exonSensitivityLimit = detectMutations(exonMatrix,studentValue,balance,quality,SDquality,windowSize,multiplicativeFactor,fileToWrite,FRcor)
                    exonSensitivity['%s_%s_%s'%(exonMatrix[0,1],exonMatrix[0,2],exonMatrix[-2,2])]=exonSensitivityLimit
                else:
                    uncovRegion = '%s_%s_%s'%(exonMatrix[0,1],exonMatrix[0,2],exonMatrix[-1,2])
                    unCoveredRegions.append(uncovRegion)
                countPosInExon = 0
                rowToAdd = [countPosInExon,mutationStats['Chromosome'],int(mutationStats['GenomicPosition']),mutationStats['refBase'],mutationStats['alt'],mutationStats['DP'],mutationStats['altReadsMutation'],int(mutationStats['totalAlt']),mutationStats['AF'],mutationStats['Qual'],mutationStats['Qual_StdDev'],mutationStats['balance'],mutationStats['WT_Bal']]
                exonMatrix = [rowToAdd]
                previousPos = genomicPos
                countPosInExon += 1
                if HSMfile != None:
                    if not not hotSpotMark:
                        for position in hotSpotMark:
                            metricsFile.write('%s\t%f\n'%(hotSpotDic[position],round(float(maxOutlier)/int(hotSpotMark[position])*100*multiplicativeFactor,2)))  
                        hotSpotMark = {}

    fileToWrite.close()
    return exonSensitivity,unCoveredRegions

def readPileupForPosition(pileup,studentValue,balance,quality,SDquality,windowSize,genomicPosition):
    previousPos = ''
    firstLine = 1
    countPosInExon = 0
    for pileupLine in pileup:
        pileupLine = pileupLine.split('\t')
        Chromosome = pileupLine[0]
        genomicPos = int(pileupLine[1])
        mutationStats = extractpileupLineInfos(pileupLine)
        rowToAdd = [countPosInExon,mutationStats['Chromosome'],int(mutationStats['GenomicPosition']),mutationStats['refBase'],mutationStats['alt'],mutationStats['DP'],mutationStats['altReadsMutation'],int(mutationStats['totalAlt']),mutationStats['AF'],mutationStats['Qual'],mutationStats['Qual_StdDev'],mutationStats['balance'],mutationStats['WT_Bal']]
        if firstLine == 1:
            exonMatrix = [rowToAdd]
            countPosInExon += 1
            firstLine = 0
            previousPos = genomicPos
        else:
            if previousPos+1 == genomicPos:
                exonMatrix.append(rowToAdd)
                previousPos = genomicPos
                countPosInExon += 1
    exonMatrix = np.asarray(exonMatrix)
    detectMutationsAtPosition(exonMatrix, studentValue, balance, quality, SDquality, windowSize, genomicPosition)

def detectMutationsAtPosition(exonMatrix,studentValue,balance,quality,SDquality,windowSize,genomicPosition):
    altReadList = exonMatrix[:,7].tolist()
    refSeqList = exonMatrix[:,3].tolist()
    for rowNum in range(0,len(exonMatrix)):
        positionToFind = '%s:%s'%(exonMatrix[rowNum,1],exonMatrix[rowNum,2])
        if positionToFind == genomicPosition:
            maxOutliers = calcConfidenceIntervalFromReads(altReadList, studentValue,int(exonMatrix[rowNum,5]))
            Chr = exonMatrix[rowNum,1]
            genPos = exonMatrix[rowNum,2]
            ref = exonMatrix[rowNum,3]
            alt = exonMatrix[rowNum,4]
            depth = exonMatrix[rowNum,5]
            AF = exonMatrix[rowNum,8]
            Qual = exonMatrix[rowNum,9]
            Qual_StdDev = exonMatrix[rowNum,10]
            balance = exonMatrix[rowNum,11]
            WTbalance = exonMatrix[rowNum,12]
            WT_F = float(exonMatrix[rowNum,12].split('/')[0])
            WT_R = float(exonMatrix[rowNum,12].split('/')[1])
            refSeqToScan = ''.join(defineWindowToScan(refSeqList, int(rowNum), 60))
            mutPos = defineMutationPosInSeq(refSeqList, int(rowNum), 60)
            motifDic = findMotif(refSeqToScan, 5, mutPos)
            if 'stretch' in motifDic.keys():
                stretchInfos = '%s (l = %i / d = %i bases from mutation)'%(motifDic['stretch'][0],len(motifDic['stretch'][0]),motifDic['stretch'][1])
            else:
                stretchInfos = 'None'
            if 'motif' in motifDic.keys():
                motifInfos = '%s (l = %i / d = %i bases from mutation)'%(motifDic['motif'][0],len(motifDic['motif'][0]),motifDic['motif'][1])
            else:
                motifInfos = 'None'
            if exonMatrix[rowNum,4] == 'WT':
                forwardPer,reversePer = 0,0
                correctedFper,correctedRper = 0,0
                overAllBalance = round((WT_F/float(depth))*100,1)
                overAllBalanceR = 100-overAllBalance
            else:
                Alt_F = float(exonMatrix[rowNum,11].split('/')[0])
                Alt_R = float(exonMatrix[rowNum,11].split('/')[1])
                if len(ref) ==1 and len(alt) ==1:
                    diffF = round((WT_F+Alt_F) / (float(depth)),3)-0.5
                    forwardPer = round((Alt_F/(Alt_F+Alt_R))*100,1)
                    reversePer = 100-forwardPer
                    overAllBalance = round(((WT_F+Alt_F)/float(depth))*100,1)
                    overAllBalanceR = 100-overAllBalance
                else:
                    diffF = round((WT_F) / (float(depth)),3)-0.5
                    forwardPer = round((Alt_F/(Alt_F+Alt_R))*100,1)
                    reversePer = 100-forwardPer
                    overAllBalance = round(((WT_F)/float(depth))*100,1)
                    overAllBalanceR = 100-overAllBalance
                if diffF>0:
                    correctedFper = round(((Alt_F-(Alt_F*diffF))/(Alt_F+Alt_R))*100,1)
                    correctedRper = 100-correctedFper
                else:
                    correctedRper = round(((Alt_R+(Alt_R*diffF))/(Alt_F+Alt_R))*100,1)
                    correctedFper = 100-correctedRper
            print altReadList
            print """
Mutation Position:         %s:%s
Reference Allele:          %s
Alternative Allele:        %s
Depth:                     %s
Allele Frequency (%%):     %s
Phred Quality:             %s
Phred Standard Deviation:  %s
Forward / Reverse alt:     %s     (%s%% / %s%%)
overAll Balance:           %s%% / %s%% 
Corrected alt F/R:         %s%% / %s%% 
Raw background Noise:      %s
Stretch nearby:            %s
Motif nearby:              %s
"""%(Chr,genPos,ref,alt,depth,AF,Qual,Qual_StdDev,balance,forwardPer,reversePer,overAllBalance,overAllBalanceR,correctedFper,correctedRper,maxOutliers,stretchInfos,motifInfos)
    if len(exonMatrix)-1 != windowSize:
        print '/!\ Warning: insufficient coverage for some positions on analyzed region.'

def writeMetricsFile(sensitivityDic,bamFile,output,bedNum):
    metricsFile = open('%s/tempAS_%i.metrics'%(output,int(bedNum)),'w')
    for interval in sensitivityDic:
        lineToWrite = '%s\t%s\n'%(interval.replace('_','\t'),sensitivityDic[interval])
        metricsFile.write(lineToWrite)
    metricsFile.close()
    
def sortVcfResults(vcfOutput,tempVcfList,tempOutput):
    vcfToWrite = open(vcfOutput,'a')
    vcfDic = {}
    for element in tempVcfList:
        vcfNum = int(element.split('.')[-2].split('_')[-1])
        vcfDic[vcfNum] = element
    sortVcfList = sorted(vcfDic.keys())
    for number in sortVcfList:
        vcfToRead = open('%s%s'%(tempOutput,vcfDic[number]),'r')
        for line in vcfToRead:
            vcfToWrite.write(line)
        vcfToRead.close()
    vcfToWrite.close()

def sortASresults(ASoutput,AStempList,tempOutput):
    mergedTempList = []
    for tempFile in AStempList:
        AStoRead = open('%s%s'%(tempOutput,tempFile),'r')
        for line in AStoRead:
            mergedTempList.append(line.replace('\n',''))
        AStoRead.close()
    sortTempFile = sorted(mergedTempList)
    AStoWrite = open(ASoutput,'w')
    for element in sortTempFile:
        AStoWrite.write(element+'\n')
    AStoWrite.close()
        
def writeHSMresults(HSMoutput,HSMtempList,tempOutput):
    HSMtoWrite = open('%s'%HSMoutput,'w')
    HSMtoWrite.write('Position_Informations\tSensitivity_Threshold(%)\n')
    for tempFile in HSMtempList:
        HSMtoRead = open('%s/%s'%(tempOutput,tempFile),'r')
        for line in HSMtoRead:
            HSMtoWrite.write(line)
        HSMtoRead.close()
    HSMtoWrite.close()
        
def findTempVcf(tempOutLyzerDir):
    allFiles = os.listdir(tempOutLyzerDir)
    tempoutLyzerResults = []
    for element in allFiles:
        if re.search('tempCalling',element):
            tempoutLyzerResults.append(element)
    return tempoutLyzerResults

def findTempAS(tempOutLyzerDir):
    allFiles = os.listdir(tempOutLyzerDir)
    tempoutLyzerResults = []
    for element in allFiles:
        if re.search('tempAS',element):
            tempoutLyzerResults.append(element)
    return tempoutLyzerResults

def findTempHSM(tempOutLyzerDir):
    allFiles = os.listdir(tempOutLyzerDir)
    tempoutLyzerResults = []
    for element in allFiles:
        if re.search('tempHSM',element):
            tempoutLyzerResults.append(element)
    return tempoutLyzerResults
    
def writeFinalResults(bamName,output,ASargs,HSMargs,commandLine):
    tempOutLyzerDir = '%s/outLyzerTemp_%s/'%(args.output,bamName)
    vcfTempList = findTempVcf(tempOutLyzerDir)
    AStempList = findTempAS(tempOutLyzerDir)
    HSMtempList = findTempHSM(tempOutLyzerDir)
    vcfOutput = '%s/%s.vcf'%(args.output,bamName)
    ASoutput = '%s/%s.metrics'%(args.output,bamName)
    HSMoutput = '%s/%s_HSM.txt'%(args.output,bamName)
    writeVcfHeader(commandLine, vcfOutput)
    sortVcfResults(vcfOutput,vcfTempList,tempOutLyzerDir)
    if 'uncoveredRegions.bed' in os.listdir(tempOutLyzerDir):
        copyfile('%s/uncoveredRegions.bed'%tempOutLyzerDir, '%s/%s_uncoveredRegions.bed'%(output,bamName))
    if ASargs:
        sortASresults(ASoutput, AStempList, tempOutLyzerDir)
    if HSMargs:
        writeHSMresults(HSMoutput, HSMtempList, tempOutLyzerDir)
    os.system('rm -r %s'%(tempOutLyzerDir))
 
def launchCallingFunction(args):
    if checkInputsCalling(args) == True:
        print 'starting outLyzer Calling Process'
        startScript = time.time()
        commandLine = ' '.join(sys.argv)
        progPath = sys.argv[0]
        bamName = args.bam.split('/')[-1].split('.')[0]
        tempOutLyzerDir = '%s/outLyzerTemp_%s/'%(args.output,bamName)
        if os.path.exists(tempOutLyzerDir):
            os.system('rm -r %s'%(tempOutLyzerDir))
        bedList = cutBedFile(args.bed, args.cut, tempOutLyzerDir)
        outLyzerCommandList = createOutLyzerCommandList(tempOutLyzerDir, args.pythonPath, args.samtools, progPath, args.cut, bedList, args.bam, args.ref, args.t, args.bal, args.Q, args.SDQ, args.WS, args.AS,args.HSM, args.x,args.verbose,args.FRcor,args.WSmin)
        manageCommandListToLaunch(outLyzerCommandList, args.core)
        if args.verbose == 1:
            print 'Formatting results for %s'%bamName
        writeFinalResults(bamName, args.output, args.AS,args.HSM,commandLine)
        endScript = time.time()
        print 'outLyzer analysis performed in %f sec.'%(endScript-startScript)
    
def launchPositionAnalysisFunction(args):
    if checkInputsPositionAnalysis(args) == True:
        print 'start analysis on position: ',args.position
        pileup = launchPileupCommandForPosition(args.bam, args.samtools, args.ref, args.position, args.WS)
        if pileup != None:
            readPileupForPosition(pileup, args.t, args.bal, args.Q, args.SDQ, args.WS, args.position)
        else:
            print 'End of Analysis'
    
def launchSubProcess(args):
    bamName = args.bam.split('/')[-1]
    bamPart = args.bed.split('/')[-1].split('.')[0].split('_')[-1]
    if args.verbose == 1:
        print "pileup conversion for %s part %s / %s"%(bamName,bamPart,args.cut)
    pileup = launchPileupCommandForBedFile(args.bam, args.samtools, args.ref, args.bed,args.output)
    bedNum = args.bed.split('.')[-2].split('_')[-1]
    if args.verbose == 1:
        print 'statistical analysis for %s part %s / %s'%(bamName,bamPart,args.cut)
    exonSensitivityLimit,uncoveredRegions = readPileup(pileup, bedNum, args.t, args.bal, args.Q, args.SDQ, args.WS, args.output, args.x, args.HSM,args.FRcor,args.WSmin)
    if args.AS:
        writeMetricsFile(exonSensitivityLimit, args.bam, args.output, bedNum)
    if len(uncoveredRegions)>0:
        uncovFile = open('%s/uncoveredRegions.bed'%args.output,'a')
        for region in uncoveredRegions:
            uncovFile.write('%s\n'%region)
        uncovFile.close()

def checkInputsCalling(args):
    checkResponse = True
    pathToCheckList = [args.bed,args.bam,args.ref,args.output]
    for pathToCheck in pathToCheckList:
        if os.path.exists(pathToCheck) == False:
            print "%s doesn't exists"%pathToCheck
            checkResponse = False
    if not args.output.endswith('/'):
        print """'/' is missing to the end of %s"""%args.output
        checkResponse = False
    if args.HSM:
        if os.path.exists('%s'%args.HSM) == False:
            print os.path.exists('%s'%args.HSM)
            print "%s doesn't exists"%args.HSM
            checkResponse = False
    return checkResponse

def checkInputsPositionAnalysis(args):
    checkResponse = True
    pathToCheckList = [args.bam,args.ref]
    for pathToCheck in pathToCheckList:
        if os.path.exists(pathToCheck) == False:
            print "%s doesn't exists"%pathToCheck
            checkResponse = False
    return checkResponse

def printLicense(args):
    print """Copyright Etienne Muller (2016)

muller.etienne@hotmail.fr

outLyzer is a computer program whose purpose is to detect low allele
frequency mutations, designed to analyze tumor samples in a clinical
context.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms."""
        
if __name__ == "__main__":
        parser = argparse.ArgumentParser(description="""
        
#               _     __                    
#    ___  _   _| |_  / / _   _ _______ _ __ 
#   / _ \| | | | __|/ / | | | |_  / _ \ '__|
#  | (_) | |_| | |_/ /__| |_| |/ /  __/ |   
#   \___/ \__,_|\__\____/\__, /___\___|_|   
#                        |___/           

OutLyzer is a variant-caller conceived for low allele-ratio mutations detection, 
based on sequencing background noise evaluation.
It evaluates if the mutation is significantly different from background noise,
using modified Thompson tau technique.""",formatter_class=RawTextHelpFormatter)
        subparsers = parser.add_subparsers(metavar='{calling,positionAnalysis,LICENSE}',help="outLyzer functions")
# create the parser for the "calling Command" command
        parser_calling = subparsers.add_parser('calling',description='outLyzer calling function, analyze the whole BAM file and displays results in a VCF (Variant Calling Format) File')
        parser_calling.add_argument('-samtools',help='Complete Samtools path if not specified in environment variable', default='samtools')
        parser_calling.add_argument('-pythonPath',help='Complete python path if different from default python version', default='python')
        parser_calling.add_argument('-core',help='define number of cores used to process analysis [1]',type = int, default = 1)
        parser_calling.add_argument('-cut',help='defines into how many parts bed file is divide [3]',type=int,default=3)
        parser_calling.add_argument('-bed',help='bed file required for analysis [REQUIRED]',required=True)
        parser_calling.add_argument('-bam',help='bam File to analyze [REQUIRED]',required=True)
        parser_calling.add_argument('-ref',help='faidx indexed reference sequence file (fasta) [REQUIRED]',required=True)
        parser_calling.add_argument('-output',help='output Path To write results [REQUIRED]',required=True)
        parser_calling.add_argument('-t',help='Student t value used in modified Thompson tau technique [0.001]',default=0.001,type=float)
        parser_calling.add_argument('-bal',help='minimum Forward / Reverse read proportion [0.3]',default=0.3,type=float)
        parser_calling.add_argument('-Q',help='minimum average Phred Score to be considered as a real mutation (only relevant for SNP) [20]',default=20,type=int)
        parser_calling.add_argument('-SDQ',help='maximum Standard deviation authorized for average Phred Score [7]',default=7,type=int)
        parser_calling.add_argument('-WS',help='Window Size: region (number of bp) around the mutation on which background noise have to be determined [200]',default=200,type=int)
        parser_calling.add_argument('-WSmin',help='Window Size Minimum Size: minimum region size (number of bp) required for analysis [10]',default=10,type=int)
        parser_calling.add_argument('-x',help='Multiplicative factor that specifies how often the mutation must be above background noise [2]',default=2,type=float)
        parser_calling.add_argument('-AS',help='Analysis sensitivity: Returns an additional file containing analysis average sensitivity for each line of bed file',action='store_true')
        parser_calling.add_argument('-FRcor',help='Forward Reverse Correction: take into account any imbalance in the Forward-Reverse reads distribution in the Forward / Reverse alternative Read Proportion (-bal option)',action='store_true')
        parser_calling.add_argument('-HSM',help='HotSpot Metrics: Produce sensitivity Threshold for HotSpot positions, in an additional file. Requires formated HotSpot File in argument (see documentation for more details).')
        parser_calling.add_argument('-verbose',help='If verbose mode is set to 1, details analysis process steps [0]',type = int, default=0)
        parser_calling.set_defaults(func=launchCallingFunction)
#create the parser for subprocess command
        parser_subprocess = subparsers.add_parser('subprocess')
        parser_subprocess.add_argument('-samtools',help='Complete Samtools path if not specified in environment variable', default='samtools')
        parser_subprocess.add_argument('-pythonPath',help='Complete python path if different from default python version', default='python')
        parser_subprocess.add_argument('-cut',help='defines into how many parts bed file is divide [3]',type=int,default=3)
        parser_subprocess.add_argument('-ref',help='faidx indexed reference sequence file (fasta)',required=True)
        parser_subprocess.add_argument('-bed',help='bed file required for analysis')
        parser_subprocess.add_argument('-bam',help='bam File to analyze',required=True)
        parser_subprocess.add_argument('-t',help='Student t value used in modified Thompson tau technique [0.001]',default=0.001,type=float)
        parser_subprocess.add_argument('-bal',help='minimum Forward / Reverse read proportion [0.3]',default=0.3,type=float)
        parser_subprocess.add_argument('-Q',help='minimum average Phred Score to be considered as a real mutation (only relevant for SNP) [20]',default=20,type=int)
        parser_subprocess.add_argument('-SDQ',help='maximum Standard deviation authorized for average Phred Score [7]',default=7,type=int)
        parser_subprocess.add_argument('-WS',help='Window Size: region (number of bp) around the mutation on which background noise have to be determined [200]',default=200,type=int)
        parser_subprocess.add_argument('-WSmin',help='Window Size Minimum Size: minimum region size (number of bp) required for analysis [10]',default=10,type=int)
        parser_subprocess.add_argument('-x',help='Multiplicative factor that specifies how often the mutation must be above background noise [2]',default=2,type=float)
        parser_subprocess.add_argument('-AS',help='Analysis sensitivity: Returns an additional file containing analysis average sensitivity for each line of bed file',action='store_true')
        parser_subprocess.add_argument('-FRcor',help='Forward Reverse Correction: take into account any imbalance in the Forward-Reverse reads distribution in the Forward / Reverse alternative Read Proportion (-bal option)',action='store_true')
        parser_subprocess.add_argument('-HSM',help='HotSpot Metrics: Produce sensitivity Threshold for HotSpot positions, in an additional file. Requires HotSpot File in argument')
        parser_subprocess.add_argument('-output',help='output Path To write results',required=True)
        parser_subprocess.add_argument('-verbose',help='If verbose mode is set to 1, details analysis process steps [0]',type = int, default=0)
        parser_subprocess.set_defaults(func=launchSubProcess)
# create the parser for the "position analysis" command
        parser_posAna = subparsers.add_parser('positionAnalysis',description='outLyzer positionAnalysis function gives an evaluation of sequencing data and local noise background for one chromosomic position')
        parser_posAna.add_argument('-samtools',help='Complete Samtools path if not specified in environment variable',default='samtools')
        parser_posAna.add_argument('-bam',help='bam File to analyze [REQUIRED]',required=True)
        parser_posAna.add_argument('-position',help='chromosomic position to analyze (ex: chr3:123456789) [REQUIRED]',required=True)
        parser_posAna.add_argument('-ref',help='faidx indexed reference sequence file (fasta) [REQUIRED]',required=True)
        parser_posAna.add_argument('-t',help='Student t value used in modified Thompson tau technique [0.001]',default=0.001,type=float)
        parser_posAna.add_argument('-bal',help='minimum Forward / Reverse read proportion [0.3]',default=0.3,type=float)
        parser_posAna.add_argument('-Q',help='minimum average Phred Score to be considered as a real mutation (only relevant for SNP) [20]',default=20,type=int)
        parser_posAna.add_argument('-SDQ',help='maximum Standard deviation authorized for average Phred Score [7]',default=7,type=int)
        parser_posAna.add_argument('-WS',help='Window Size: region (number of bp) around the mutation on which background noise have to be determined [200]',default=200,type=int)
        parser_posAna.set_defaults(func=launchPositionAnalysisFunction)
# LICENSE
        parser_LICENSE = subparsers.add_parser('LICENSE')
        parser_LICENSE.set_defaults(func=printLicense)
        args=parser.parse_args()
        args.func(args)


        