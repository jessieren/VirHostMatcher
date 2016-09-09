#!/usr/bin/env python
import os,sys
import optparse
import subprocess
from subprocess import call
import time
import platform

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-v", "--virusFaDir", action = "store", type = "string", dest = "virusFaDir",
									help = "the directory to the folder containing bacteria virus fasta files")
parser.add_option("-b", "--hostFaDir", action = "store", type = "string", dest = "hostFaDir",
									help = "the directory to the folder containing bacteria host fasta files")
parser.add_option("-o", "--out", action = "store", type = "string", dest = "outDir",
									default='./', help = "output directory")
parser.add_option("-t", "--taxa", action = "store", type = "string", dest = "hostTaxaFile",
									help = "the host taxa file (including the path) ")
parser.add_option("-d", "--d2star", action = "store", type = "string", dest = "onlyD2star",
									default=0, help = "compute only d2star dissimilarity? 1 for yes, 0 for no")
#parser.add_option("-k", "--kLen", action = "store", type = "string", dest = "kLen",
#									help = "the length of k-tuple")

(options, args) = parser.parse_args()
if (options.virusFaDir is None or
		options.hostFaDir is None) :
	sys.stderr.write(prog_base + ": Error: missing required command-line argument")
	parser.print_help()
	sys.exit(0)


nameLen = 93 - len(options.outDir)
## possibly because of the kmercount folder name for each contig is too long?

#################### 0: preparation ############################
## path to the programs
vhmPath  = os.path.dirname(sys.argv[0])
if len(vhmPath) == 0 :
	vhmPath="./"
#print vhmPath

## kmer length and MC order
kmax = 6
order = 2

## compile c++ code if not
optSys = platform.system()
if optSys == 'Linux' :
	exePath = os.path.join(vhmPath, "bin", "linux64")
elif optSys == 'Darwin' :
	exePath = os.path.join(vhmPath, "bin", "macDarwin")
elif optSys == 'Windows' :
	exePath = os.path.join(vhmPath, "bin", "windows64")
else :
	sys.stderr.write("Error: can't recognize the operating system" + optSys + " \n")

## copy exe files to vhmPath
os.system("cp " + os.path.join(exePath, "*") + " " + vhmPath)


## countKmer c++ code
countKmerCpp = os.path.join(vhmPath, "countKmer.cpp")
countKmerOut = os.path.join(vhmPath, "countKmer.out")
if not os.path.exists(countKmerOut) :
	sys.stderr.write( "Warning: can't find file " + countKmerOut + "\n")
	sys.stderr.write( "Trying to compile..." + "\n")
	if os.path.exists(countKmerCpp) :
		os.system("g++ " + countKmerCpp + " -o " + countKmerOut)
	else :
		sys.stderr.write( "Error: can't find file " + countKmerCpp + "\n")
		sys.exit(0)
else :
	os.system("chmod 777 " + countKmerOut)

## computeMeasure c++ code
computeMeasureCpp = os.path.join(vhmPath, "computeMeasure.cpp")
computeMeasureOut = os.path.join(vhmPath, "computeMeasure.out")
if not os.path.exists(computeMeasureOut) :
	sys.stderr.write( "Warning: can't find file " + computeMeasureOut + "\n")
	sys.stderr.write( "Trying to compile..." + "\n")
	if os.path.exists(computeMeasureCpp) :
		os.system("g++ " + computeMeasureCpp + " -o " + computeMeasureOut)
	else :
		sys.stderr.write( "Error: can't find file " + computeMeasureCpp + "\n")
		sys.exit(0)
else :
	os.system("chmod 777 " + computeMeasureOut)

## computeMeasure c++ code
computed2starCpp = os.path.join(vhmPath, "computeMeasure_onlyd2star.cpp")
computed2starOut = os.path.join(vhmPath, "computeMeasure_onlyd2star.out")
if not os.path.exists(computed2starOut) :
	sys.stderr.write( "Warning: can't find file " + computed2starOut + "\n")
	sys.stderr.write( "Trying to compile..." + "\n")
	if os.path.exists(computed2starCpp) :
		os.system("g++ " + computed2starCpp + " -o " + computed2starOut)
	else :
		sys.stderr.write( "Error: can't find file " + computed2starCpp + "\n")
		sys.exit(0)
else :
	os.system("chmod 777 " + computed2starOut)




#countKmerCpp = os.path.join(vhmPath, "countKmer.cpp")
#countKmerOut = os.path.join(vhmPath, "countKmer.out")
##print countKmerOut
#
##if not os.path.exists(countKmerOut) and os.path.exists(countKmerCpp) :
##	os.system("g++ " + countKmerCpp + " -o " + countKmerOut)
##elif not os.path.exists(countKmerCpp) :
##	sys.stderr.write( "can't find countKmer.cpp file in " + vhmPath + "\n")
##	sys.exit(0)
#if os.path.exists(countKmerCpp) :
#	os.system("g++ " + countKmerCpp + " -o " + countKmerOut)
#else :
#	sys.stderr.write( "Error: can't find countKmer.cpp file in " + vhmPath + "\n")
#	sys.exit(0)


#computeMeasureCpp = os.path.join(vhmPath, "computeMeasure.cpp")
#computeMeasureOut = os.path.join(vhmPath, "computeMeasure.out")
##if not os.path.exists(computeMeasureOut) and os.path.exists(computeMeasureCpp) :
##  os.system("g++ " + computeMeasureCpp + " -o " + computeMeasureOut)
##elif not os.path.exists(computeMeasureCpp) :
##	sys.stderr.write( "can't find computeMeasure.cpp file in " + vhmPath + "\n")
##	sys.exit(0)
#if os.path.exists(computeMeasureCpp) :
#	os.system("g++ " + computeMeasureCpp + " -o " + computeMeasureOut)
#else :
#	sys.stderr.write( "Error: can't find computeMeasure.cpp file in " + vhmPath + "\n")
#	sys.exit(0)
#
#computed2starCpp = os.path.join(vhmPath, "computeMeasure_onlyd2star.cpp")
#computed2starOut = os.path.join(vhmPath, "computeMeasure_onlyd2star.out")
##if not os.path.exists(computed2starOut) and os.path.exists(computed2starCpp) :
##	os.system("g++ " + computed2starCpp + " -o " + computed2starOut)
##elif not os.path.exists(computed2starCpp) :
##	sys.stderr.write( "can't find computeMeasure_onlyd2star.cpp file in " + vhmPath + "\n")
##	sys.exit(0)
#if os.path.exists(computed2starCpp) :
#	os.system("g++ " + computed2starCpp + " -o " + computed2starOut)
#else :
#	sys.stderr.write( "Error: can't find computeMeasure_onlyd2star.cpp file in " + vhmPath + "\n")
#	sys.exit(0)
#

if int(options.onlyD2star) == 1 :
	computeMeasureOut=computed2starOut
else :
	computeMeasureOut=computeMeasureOut

## tmp file directory
if not os.path.exists(options.outDir) :
	os.makedirs(options.outDir)
tmpDir = os.path.join(options.outDir, "tmp")
if not os.path.exists(tmpDir) :
	os.makedirs(tmpDir)

## kmer count directory
kmerCountPath = os.path.join(tmpDir, "KC")
if not os.path.exists(kmerCountPath) :
	os.makedirs(kmerCountPath)

## virusFaList, hostFaList
virusFaList = os.listdir(options.virusFaDir)
#virusFaList = os.path.join(options.virusFaDir, os.listdir(options.virusFaDir))
hostFaList = os.listdir(options.hostFaDir)

## virus list file, host list file
virusListFile = os.path.join(tmpDir, "virusList")
hostListFile = os.path.join(tmpDir, "hostList")

virusListFileWrite = open(virusListFile, 'w') ## make file blank
virusListFileWrite.close()
virusListFileWrite = open(virusListFile, 'a')

hostListFileWrite = open(hostListFile, 'w') ## make file blank
hostListFileWrite.close()
hostListFileWrite = open(hostListFile, 'a')

#time.sleep(6) # delays for 10 seconds

#################### 00: if hostTaxa missing ###################
if options.hostTaxaFile is None :
	sys.stdout.write("no hostTaxa file provided, creating a dummy one \n")
	hostTaxaFile = os.path.join(options.outDir, "hostTaxa.txt")
	hostTaxaFileWrite = open(hostTaxaFile, 'w') ## make file blank
	hostTaxaFileWrite.close()
	hostTaxaFileWrite = open(hostTaxaFile, 'a')
	
	hostTaxaFileWrite.write("hostNCBIName	hostName	hostSuperkingdom	hostPhylum	hostClass	hostOrder	hostFamily	hostGenus	hostSpecies\n")
	for currentFileName in hostFaList :
		if currentFileName.startswith('.') :
			continue
		if len(currentFileName) > nameLen :
			sys.stdout.write( " the file name has more than " + str(nameLen) + " letters! Use the first " + str(nameLen) + " letters as the name \n")
			currentFileNameS = currentFileName[:nameLen]
		else :
			currentFileNameS = currentFileName
		hostTaxaFileWrite.write(currentFileNameS)
		for i in range(1,9) :
			hostTaxaFileWrite.write("\t" + "unknown")
		hostTaxaFileWrite.write("\n")
	hostTaxaFileWrite.close()
	options.hostTaxaFile = hostTaxaFile



#################### 1: count kmer and prepare list files ############################
#sys.stdout.write("Step 1: counting kmers \n")
for currentFileName in virusFaList :
	if currentFileName.startswith('.') :
		continue
	if len(currentFileName) > nameLen :
		currentFileNameS = currentFileName[:nameLen]
	else :
		currentFileNameS = currentFileName
	sys.stdout.write("Step 1: counting kmers for virus " + currentFileNameS + "\n")
	for w in range(1, (kmax+1)) :
		currentFilePath = os.path.join(options.virusFaDir, currentFileName)
		currentKmerCountPath = os.path.join(kmerCountPath, currentFileNameS)
		cmdKmer = countKmerOut + " -l -k " + str(w) + \
			" -i " + currentFilePath +\
				" -o " + currentKmerCountPath +\
					" -s " + currentFileNameS
		cmdKmerOut = subprocess.Popen(cmdKmer, shell=True, \
																	stderr = subprocess.PIPE, \
																	stdout = subprocess.PIPE)
cmdKmerOut.wait()
	
	if len(os.listdir(currentKmerCountPath)) == ( kmax + 1 ):
		virusListFileWrite.write(currentFileNameS + " " + \
														 currentKmerCountPath + " " +\
														 str(2) + "\n")
	else :
		sys.stderr.write( "Error in counting kmers for " + currentFileNameS + "\n")
		sys.exit(0)

virusListFileWrite.close()


for currentFileName in hostFaList :
	if currentFileName.startswith('.') :
		continue
	if os.path.isdir(os.path.join(options.hostFaDir, currentFileName)) :
		sys.stderr.write( "Error: zero bytes of file " + currentFileNameS + "\n")
		sys.exit(0)
	if len(currentFileName) > nameLen :
		currentFileNameS = currentFileName[:nameLen]
	else :
		currentFileNameS = currentFileName
	sys.stdout.write("Step 1: counting kmers for host " + currentFileNameS + "\n")
	for w in range(1, (kmax+1)) :
		currentFilePath = os.path.join(options.hostFaDir, currentFileName)
		currentKmerCountPath = os.path.join(kmerCountPath, currentFileNameS)
		cmdKmer = countKmerOut + " -l -k " + str(w) + \
			" -i " + currentFilePath +\
				" -o " + currentKmerCountPath +\
					" -s " + currentFileNameS
		#print cmdKmer
		cmdKmerOut = subprocess.Popen(cmdKmer, shell=True, \
																	stderr = subprocess.PIPE, \
																	stdout = subprocess.PIPE)
cmdKmerOut.wait()
	
	if len(os.listdir(currentKmerCountPath)) == ( kmax + 1 ) :
		hostListFileWrite.write(currentFileNameS + " " + \
														currentKmerCountPath + " " +\
														str(2) + "\n")
	else :
		sys.stderr.write( "Error in counting kmers for " + currentFileNameS + "\n")
		sys.exit(0)
hostListFileWrite.close()



#time.sleep(6) # delays for 10 seconds

################### 2: compute measures #####################
sys.stdout.write("Step 2: compute distance/dissimialrity measures \n")
cmdCptMeasure = computeMeasureOut + " -k " + str(kmax) + \
	" -i " + virusListFile + " -j " + hostListFile + " -o " + options.outDir + " -t " + options.hostTaxaFile
#print cmdCptMeasure
#print cmdCptMeasure
cmdCptMeasureOut = subprocess.Popen(cmdCptMeasure, shell=True, \
																		stderr = subprocess.PIPE, \
																		stdout = subprocess.PIPE)
for line in iter(cmdCptMeasureOut.stderr.readline, b''):
	sys.stdout.write(line)
cmdCptMeasureOut.wait()






