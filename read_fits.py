import logging
import os
import re
import shutil
import sys
import time
from decimal import *
import astropy
import numpy
from astropy.io import fits
from astropy.table import Table

bkgradius_ra = 0.00361111
bkgradius_dec = 0.00361111
RA_DELTA = 0.000361111
DEC_DELTA = 0.000361111

def isLikelySameSource(pair1, pair2, radelta, decdelta):
	if abs(pair2[0] - pair1[0]) < radelta and abs(pair2[1] - pair1[1]) < decdelta:
		return True
	else:
		return False

# Generate a list of coordinates that starts from center point (5,5) and spirals out
def getSpiralCoordinates(X, Y):
	pointList = []
	x = y = 0
	dx = 0
	dy = -1
	for i in range(max(X,Y)**2):
		if (-X/2 < x <= X/2) and (-Y/2 < y <= Y/2):
			pointList.append((x+5, y+5))
		if (x==y or (x<0 and x == -y) or (x>0 and x == 1-y)):
			dx, dy = -dy, dx
		x, y = x+dx, y+dy
	return pointList

# extract data from fits file of UV sources and put RA/DEC (decimal degrees)
# into 2D array; 2D array is basically the same thing as first two data
# columns on fits file
# @param fits	name of fits file to extract RA/DEC values from
# @return	an array with two columns, RA and DEC, in decimal degrees
def FitsLocations(filename):
	hdulist = fits.open(filename)
	hdu = hdulist[1]
	data = hdu.data
	t = Table(data)
	rows = len(t)
	columns = 2
	ra_dec = [[0 for x in range(columns)] for y in range(rows)]
	for i in range(rows):
		for j in range(columns):
			ra_dec[i][j] = t[i][j]
	return ra_dec

# converts decimal degrees values to dms values for RA
def ddToDms_RA(dd):
	dd = dd/15
	degrees = int(dd)
	leftover = Decimal(dd)-Decimal(degrees)
	minutes = int(leftover*60)
	seconds = (leftover*60-minutes)*60
	return str(degrees) + ":" + str(minutes) + ":" + str(seconds)

#converts decimal degrees values to dms values for DEC
def ddToDms_DEC(dd):
	degrees = int(dd)
	leftover = Decimal(dd)-Decimal(degrees)
	minutes = int(leftover*60)
	seconds = (leftover*60-minutes)*60
	return str(degrees) + ":" + str(minutes) + ":" + str(seconds)

# Create source region file given coordinates of sources in an array
def createOneSrcRegFileFromArray(array, filename):
	totalsrcreg = open(filename, "w+")
	totalsrcreg.write("global color=green dashlist=8 3 width=3 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1")
	for pair in array:
		totalsrcreg.write("\nfk5\n")
		ra = ddToDms_RA(pair[0])
		dec = ddToDms_DEC(pair[1])
		currentregion = "circle(" + ra + "," + dec + ",3\")\n"
		totalsrcreg.write(currentregion)
	totalsrcreg.close()

# Create one source region file for one pair of coordinates of sources
def createOneSrcRegFile(pair, index, dir, sum, whichImage):
	srcreg = "srcreg" + sum + whichImage + "." + dir + "u_" + str(index) + ".reg"
	filename = os.path.join("regionresults", srcreg)
	srcreg = open(filename, "w+")
	srcreg.write("global color=green dashlist=8 3 width=3 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1")
	srcreg.write("\nfk5\n")
	ra = ddToDms_RA(pair[0])
	dec = ddToDms_DEC(pair[1])
	currentregion = "circle(" + ra + "," + dec + ",3\")\n"
	srcreg.write(currentregion)
	srcreg.close()

# Create source region files for pairs of coordinates of sources in an array
def createSrcRegFiles(array, dir):
	regcounter = 1
	for pair in array:
		createOneSrcRegFile(pair, regcounter, dir, "", "")
		regcounter = regcounter + 1

# Create background region file given coordinates of sources in an array
# Divide the region into 10x10 grids, try top to down and left to right
def createBkgRegFile(array, filename):
	bkgreg = open(filename, "w+")
	bkgreg.write("global color=green dashlist=8 3 width=3 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1")
	bkgreg.write("\nfk5\n")
	bkgradius_ra = 0.00361111
	bkgradius_dec = 0.00361111
	ra_max=array[0][0]
	ra_min=array[0][0]
	dec_max=array[0][1]
	dec_min=array[0][1]
	for pair in array:
		if pair[0]>ra_max:
			ra_max=pair[0]
		if pair[0]<ra_min:
			ra_min=pair[0]
		if pair[1]>dec_max:
			dec_max=pair[1]
		if pair[1]<dec_min:
			dec_min=pair[1]
	ra_unit=(ra_max-ra_min)/10
	dec_unit=(dec_max-dec_min)/10
	ra_bkg = 0
	dec_bkg = 0
	found = 0
	for i in range(2,9):
		for j in range(2,9):
			ra_current = ra_unit*i + ra_min
			dec_current = dec_unit*j + dec_min
			print(str(i))
			print(str(ddToDms_RA(ra_current)))
			for pair in array:
				if abs(ra_current-pair[0])>= bkgradius_ra and abs(dec_current-pair[1]) >= bkgradius_dec:
					ra_bkg = ra_unit*i + ra_min
					dec_bkg = dec_unit*j + dec_min
					found = 1
					break	
			if found == 1:
				break
		if found == 1:
			break
	ra_bkg = ddToDms_RA(ra_bkg)
	dec_bkg = ddToDms_DEC(dec_bkg)
	currentregion = "circle(" + str(ra_bkg) + "," + str(dec_bkg) + ",15\")\n"
	bkgreg.write(currentregion)
	bkgreg.close()

# Create background region file given coordinates of sources in an array
# Divide the region into 10x10 grids, try center frst, then spirals out
def createBkgRegFileFromCenter(array, filename):
	bkgreg = open(filename, "w+")
	bkgreg.write("global color=green dashlist=8 3 width=3 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1")
	bkgreg.write("\nfk5\n")
	ra_max=array[0][0]
	ra_min=array[0][0]
	dec_max=array[0][1]
	dec_min=array[0][1]
	for pair in array:
		if pair[0]>ra_max:
			ra_max=pair[0]
		if pair[0]<ra_min:
			ra_min=pair[0]
		if pair[1]>dec_max:
			dec_max=pair[1]
		if pair[1]<dec_min:
			dec_min=pair[1]
	ra_unit=(ra_max-ra_min)/10
	dec_unit=(dec_max-dec_min)/10
	ra_bkg = 0
	dec_bkg = 0
	found = 0
	pointList = getSpiralCoordinates(8, 8)
	for i, j in pointList:
		ra_current = ra_unit*i + ra_min
		dec_current = dec_unit*j + dec_min
		print(str(i))
		print(str(ddToDms_RA(ra_current)))
		for pair in array:
			if abs(ra_current-pair[0])>= bkgradius_ra and abs(dec_current-pair[1]) >= bkgradius_dec:
				ra_bkg = ra_unit*i + ra_min
				dec_bkg = dec_unit*j + dec_min
				found = 1
				break
		if found == 1:
			break
	ra_bkg = ddToDms_RA(ra_bkg)
	dec_bkg = ddToDms_DEC(dec_bkg)
	currentregion = "circle(" + str(ra_bkg) + "," + str(dec_bkg) + ",15\")\n"
	bkgreg.write(currentregion)
	bkgreg.close()

def runUvotdetect(rootDir, dir):
	infile = "sw" + dir + "u_sk.img.gz";
	expfile = "sw" + dir + "u_ex.img.gz";
	outfile = "sw" + dir + "u.fits";
	infilename = os.path.join(rootDir, dir, "uvot", "products", infile)
	if os.path.isfile(infilename):
		command = "uvotdetect" + " infile=" + infilename + " expfile=" \
				  + os.path.join(rootDir, dir, "uvot", "products", expfile)\
				  + " outfile=" + os.path.join("fitsresults", outfile) + " threshold=3"
		os.system(command)
		currentFits = os.path.join("fitsresults", outfile)
		array = FitsLocations(currentFits)
		return array
	else:
		logging.warn("IMAGE FILE NOT FOUND in " + dir)
		return []

def createDirectories():
	if os.path.isdir("fitsresults"):
		shutil.rmtree("fitsresults")
	os.makedirs("fitsresults")

	if os.path.isdir("regionresults"):
		shutil.rmtree("regionresults")
	os.makedirs("regionresults")

	if os.path.isdir("uvotsourceresults"):
		shutil.rmtree("uvotsourceresults")
	os.makedirs("uvotsourceresults")

	if os.path.isdir("uvotsourcefits"):
		shutil.rmtree("uvotsourcefits")
	os.makedirs("uvotsourcefits")

	if os.path.isdir("summedskyimages"):
		shutil.rmtree("summedskyimages")
	os.makedirs("summedskyimages")

	if os.path.isdir("summedexpimages"):
		shutil.rmtree("summedexpimages")
	os.makedirs("summedexpimages")

def oneImageToSrcBkgRegionFile(rootDir, dir):
	# Run uvotdetect on image file in dir to generate a fits file
	# Extract array of pairs of coordinates from the fits file.
	# These are coordinates of sources detected by uvotdetect
	array = runUvotdetect(rootDir, dir)
	if (len(array) == 0):
		return

	# Create a region file that has regions of all sources.
	# This file can be used by DS9 for visually check the regions
	totalsrcreg = "srcreg." + dir + "u_" + ".reg"
	filename = os.path.join("regionresults", totalsrcreg)
	createOneSrcRegFileFromArray(array, filename)

	# For each of the sources, create a region file in regionresults directory
	createSrcRegFiles(array, dir)

	# Create a background region file.
	# The createBkgRegFile searches regions top to down, left to right
	# The createBkgRegFileFromCenter searches regions from center spiraling out
	bkgreg = "bkgreg." + dir + "u.reg"
	# createBkgRegFile(array,os.path.join("regionresults", bkgreg))
	createBkgRegFileFromCenter(array, os.path.join("regionresults", bkgreg))

def imageToSrcBkgRegionFiles(rootDir, dirList):
	# Identify all sub-directory name NOT ending with '1' and
	# store the prefix in a set for later lookup. For example,
	# if sub-directory name is 00082063002, add 0008206300 to set
	excludeDirs = set()
	for dir in dirList:
		if not dir.endswith('1'):
			prefix = dir[:len(dir)-1]
			excludeDirs.add(prefix)

	# Run uvotdetect for image file in each of the sub-directories
	excludeDirsAll = set()
	dirIndex = 1
	for dir in dirList:
		# Skip sub-directories that should be excluded
		prefix = dir[:len(dir)-1]
		if (prefix in excludeDirs):
			excludeDirsAll.add(dir)
			continue
		logging.info(str(dirIndex) + " : uvotdetect " + dir)
		oneImageToSrcBkgRegionFile(rootDir, dir)
		dirIndex = dirIndex + 1
#		if (dirIndex == 3):
#			break

	return excludeDirsAll

def genUvotsourceData(rootDir):
	# Run uvotsource for all region files
	# Gather all region file sub-directory names in regdirList
	regiondir = "regionresults"
	regdirList = os.listdir(regiondir)

	# Identify all source region files and add to set
	srcregfiles = set()
	for srcFileName in regdirList:
		match = re.search("srcreg.*u_\d+\.reg", srcFileName)
		if match:
			srcregfiles.add(srcFileName)

	# Run uvotsource for all source region files and create two arrays:
	# coordinates contains RA and DEC only for all sources
	# uvotsource_data contains coordinates, flux density, and magnitude(AB)
	dirIndex = 0
	rows=len(srcregfiles)
	print(rows)
	columns=8
	uvotsource_data = [[0 for x in range(columns)] for y in range(rows)]
	coordinates = [[0 for x in range(2)] for y in range(rows)]
	for dir in srcregfiles:
		m = re.match("(\w+)\.(\w+)_(\w+)\.(\w+)", dir) # m.group(1) = srcreg, m.group(2) = 00048121001u, m.group(3) = 3, m.group(4) = reg
		swiftdir = m.group(2)  # swiftdir = dir[7:19]
		if (m.group(1).startswith("srcregsum")):
			# For images need to be summed
			n = re.match("(\D+)(\d+)",
						 m.group(1))  # m.group(1) = srcreg, m.group(2) = 00048121001u, m.group(3) = 3, m.group(4) = reg
			infile = "sum" + n.group(2) + "_" + swiftdir + "_sk.fits"
			srcfile = dir
			bkgfile = "bkgreg" + n.group(2) + "." + swiftdir + ".reg"
			imagefilename = os.path.join(rootDir, "..", "summedskyimages", infile)
		else:
			# For images without need of sum
			infile = "sw" + swiftdir + "_sk.img.gz"
			srcfile = dir
			bkgfile = "bkgreg."+swiftdir+".reg"
			imagefilename = os.path.join(rootDir, swiftdir[:len(swiftdir) - 1], "uvot", "products", infile)

		outfile = "uvotsource." + swiftdir + "_" + m.group(3) + ".fits"
		srcregfilename = os.path.join(regiondir, srcfile)
		bkgregfilename = os.path.join(regiondir, bkgfile)
		outfilename = os.path.join("uvotsourcefits", outfile)
		command = "uvotsource " + "image=" + imagefilename + " srcreg=" + srcregfilename + " bkgreg=" + bkgregfilename + " sigma=5 outfile=" + outfilename
		logging.info(str(dirIndex) + " : " + command)
		os.system(command)

		hdulist = fits.open(outfilename)
		hdu = hdulist[1]
		data = hdu.data
		t = Table(data)

		coordinates[dirIndex][0]=t[0][110] #column 1 is ra coordinate
		coordinates[dirIndex][1]=t[0][111] #column 2 is dec coordinate

		uvotsource_data[dirIndex][0] = t[0][110] #column 1 is ra coordinate
		uvotsource_data[dirIndex][1] = t[0][111] #column 2 is dec coordinate
		uvotsource_data[dirIndex][2] = t[0][45] #column 3 is AB magnitude
		uvotsource_data[dirIndex][3] = t[0][47] #column 4 is AB magnitude stat err
		uvotsource_data[dirIndex][4] = t[0][48] #column 5 is AB magnitude sys err
		uvotsource_data[dirIndex][5] = t[0][87] #column 6 is flux density
		uvotsource_data[dirIndex][6] = t[0][88] #column 7 is flux density stat err
		uvotsource_data[dirIndex][7] = t[0][89] #column 8 is flux density sys err

		dirIndex=dirIndex+1

	# Create file with all ra and dec coordinates
	# Each line has one pair of coordinates separated by a comma
	coordinate_file="coordinates.txt"
	filename = os.path.join(rootDir, "..", coordinate_file)
	coordinate_file = open(filename, "w+")
	for coordinate in coordinates:
		coordinate_file.write(str(coordinate[0]))
		coordinate_file.write(",")
		coordinate_file.write(str(coordinate[1]))
		coordinate_file.write("\n")
	coordinate_file.close()

	# Create file with all information
	# Each line has all information, separated by commas, for one source
	uvotsource_file = "uvotsource.txt"
	filename = os.path.join(rootDir, "..", uvotsource_file)
	uvotsource_file = open(filename, "w+")
	for sourceinfo in uvotsource_data:
		for i in range(0,8):
			uvotsource_file.write(str(sourceinfo[i]))
			if not i == 7:
				uvotsource_file.write(",")
		uvotsource_file.write("\n")
	uvotsource_file.close()

def sourceInArray(pair, array):
	count = 0
	for tmpPair in array:
		if (isLikelySameSource(pair, tmpPair, RA_DELTA, DEC_DELTA)):
			count += 1
			logging.warning("Found Match "+ str(count) + " : " + str(pair)+str(tmpPair)+str(abs(tmpPair[0] - pair[0]))+","+str(abs(tmpPair[1] - pair[1])))

	if (count > 0):
		return True
	else:
		return False

def getImageFrom3(pair, oneImageArrayOfSrcArray, twoImageArrayOfSrcArray, array):
	inImage1 = sourceInArray(pair, oneImageArrayOfSrcArray[0])
	inImage2 = sourceInArray(pair, oneImageArrayOfSrcArray[1])
	inImage3 = sourceInArray(pair, oneImageArrayOfSrcArray[2])
#	inImage12 = sourceInArray(pair, twoImageArrayOfSrcArray[0])
#	inImage23 = sourceInArray(pair, twoImageArrayOfSrcArray[1])
#	inImage31 = sourceInArray(pair, twoImageArrayOfSrcArray[2])
#	inImage123 = sourceInArray(pair, array)

	if (inImage1 and inImage2 and inImage3):
		return "123"
	elif (inImage1 and inImage2):
		return "12"
	elif (inImage1 and inImage3):
		return "31"
	elif (inImage3 and inImage2):
		return "23"
	elif (inImage1):
		return "1"
	elif (inImage2):
		return "2"
	elif (inImage3):
		return "3"
	else:
		#logging.warning("Source does not belong to any image!")
		return ""

def getImageFrom2(pair, oneImageArrayOfSrcArray):
	inImage1 = sourceInArray(pair, oneImageArrayOfSrcArray[0])
	inImage2 = sourceInArray(pair, oneImageArrayOfSrcArray[1])

	if (inImage1 and inImage2):
		return "12"
	elif (inImage1):
		return "1"
	elif (inImage2):
		return "2"
	else:
		#logging.warning("Source does not belong to any image!")
		return ""

# Special case when there are 3 images with the same prefix
# Get image 1, 2, 3, 12, 13, 23, and 123 ready, then
# go through the sources in summed image 123.
# If the source is in 1, 2 and 3, then, it belongs to 123
def processThreeImageSum(rootDir, prefix, dirSet):
	arrays = []
	oneImageFileNames = []
	oneExpFileNames = []
	oneUvotdetectFitsFileNames = ("1", "2", "3")
	index = 0
	oneImageArrayOfSrcArray = []
	for dir in dirSet:
		# Copy image file
		# srcFileName is used for uvotimsum, it is copied to summedskyimages
		# srcFileNameProduct is used for uvotdetect, as uvotdetect cannot
		# consume some files from image directory
		srcFileNameProduct = rootDir+"/"+dir+"/uvot/products/sw"+dir+"u_sk.img.gz"
		imageFileName = "sw"+dir+"um2_sk.img.gz"
		srcFileName = rootDir+"/"+dir+"/uvot/image/"+imageFileName
		oneImageFileNames.append(srcFileName)
		dstPath = rootDir+"/../summedskyimages/"+"sum"+str(index+1)+"_"+prefix+"u_sk.fits"
		command = "cp -f " + srcFileName + " " + dstPath
		logging.info(command)
		os.system(command)

		# Copy exposure file
		srcExpFileNameProduct = rootDir+"/"+dir+"/uvot/products/sw"+dir+"u_ex.img.gz"
		expFileName = "sw"+dir+"um2_ex.img.gz"
		srcExpFileName = rootDir+"/"+dir+"/uvot/image/"+expFileName
		oneExpFileNames.append(srcExpFileName)
		dstExpPath = rootDir+"/../summedexpimages/"+"sum"+str(index+1)+"_"+prefix+"u_ex.fits"
		command = "cp -f " + srcExpFileName + " " + dstExpPath
		logging.info(command)
		os.system(command)

		# Run uvotdetect to get fits file, then extract sources into array
		dstFitsFileName = rootDir+"/../summedskyimages/uvodetect"+oneUvotdetectFitsFileNames[index]+"."+prefix+".fits"
		command = "uvotdetect" + " infile=" + srcFileNameProduct + " expfile="\
				  + srcExpFileNameProduct + " outfile=" + dstFitsFileName + " threshold=3"
		logging.info(command)
		os.system(command)
		array = FitsLocations(dstFitsFileName)
		oneImageArrayOfSrcArray.append(array)

		# Extract the background region
		bkgRegFileName = rootDir+"/../regionresults/bkgreg"+str(index+1)+"."+prefix+"u.reg"
		createBkgRegFileFromCenter(array, bkgRegFileName)

		index = index + 1

	twoImageFileNames = ("sum12_", "sum23_", "sum31_")
	twoExpFileNames = ("sum12_", "sum23_", "sum31_")
	twoUvotdetectFitsFileNames = ("12", "23", "31")
	twoBkgRegFileNames = ("12", "23", "31")
	twoImageArrayOfSrcArray = []
	for index in range(0, 3):
		# Sum two image files
		inputFileName = oneImageFileNames[index]+","+oneImageFileNames[(index+1)%3]
		outputFileName = rootDir+"/../summedskyimages/"+twoImageFileNames[index]+prefix+"u_sk.fits"
		command = "uvotimsum infile=" + inputFileName + " outfile=" + outputFileName
		logging.info(command)
		os.system(command)

		# Sum two exposure files
		inputExpFileName = oneExpFileNames[index] + "," + oneExpFileNames[(index + 1) % 3]
		outputExpFileName = rootDir + "/../summedexpimages/" + twoExpFileNames[index]+prefix+"u_ex.fits"
		command = "uvotimsum infile=" + inputExpFileName + " outfile=" + outputExpFileName
		logging.info(command)
		os.system(command)

		# Run uvotdetect on the summed image to get fits file, then extract sources into array
		dstFitsFileName = rootDir+"/../summedskyimages/uvodetect"+twoUvotdetectFitsFileNames[index]+"."+prefix+".fits"
		command = "uvotdetect" + " infile=" + outputFileName + " expfile="\
				  + outputExpFileName + " outfile=" + dstFitsFileName + " threshold=3"
		logging.info(command)
		os.system(command)
		array = FitsLocations(dstFitsFileName)
		twoImageArrayOfSrcArray.append(array)

		# Extract the background region
		bkgRegFileName = rootDir+"/../regionresults/bkgreg"+twoBkgRegFileNames[index]+"."+prefix+"u.reg"
		createBkgRegFileFromCenter(array, bkgRegFileName)

	# Sum three image and exposure files
	inputFileName = oneImageFileNames[0]+","+oneImageFileNames[1]+","+oneImageFileNames[2]
	outputFileName = rootDir+"/../summedskyimages/sum123_"+prefix+"u_sk.fits"
	command = "uvotimsum infile=" + inputFileName + " outfile=" + outputFileName
	logging.info(command)
	os.system(command)
	inputExpFileName = oneExpFileNames[0]+","+oneExpFileNames[1]+","+oneExpFileNames[2]
	outputExpFileName = rootDir+"/../summedexpimages/sum123_"+prefix+"u_ex.fits"
	command = "uvotimsum infile=" + inputExpFileName + " outfile=" + outputExpFileName
	logging.info(command)
	os.system(command)

	# Run uvotdetect on the summed image to get fits file, then extract sources into array
	dstFitsFileName = rootDir + "/../summedskyimages/uvodetect123."+prefix+".fits"
	command = "uvotdetect" + " infile=" + outputFileName + " expfile=" \
			  + outputExpFileName + " outfile=" + dstFitsFileName + " threshold=3"
	logging.info(command)
	os.system(command)
	array = FitsLocations(dstFitsFileName)

	# Extract the background region
	bkgRegFileName = rootDir + "/../regionresults/bkgreg123." + prefix + "u.reg"
	createBkgRegFileFromCenter(array, bkgRegFileName)

	# Create a region file that has regions of all sources from the summed image.
	# This file can be used by DS9 for visually check the regions
	totalsrcreg = "srcregsum." + prefix + "u_" + ".reg"
	filename = os.path.join("regionresults", totalsrcreg)
	createOneSrcRegFileFromArray(array, filename)

	# Go through each source in summed image 123 to determine which background
	# region should be used.
	# For each of the sources, create a region file in regionresults directory
	sourceNotInImage = False
	index = 1
	for pair in array:
		# whichImage is a string indicating in which images the pair appears.
		# It could be "1", "2", "3", "12", "23", "31", or "123"
		#
		whichImage = getImageFrom3(pair, oneImageArrayOfSrcArray, twoImageArrayOfSrcArray, array)
		if (len(whichImage) > 0):
			createOneSrcRegFile(pair, index, prefix, "sum", whichImage)
			index = index + 1
		else:
			sourceNotInImage = True
	if sourceNotInImage:
		logging.warning(prefix + " : Source does not belong to any image!")

# Special case when there are 2 images with the same prefix
# Get image 1, 2, and 12 ready, then go through the sources in summed image 12
# If the source is in 1 and 2 then it belongs to 12
def processTwoImageSum(rootDir, prefix, dirSet):
	arrays = []
	oneImageFileNames = []
	oneExpFileNames = []
	oneUvotdetectFitsFileNames = ("1", "2")
	index = 0
	oneImageArrayOfSrcArray = []
	for dir in dirSet:
		# Copy image file
		# srcFileName is used for uvotimsum, it is copied to summedskyimages
		# srcFileNameProduct is used for uvotdetect, as uvotdetect cannot
		# consume some files from image directory
		srcFileNameProduct = rootDir+"/"+dir+"/uvot/products/sw"+dir+"u_sk.img.gz"
		imageFileName = "sw"+dir+"um2_sk.img.gz"
		srcFileName = rootDir+"/"+dir+"/uvot/image/"+imageFileName
		oneImageFileNames.append(srcFileName)
		dstPath = rootDir+"/../summedskyimages/"+"sum"+str(index+1)+"_"+prefix+"u_sk.fits"
		command = "cp -f " + srcFileName + " " + dstPath
		logging.info(command)
		os.system(command)

		# Copy exposure file
		srcExpFileNameProduct = rootDir+"/"+dir+"/uvot/products/sw"+dir+"u_ex.img.gz"
		expFileName = "sw"+dir+"um2_ex.img.gz"
		srcExpFileName = rootDir+"/"+dir+"/uvot/image/"+expFileName
		oneExpFileNames.append(srcExpFileName)
		dstExpPath = rootDir+"/../summedexpimages/"+"sum"+str(index+1)+"_"+prefix+"u_ex.fits"
		command = "cp -f " + srcExpFileName + " " + dstExpPath
		logging.info(command)
		os.system(command)

		# Run uvotdetect to get fits file, then extract sources into array
		dstFitsFileName = rootDir+"/../summedskyimages/uvodetect"+oneUvotdetectFitsFileNames[index]+"."+prefix+".fits"
		command = "uvotdetect" + " infile=" + srcFileNameProduct + " expfile="\
				  + srcExpFileNameProduct + " outfile=" + dstFitsFileName + " threshold=3"
		logging.info(command)
		os.system(command)
		array = FitsLocations(dstFitsFileName)
		oneImageArrayOfSrcArray.append(array)

		# Extract the background region
		bkgRegFileName = rootDir+"/../regionresults/bkgreg"+str(index+1)+"."+prefix+"u.reg"
		createBkgRegFileFromCenter(array, bkgRegFileName)

		index = index + 1

	twoImageArrayOfSrcArray = []

	# Sum two image and exposure files
	inputFileName = oneImageFileNames[0]+","+oneImageFileNames[1]
	outputFileName = rootDir+"/../summedskyimages/sum12_"+prefix+"u_sk.fits"
	command = "uvotimsum infile=" + inputFileName + " outfile=" + outputFileName
	logging.info(command)
	os.system(command)
	inputExpFileName = oneExpFileNames[0]+","+oneExpFileNames[1]
	outputExpFileName = rootDir+"/../summedexpimages/sum12_"+prefix+"u_ex.fits"
	command = "uvotimsum infile=" + inputExpFileName + " outfile=" + outputExpFileName
	logging.info(command)
	os.system(command)

	# Run uvotdetect on the summed image to get fits file, then extract sources into array
	dstFitsFileName = rootDir + "/../summedskyimages/uvodetect12."+prefix+".fits"
	command = "uvotdetect" + " infile=" + outputFileName + " expfile=" \
			  + outputExpFileName + " outfile=" + dstFitsFileName + " threshold=3"
	logging.info(command)
	os.system(command)
	array = FitsLocations(dstFitsFileName)

	# Extract the background region
	bkgRegFileName = rootDir + "/../regionresults/bkgreg12." + prefix + "u.reg"
	createBkgRegFileFromCenter(array, bkgRegFileName)

	# Create a region file that has regions of all sources from the summed image.
	# This file can be used by DS9 for visually check the regions
	totalsrcreg = "srcregsum." + prefix + "u_" + ".reg"
	filename = os.path.join("regionresults", totalsrcreg)
	createOneSrcRegFileFromArray(array, filename)

	# Go through each source in summed image 12 to determine which background
	# region should be used.
	# For each of the sources, create a region file in regionresults directory
	sourceNotInImage = False
	index = 1
	for pair in array:
		# whichImage is a string indicating in which images the pair appears.
		# It could be "1", "2", or "12"
		#
		whichImage = getImageFrom2(pair, oneImageArrayOfSrcArray)
		if (len(whichImage) > 0):
			createOneSrcRegFile(pair, index, prefix, "sum", whichImage)
			index = index + 1
		else:
			sourceNotInImage = True
	if sourceNotInImage:
		logging.warning(prefix + " : Source does not belong to any image!")

def incompleteImageToSrcBkgRegionFiles(rootDir, excludeDirsAll):
	# Dictionary excludeDirSets is like {prefixA : (prefixA1, prefixA2, prefixA3), prefixB : (prefixB1, prefixB2)}
	excludeDirSets = {}
	for dir in excludeDirsAll:
		# Skip directory if the files don't exist
		filename = os.path.join(rootDir, dir, "uvot", "image", "sw"+dir+"um2_sk.img.gz")
		if not os.path.isfile(filename):
			continue
		prefix = dir[:len(dir)-1]
		if (prefix in excludeDirSets):
			tmpSet = excludeDirSets[prefix]
			tmpSet.add(dir)
		else:
			tmpSet = set()
			tmpSet.add(dir)
			excludeDirSets[prefix] = tmpSet

	# Sum up images.  Skip image sets for which uvotdetect errors out
	prefixToSkip = set(["0004814800"])
	count3 = count2 = count1 = 0
	for prefix, dirSet in excludeDirSets.items():
		print(prefix)
		print(dirSet)
		if (prefix in prefixToSkip):
			continue
		if (len(dirSet) == 3):
			processThreeImageSum(rootDir, prefix, dirSet)
			count3 += 1
		elif (len(dirSet) == 2):
			processTwoImageSum(rootDir, prefix, dirSet)
			count2 += 1
		else:
			count1 += 1
			for dir in dirSet:
				oneImageToSrcBkgRegionFile(rootDir, dir)
	logging.info("Set of 3 images : "+str(count3))
	logging.info("Set of 2 images : "+str(count2))
	logging.info("Set of 1 images : "+str(count1))

def generateData(rootDir):
	# Create all result directories
	createDirectories()

	# Gather names of all sub-directory under rootDir in dirList
	dirList = os.listdir(rootDir)

	# Identify all sub-directory name NOT ending with '1' and
	# store the prefix in excludeDirs for later lookup. For example,
	# if sub-directory name is 00082063002, add 0008206300 to excludeDirs
	# Run uvotdetect on images to get fits file.  Extract coordinates of
	# sources and create source region files and background region file
	excludeDirsAll = imageToSrcBkgRegionFiles(rootDir, dirList)

	# Special handling for incomplete images.  However, if there is only one good image among
	# the directories sharing same prefix, treat it as a single complete image by adding the
	# source region and background region files to regionresults directory like imageTosrcBkgRegionFiles
	incompleteImageToSrcBkgRegionFiles(rootDir, excludeDirsAll)

	# Run uvotsource for all region files
	# Gather all region file sub-directory names in regdirList
	genUvotsourceData(rootDir)

# Main script
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

if len(sys.argv) != 2:
	print(sys.argv)
	sys.exit("Usage : python run.py <Swift data directory name>")

if not os.path.isdir(sys.argv[1]):
	message = str(sys.argv[1]) + " is not a directory"
	sys.exit(message)

startTime = time.time()
generateData(sys.argv[1])
endTime = time.time()
print("wall time : " + str(int(endTime - startTime)) + "s")
