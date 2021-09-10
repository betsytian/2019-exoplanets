import csv
import os
import re
import sys

# Following are columns from Kepler Target Search result in given order
# Kepler_ID,RA (J2000),Dec (J2000),GALEX FUV,GALEX NUV,g,r,i,z,J,H,K,Teff,Data Availability,Ang Sep (')
# integer,ra,dec,float,float,float,float,float,float,float,float,float,integer,integer,float
# If columns are added/removed, change the value of KIC_NUM_COLUMNS accordingly
SWIFT_NUM_COLUMNS = 9
KIC_NUM_COLUMNS = 15
TOTAL_NUM_COLUMNS = SWIFT_NUM_COLUMNS + KIC_NUM_COLUMNS - 1

def hmsToDegree_RA(hours, minutes, seconds):
	return (float(hours)+float(minutes)/60+float(seconds)/3600)*15

def daaToDegree_DEC(degree, arcmin, arcsec):
	return float(degree)+float(arcmin)/60+float(arcsec)/3600

# In Kepler Target Search webpage, put
# input file as coordinates1.txt, select "File:comman-separated values"as "Search Output Format".
# Assume it's 1.txt.  Corresponding UV data is in uvotsource1.txt.  Based on the #.txt file's name,
# load corresponding uvotsource#.txt file.
inputFileName = sys.argv[1]
if not os.path.isfile(inputFileName):
    sys.exit("File "+inputFileName+" can't be found.")
m = re.match("kepler_fov_search_(\d+).txt", inputFileName)
fileIndex = m.group(1)
if m:
    uvotsourceFileName = "uvotsource"+fileIndex+".txt"
else:
    sys.exit("File name must be kepler_fov_search_#.txt")
if not os.path.isfile(uvotsourceFileName):
    sys.exit("File "+uvotsourceFileName+" can't be found.")

# Input line number indices in web search result start from 1.
# Add dummy row 0 in uvotsourceData to align to web search result
uvotsourceData = []
columns = [0 for x in range(TOTAL_NUM_COLUMNS)]
uvotsourceData.append(columns)
with open(uvotsourceFileName) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        columns = [0 for x in range(TOTAL_NUM_COLUMNS)]
        for i in range(len(row)):
            columns[i] = row[i]
        uvotsourceData.append(columns)

# Parse file saved from MAST Kepler Target Search
# Currently there are only 6 type of lines:
# - "Kepler_ID,...", header line lists out columns.  Skip it.
# - "integer,...", header line lists out data type of the columns.  Skip it.
# - "Input line ####:...", it marks the start of data for one source given in coordinates#.txt
# - "no rows found", there is no match of the given source
# - An empty line, this follows the "no row found" line.  Skip it.
# - comma separated data for one source that matched the given source
# sourcesNotMatched is an array storing indices of sources that do not get a match
numMatchedSources = 0
sourcesNotMatched = []
kicInfo = [[0 for x in range(KIC_NUM_COLUMNS)] for y in range(100)]
with open(inputFileName) as f:
    for linetmp in f:
        line = linetmp.rstrip()  # remove \n from the line
        if not re.search("[A-Za-z0-9]+", line):
            continue # skip empty lines
        if re.search("Kepler_ID,", line):
            continue # skip header line
        if re.search("integer,", line):
            continue # skip header line
        m1 = re.match(".*Input line\s*(\d+).*:\s*(\d+\.\d+),(\d+\.\d+).*", line)
        m2 = re.search("no rows found", line)
        if m1:
            if numMatchedSources > 0:
                numMatchedSources = 0

                # Populate the big table uvotsourceData
                # Pick the match with smallest angular difference.
                uvotsourceData[lineId][0] = kicInfo[0][0]
                for j in range(SWIFT_NUM_COLUMNS, TOTAL_NUM_COLUMNS):
                    uvotsourceData[lineId][j] = kicInfo[0][j+1-SWIFT_NUM_COLUMNS]
            lineId = int(m1.group(1), 10)

            # Sanity check : the RA and DEC in kepler_fov_search_1.txt and uvotsource1.txt should match
            lineRA = float(m1.group(2))
            lineDEC = float(m1.group(3))
            # print("line_id="+m1.group(1)+" RA="+m1.group(2)+" DEC="+m1.group(3))
            # print("RA_diff="+str(float(uvotsourceData[lineId][1]) - lineRA)+" DEC_diff="+str(float(uvotsourceData[lineId][2]) - lineDEC))
            if float(uvotsourceData[lineId][1]) - lineRA > 0.0000001 or float(uvotsourceData[lineId][2]) - lineDEC > 0.0000001:
                exit("input not matching!")
        elif m2:
            sourcesNotMatched.append(lineId)
            # Fill in blanks for row where match was not found
            numMatchedSources += 1
            for i in range(0, len(kicInfo[0])):
                kicInfo[0][i] = ""
        else:
            # Only take first 100 matches
            if (numMatchedSources < 100) :
                cols = line.split(",")
                kicInfo[numMatchedSources][0] = cols[0]
                units = cols[1].split(" ")
                kicInfo[numMatchedSources][1] = str(hmsToDegree_RA(units[0], units[1], units[2]))
                units = cols[2].split(" ")
                kicInfo[numMatchedSources][2] = str(daaToDegree_DEC(units[0], units[1], units[2]))
                for i in range(3, KIC_NUM_COLUMNS):
                    kicInfo[numMatchedSources][i] = cols[i]
                numMatchedSources += 1

    # Write out the last source
    if numMatchedSources > 0:
        uvotsourceData[lineId][0] = kicInfo[0][0]
        for j in range(SWIFT_NUM_COLUMNS, TOTAL_NUM_COLUMNS):
            uvotsourceData[lineId][j] = kicInfo[0][j + 1 - SWIFT_NUM_COLUMNS]

resultFileName = "uvotsource_result_"+fileIndex+".csv"
resultFile = open(resultFileName, "w+")
skippedFirstRow = False
count = 0
for row in uvotsourceData:
    if not skippedFirstRow:
        skippedFirstRow = True
        continue
    line = ""
    for i in range(0, len(row)):
        line += row[i]
        line += ","
    line += "\n"
    resultFile.write(line)
    count += 1
resultFile.close()
