
import numpy as np
import os

species = 'Nitrogen'
sorbent = 'IRMOF1'
path = './'
tag = 'postfile'

def isothermextract(page, species): #this function scrapes subdirectories of name for isotherm data for a given isotherm file name
    isotherms = np.zeros((17, 20)) # this is a matrix of a 16x 20-point isotherms with identical pressure values
    z = 0
    z = int(page)
    print(z)
    drx = str("./{0:02d}/isotherm.{1}".format(page,species)) # defines the isotherm
    print(drx)
    f = open(drx, "r")
    x = []
    page = []
    for curline in f:
        if "#" not in curline:
            fish = curline.split()
            if len(fish) == 2:
                x.append(fish[0])
                page.append(fish[1])
    f.close()
    isotherms[0] = x
    isotherms[z] = page
    return isotherms
#####################################################################################################################################################	


def FindFiles(path = './', tag = 'postfile'):
    postfiles = []
	for filename in os.listdir(path): ### find your post files
        print(filename)
        if filename.endswith(tag):
            print('Identified post file ', filename, ' for analysis.')
            postfiles.append(filename)
	return postfiles

def DetectNSteps(filetext, Template):
    alldata = []
	confiles = set()
    for line in filetext:
        if Template in line:
            confiles.add(line.split(' ')[-1].split('.')[-1])
	return confiles

def DetectStepBegin(filetext, confile):
    startlines = []
    for count,line in enumerate(filetext):
        if confile in line:
            startlines.append(count)
	return startlines[0:1] 			
	### WARNING: this line assumes the first 2 mentions of your con file name are immediately before iter/Energy and iter/N data

def FindEnergyDataBegin(filetext, line):
    if 'Total Energy averages' in filetext(line-1):
	    return line+2
    else:
	    print('Oh no! I failed to find your energy data for this step!')
		return False

def FindIntDataBegin(filetext, line):
    while True:
        line +=1
        if len(alldata[line]) == 0:
            break
	line +=1
	if 'All the averages below are based on cum avg' in filetext(line):
		return line+3
	else:
	    print('Oh no! I failed to find your interactions data for this step!')
		return False

def FindNDataBegin(filetext,line):
    if 'Per Unitcell Loading averages' in filetext(line-1):
	    return line+2
    else:
	    print('Oh no! I failed to find your loading data for this step!')
		return False

def RawXvsIter(filetext,line):
    iter = []
	val = []
	while True:
	    if len(filetext[line]) ==0:
		    break
		linedata = filetext[line].split(' ')
		iter.append(linedata[0])
		val.append(linedata[2])
		line+=1
	return iter, val
	
def RawInt(filetext,line, species, sorbent):
    NCFluFlu = None
	CoulFluFlu = None
	NCFluFra = None
	CoulFluFra = None
	while True:
	    rawdata = []
	    if len(filetext[line]) == 0:
		    break
		elif 'No Molecs' in filetext[line]:
		    line +=1
		elif len(filetext[line].split' ') !=6:
		    print('!!!\nSomething weird just happened, and I don\'t recognise this line. Looks like the RawInt function needs to be fixed!')
			print('Linedata:\n',filetext[line],'\n!!!')
		else:
		    rawdata.append(i for i in filetext[line].split' ')
			if rawdata[1] == 'Coulombic':
			    if rawdata[0].split('--')[1] in species:
				    CoulFluFlu = rawdata[5]
				elif rawdata[0].split('--')[1] in sorbent:
				    CoulFluFra = rawdata[5]
				else:
				    print('!!!\nSomething weird just happened, and I don\'t recognise this line. Looks like the RawInt function needs to be fixed!')
			        print('Data read in:\n',rawdata,'\n!!!')
			elif rawdata[1] == 'Coulombic':
			    if rawdata[0].split('--')[1] in species:
				    CoulFluFlu = rawdata[5]
				elif rawdata[0].split('--')[1] in sorbent:
				    CoulFluFra = rawdata[5]
				else:
				    print('!!!\nSomething weird just happened, and I don\'t recognise this line. Looks like the RawInt function needs to be fixed!')
			        print('Data read in:\n',rawdata,'\n!!!')
            else:
			    print('!!!\nSomething weird just happened, and I don\'t recognise this line. Looks like the RawInt function needs to be fixed!')
			    print('Data read in:\n',rawdata,'\n!!!')
		line +=1
	return NCFluFlu,CoulFluFlu,NCFluFra,CoulFluFra	

def directorymaker(dxout = "./"):
    filename = "{0}test.txt" .format(dxout) #Test file name
    if not os.path.exists(os.path.dirname(filename)): #Checks if the test file exists
        try:
            os.makedirs(os.path.dirname(filename)) #Makes the file
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    with open(filename, "w") as f:
        f.write("FOOBAR") #Writes something in the test file
    print("Directory {0} written!".format(dxout))

