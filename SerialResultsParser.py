#### Hello! I'm a python script who's here to post-process your Music Post output files!
#### My goal is to read in your Post file (identified by its suffix e.g. 'postfile') to an array;
#### Identify the number of simulation steps, then find the start point of each in the file;
#### I'll then read in the three raw data outputs for each step: energy vs iteration, number vs iteration, and specific interation;
#### I'll analyse these inputs to work out the step statistics and some other useful bits
#### Finally, I'll make a directory with .csv files ready to be graphed/further processed!
####
#### To do this, your .post files will need to follow a rigid data pattern:
#### ## Total Energy averages , for the ensemble
#### ## configfile : ./{{Tag}}.{{Step}}
#### ~~~Data~~~ (iter, inst. energy, block avg energy, cum avg energy)
#### (blank line)
#### ## All the averages below are based on cum avg
#### ## Unstored/not used nrgs might appear as zero
#### Sorbs: Nrg-Total mol/uc Nrgs(KJ/mol)
#### ~~~Data~~~ (i--j pair (max 20 chars, just i if intramolecular), interaction type, ':', cum avg E, cum avg N, kJ/mol)
#### (NB data can also be: i--j pair, '-- No Molecs --' if N =0)
#### (blank line)
#### ## Per Unitcell Loading averages ,  for the ensemble
#### ## configfile : ./{{Tag}}.{{Step}}
####   Iter. No    Inst.(molec/uc)   Block (molec/uc)   Cumul(molec/uc
#### Molecule : {{Sorbent}}
#### ~~~Data~~~
#### ## Per Unitcell Loading averages ,  for the ensemble
#### ## configfile : ./{{Tag}}.{{Step}}
####  Iter. No    Inst.(molec/uc)   Block (molec/uc)   Cumul(molec/uc 
#### Molecule : {{Sorbate}}
#### FIXED CONFIG, SO NO LOADING AVG REPORTED FOR THIS SORB
#### (blank line)
#### (Next step begins)

import os

species = 'Nitrogen'
sorbent = 'IRMOF1'
path = './'
tag = 'postfile'



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


######
The script bit
######

filelist = Findfiles(path, tag)
for i in filelist:
    with open(i, 'r') as f:
	    alldata = f.readlines()
	confiles = DetectNSteps(alldata, tag)
	lines = []
	Elines = []
	Ilines = []
	Nlines = []
	for j in confiles:
	    val1 = DetectStepBegin(alldata, j)
		lines.append(val1)
	for k in lines:
        val2 = FindEnergyDataBegin(alldata, k[0])
		Elines.append(val2)
		val3 = FindIntDataBegin(alldata, k[0])
		Ilines.append(val3)
		val4 = FindNDataBegin(alldata, k[1])
		Nlines.append(val4)
    directorymaker('./processed_results/')
	for count,l in Elines:
	    iter,nrg = RawXvsIter(alldata,l)
		with open('./processed_results/{0}.{1}.{2:02d}.energytraj.csv'.format(species,sorbent,count), 'w') as f:
		    f.write(
'''##Raw energy vs iteration data output data from simulation {0} on {1}
Iter, Etot (kJ/mol)'''.format(species, sorbent))
            for i in range(len(iter)):
			    f.write('{0}, {1}\n'.format(iter[i],nrg[i]))
	for count,m in Nlines:
	    iter,N = RawXvsIter(alldata,l)
		with open('./processed_results/{0}.{1}.{2:02d}.occupancytraj.csv'.format(species,sorbent,count), 'w') as f:
		    f.write(
'''##N vs iteration data output data from simulation {0} on {1}
Iter, Etot (kJ/mol)'''.format(species, sorbent))
            for i in range(len(iter)):
			    file.write('{0}, {1}\n'.format(iter[i],N[i]))
	NCFluFlu = []
	CoulFluFlu = []
	NCFluFra = []
	CoulFluFra = []
	for n in Ilines:
	    a,b,c,d = RawInt(alldata,n,species, sorbent)
		NCFluFlu.append(a)
		CoulFluFlu.append(b)
		NCFluFra.append(c)
		CoulFluFra.append(d)
	with open('./processed_results/{0}.{1}.all.SimVsInt.csv'.format(species,sorbent), 'w') as f:
	    f.write(
'''## Simulation step vs cum avg. interaction energy from simulation of {0} on {1}
Step, NC fluid-fluid (kJ/mol), Coul fluid-fluid, NC fluid-framework, Coul fluid-framework'''.format(species, sorbent))
        for o in range(len(NcFluFlu)):
            f.write('{0}, {1}, {2}, {3}, {4}\n'.format(o, NCFluFlu[o], CoulFluFlu[o], NcFluFra[o], CoulFluFra[o]))