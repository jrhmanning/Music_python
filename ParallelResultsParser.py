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
#### Molecule : {{framework}}
#### ~~~Data~~~
#### ## Per Unitcell Loading averages ,  for the ensemble
#### ## configfile : ./{{Tag}}.{{Step}}
####  Iter. No    Inst.(molec/uc)   Block (molec/uc)   Cumul(molec/uc 
#### Molecule : {{Sorbate}}
#### FIXED CONFIG, SO NO LOADING AVG REPORTED FOR THIS SORB
#### (blank line)
#### (Next step begins)

import os
import sys
import numpy as np
from pymbar import timeseries
from collections import defaultdict

species = 'Nitrogen'
framework = 'IRMOF1step3'
path = './postfiles/'
tag = 'postoutput'
Template='{0}.{1}.con'.format(framework, species)
outputpath='./processed_results/'

def FindFiles(path = './', tag = 'postfile'):
    postfiles = []
    for filename in os.listdir(path): ### find your post files
        print(filename)
        if filename.endswith(tag):
            print('Identified post file ', filename, ' for analysis.')
            postfiles.append(filename)
    return postfiles

def DetectNSteps(filetext, Template):
    confiles = set()
    for line in filetext:
        if Template in line:
            confiles.add(line.split(' ')[-1].split('.')[-1])
    print('{0} steps'.format(len(list(confiles))))
    #print(confiles)
    return confiles

def DetectStepBegin(filetext, confile):
    startlines = []
    print(confile)
    for count,line in enumerate(filetext):
        if confile in line:
            startlines.append(count)
    print(startlines)
    return startlines[0:2]
    ### WARNING: this line assumes the first 2 mentions of your con file name are immediately before iter/Energy and iter/N data

def FindEnergyDataBegin(filetext, line):
    if 'Total Energy averages' in filetext[line-1]:
        return line+2
    else:
        print('Oh no! I failed to find your energy data for this step!')
    return False

def FindIntDataBegin(filetext, line):
    while True:
        line +=1
        if len(filetext[line].strip()) == 0:
            break
    line +=1
    if 'All the averages below are based on cum avg' in filetext[line]:
        return line+3
    else:
        print('Oh no! I failed to find your interactions data for this step!')
    return False

def FindNDataBegin(filetext,line):
    if 'Per Unitcell Loading averages' in filetext[line-1]:
        return line+2
    else:
        print('Oh no! I failed to find your loading data for this step!')
        return False

def RawXvsIter(filetext,line):
    XvIDict={}
    while True:
        linedata=filetext[line].strip()
        if len(linedata) ==0:
            break
        if linedata.startswith('##'):
            break
        linedata2=linedata.split()
        #print(linedata2)
        XvIDict[int(linedata2[0])]=float(linedata2[2])
        line+=1
    #print(XvIDict)
    return XvIDict

def IntReader(linedata, species, framework):
    FluFlu=False
    FluFra=False
    Coul=False
    NCoul=False
    Value=None
    rawdata=linedata.split()
    #print(rawdata)
    if '--' in rawdata[0]:
        tag1=rawdata[0].split('--')[0]
        tag2=rawdata[0].split('--')[1]
    else:
        return FluFlu, FluFra, Coul, NCoul, Value
    if tag1 not in species:
        return FluFlu, FluFra, Coul, NCoul, Value 
    if not tag2:
        print('!!!\nI couldn\'t find a second species in your interaction at line {0}. \nCheck the post file at line {0}, and if that\'s not the problem, check the IntReader function!\n!!!'.format(line))
        sys.exit()
    if tag2 in species:
        FluFlu=True
    elif tag2 in framework:
        FluFra=True
    else:
        print('!!!\nI couldn\'t identify the second species in your interaction at line {0}. \nCheck the post file at line {0}, and if that\'s not the problem, check your variable definitions, then the IntReader function!\n!!!'.format(line))
        sys.exit()
    if rawdata[1] == 'Coulombic':
        Coul=True
    elif rawdata[1] =='NonCoulom':
        NCoul=True
    else:
        return FluFlu, FluFra, Coul, NCoul, Value
    Value=rawdata[5]
    return FluFlu, FluFra, Coul, NCoul, Value
	
def RawInt(filetext,line, species, framework):
    ResultDict={
	'NCFluFlu':None,
	'CoulFluFlu':None,
	'NCFluFra':None,
	'CoulFluFra':None
	}
    Value=None
    while True:
        linedata=filetext[line].strip()
        #print(linedata)
        if len(linedata) == 0:
            break
        elif 'No Molecs' in linedata:
            line +=1
            continue
        elif len(linedata.split()) <5:
            print('!!!\nSomething weird just happened, and I don\'t recognise this line. Looks like the RawInt function needs to be fixed!')
            print('Linedata:\n',linedata,'\n!!!')
        else:
            FluFlu, FluFra, Coul, NCoul, Value = IntReader(linedata, species, framework)
        if Value is None:
            line+=1
            continue
        if FluFlu:
            if Coul:
                ResultDict['CoulFluFlu']=[float(Value)]
            if NCoul:
                ResultDict['NCFluFlu']=[float(Value)]
        if FluFra:
            if Coul:
                ResultDict['CoulFluFra']=[float(Value)]
            if NCoul:
                ResultDict['NCFluFra']=[float(Value)]
        line +=1
    return ResultDict	

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

def collatedata(dictionary):
#I want this function to read in my huge data dictionaries 
#Then I want it to work out the approximate final value (U or N) for each markov chain and simulation step
#Then i want to average each of these approximate final values across parallel markov chains 
    simsteps=list(dictionary.keys())
    print(simsteps)
    timesteps=list(dictionary[simsteps[0]].keys())
    print(timesteps)
    parallelsims=len(dictionary[simsteps[0]][timesteps[0]])
    print(dictionary[simsteps[0]][timesteps[0]])
    resultstot={}
#    datalist=np.zeros(len(EDicts[simsteps[0]].keys()))
    for t in simsteps:
        resultstot[t]={'Ave':[], 'Std':[]}
    for q in range(parallelsims):
        print(q)
        for r in simsteps: 
            #print(len(datalist))
            placeholder=[]
            for s in timesteps:
                placeholder.append(dictionary[r][s][q])
            #print(placeholder)
            datalist=np.asarray(placeholder, dtype=float)
            #print(datalist)
            [t0, g, Neff_max] = timeseries.detectEquilibration(datalist)
            #print('t0={0}'.format(t0))
            #print(datalist[t0:])
            avg=np.mean(datalist[t0:])
            sdv=np.std(datalist[t0:])
            resultstot[r]['Ave'].append(round(avg,3))
            resultstot[r]['Std'].append(round(sdv,3))
    return resultstot


######
#The data gathering bit
######

filelist = FindFiles(path, tag)
print(filelist)
EDicts=defaultdict(list)
NDicts=defaultdict(list)
IntDicts=defaultdict(list)
print(IntDicts)
for i in range(len(filelist)):
    with open('{0}{1}'.format(path, filelist[i]), 'r') as f:
        print('~~~~~reading file {0}~~~~~'.format(filelist[i]))
        alldata = f.readlines()
    confiles = DetectNSteps(alldata, Template)
    print('I\'ve found {0} steps in your simulation data - please check this is right!'
          .format(len(list(confiles))))
    lines = []
    Elines = []
    Ilines = []
    Nlines = []
    for j in range(1,len(confiles)+1):
        val1= DetectStepBegin(alldata, '{0}.{1}'.format(Template,j))
        lines.append(val1)
    #print(lines)
    for k in lines:
        #print(k)
        val2 = FindEnergyDataBegin(alldata, k[0])
        Elines.append(val2)
        val3 = FindIntDataBegin(alldata, k[0])
        Ilines.append(val3)
        val4 = FindNDataBegin(alldata, k[1])
        Nlines.append(val4)

        for count,l in enumerate(Elines, 1):
            if count not in EDicts:
                EDicts[count]=defaultdict(list)
            NEWEDict = RawXvsIter(alldata,l)
        for key in NEWEDict.keys():
                #print(key)
            EDicts[count][key]
            EDicts[count][key].append(NEWEDict[key])
        #print(EDicts)
        for count, m in enumerate(Nlines, 1):
            if count not in NDicts:
                NDicts[count]=defaultdict(list)
            NEWNDict = RawXvsIter(alldata,m+1)
        for key in NEWNDict.keys():
                #print(key)
            NDicts[count][key]
            NDicts[count][key].append(NEWNDict[key])

        for count, n in enumerate(Ilines, 1):
            if count not in IntDicts:
            #print(IntDicts)
                IntDicts[count]=defaultdict(list)
            NEWIntDict=RawInt(alldata, n, species, framework)
            #print(NEWIntDict)
        #print(IntDicts)
        for key in NEWIntDict.keys():
            IntDicts[count][key]
            IntDicts[count][key].append(NEWIntDict[key][0])		
            #print(IntDicts)

print('~~~~~Finished reading all files, beginning to process data.~~~~~')
####
#The data analysis bit
#####

#print(EDicts[1])

Ntot = collatedata(NDicts)
Etot = collatedata(EDicts)
Inttot={}
simsteps = list(IntDicts.keys())
inttypes = list(IntDicts[1].keys())
#print(Ntot)
    
for u in IntDicts.keys():
    Inttot[u]={}
    for v in IntDicts[u].keys():
        #print(IntDicts[u][v])
        Inttot[u][v]=[np.mean(IntDicts[u][v]), np.std(IntDicts[u][v])]



####
#The data outputting bit
####			
print('~~~~~Finished processing all data, beginning to write files.~~~~~')

directorymaker(outputpath)
with open('{2}{0}.{1}.energytraj.csv'.format(species,framework, outputpath), 'w') as f:
    f.write(
'''##Raw energy vs iteration data output data from simulation {0} on {1}\n'''.format(species, framework))
    f.write('''N iterations, Etot (kJ/mol)\n''')

    for count,o in enumerate(EDicts.keys()):
        f.write('Isotherm step {0}\n'.format(count+1))
        f.write('n, ')
        f.write(', '.join([x.split('.')[0] for x in filelist]))
        f.write(', Average, Standard deviation')
        f.write('\n')
        for p in EDicts[o].keys():
            f.write('{0}, '.format(p))
            f.write(', '.join([str(x) for x in EDicts[o][p]]))
            ave=round(np.mean(EDicts[o][p]), 3)
            sdv=round(np.std(EDicts[o][p]), 3)
            f.write(', {0}, {1}'.format(ave, sdv))
            f.write('\n')
        f.write('Average stationary energy')
        for q in Etot[o]['Ave']:
            f.write(', {0}'.format(str(round(q,3))))
        f.write('\n')
        f.write('Standard deviation')
        for q in Etot[o]['Std']:
            f.write(', {0}'.format(str(round(q,3))))
        f.write('\n')

with open('{2}{0}.{1}.occupancytraj.csv'.format(species,framework, outputpath), 'w') as f:
    f.write(
'''##N vs iteration data output data from simulation {0} on {1}\n'''.format(species, framework))
    f.write('''Iteration, Ntot (mol/uc):\n''')
    for count,o in enumerate(NDicts.keys()):
        f.write('Isotherm step {0}\n'.format(count+1))
        f.write('n, ')
        f.write(', '.join([x.split('.')[0] for x in filelist]))
        f.write(', Average, Standard deviation')
        f.write('\n')
        for p in NDicts[o].keys():
            f.write('{0}, '.format(p))
            f.write(', '.join([str(x) for x in NDicts[o][p]]))
            ave=round(np.mean(NDicts[o][p]), 3)
            sdv=round(np.std(NDicts[o][p]), 3)
            f.write(', {0}, {1}'.format(ave, sdv))
            f.write('\n')
        f.write('Average stationary occupancy')
        for q in Ntot[o]['Ave']:
            f.write(', {0}'.format(str(round(q,3))))
        f.write('\n')
        f.write('Standard deviation')
        for q in Ntot[o]['Std']:
            f.write(', {0}'.format(str(round(q,3))))
        f.write('\n')#    for o in NDicts.keys():


with open('{2}{0}.{1}.Interactions.csv'.format(species,framework, outputpath), 'w') as f:
    f.write(
'''##interaction strength vs step output data from simulation {0} on {1}\n'''.format(species, framework))
    f.write('''Iteration type, strenght (kJ/mol)\n''')
    for count, o in enumerate(IntDicts.keys()):
        f.write('Isotherm step {0}\n'.format(count+1))
        f.write('Interaction type, ')
        f.write(', '.join([x.split('.')[0] for x in filelist]))
        f.write(', Average, Standard deviation')
        f.write('\n')
        for p in IntDicts[o].keys():
            f.write('{0}, '.format(p))
            f.write(', '.join([str(x) for x in IntDicts[o][p]]))
            ave=round(np.mean(IntDicts[o][p]), 5)
            sdv=round(np.std(IntDicts[o][p]), 5)
            f.write(', {0}, {1}'.format(ave, sdv))
            f.write('\n')
        f.write('\n')


with open('{2}{0}.{1}.Alldata.csv'.format(species,framework, outputpath), 'w') as f:
    f.write('''## Simulation step vs: occupancy, energy, interaction energies (times 4). All data is from a simulation of {0} on {1}\n'''.format(species, framework))
    f.write('''Simulation step, Ntot average(mol/uc), stdev, Etot (kJ/mol), stdev, ''')
    f.write(', stdev, '.join([str(x) for x in IntDicts[1].keys()]))
    f.write('\n')
    for o in Ntot.keys():
        f.write('{0}, '.format(o))
        f.write('{0}, {1}, '.format(np.mean(Ntot[o]['Ave']), np.mean(Ntot[o]['Std'])))
        f.write('{0}, {1}, '.format(np.mean(Etot[o]['Ave']), np.mean(Etot[o]['Std'])))
        for p in Inttot[o].keys():
            f.write(', '.join([str(round(x,5)) for x in Inttot[o][p]]))
            f.write(', ')
        f.write('\n')
		
print('~~~~~Finished writing files, your data should be in {0}~~~~~'.format(outputpath))
