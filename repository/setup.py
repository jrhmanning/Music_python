
import os
import re
import random
import numpy as np
from math import sqrt, log
from pathlib import Path
import logging
import datetime
import argparse

import Antoine.py
import Forcefield.py
now = datetime.datetime.now()

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler('{0}/umbrella.log'.format(xptpath))
handler.setLevel(logging.DEBUG)
handler2 = logging.StreamHandler()
handler2.setLevel(logging.INFO)
logger.addHandler(handler)
logger.addHandler(handler2)
logger.debug('-----------------------------------')
logger.debug('Generating new files on {0}-{1}-{2} at {3}:{4}'.format(now.year, now.month, now.day, now.hour, now.minute))
logger.debug('-----------------------------------')



#############
#Functions used here:

#Reads a pdb file to get atom types, positions, and conneciton information. Does it by finding the phrases 'HETATM' and 'CONECT', and required the lines ot have a specific number of strings in them
def dataextract(species, sourcedir = './'):
    atoms = []
    connections = []
    f = open('{0}{1}.pdb' .format(sourcedir, species), "r")
    for curline in f:
        fish = curline.split()
        #print(fish)
        if fish[0] == "HETATM": #I expect HETATM lines to look like this: HETATM index symbol residue ? x y z ? q symbol
            del(fish[8:10]) #deletes "? q " from the end fo the above
            del(fish[2:5]) #deletes "symbol residue ?" from above
            del(fish[0]) #deletes 'HETATM' from above
            atoms.append(fish)
        elif fish[0] == "CONECT": #I expect CONECT lines to be 'CONECT index1 index2 index3... where index1 is an atom and intex2 etc. are atoms bound to it
            del(fish[0]) #deletes 'CONECT from the above
            connections.append(fish)
    for thing in atoms:
        if thing[4] in element_names:
            thing[4] = element_names[thing[4]] #replaces the atom symbol from above with the element name
    return atoms, connections

#This sorts the connections detected from apdb into a list for a molfile
def connectiontypeswrite(atoms, connections):
    atomnames = []
    pairset = set()
    bondlengths = {}
    for i in atoms:
        atomnames.append(i[4]) #creates a list of atom names from the stom position information from dataxtract
    logger.debug(atomnames)
    #atomstypes = dict.fromkeys(atomnames)
    indices = []
    for i in connections:
        for j in i:
            if j != i[0]:
                pairset.add(tuple(sorted([atomnames[int(i[0])-1], atomnames[int(j)-1]]))) #creates a set of tuples of bonded atom pairs, so there;s no duplicates
                for thing in pairset:
                    bondlengths[thing] = [] #creates a dictionary of the above tuples, so you can load in your bond lengths
    for i in connections:
        for j in i:
            key = tuple(sorted([atomnames[int(i[0])-1], atomnames[int(j)-1]])) #detects the atom bonding pair from each 'i j' pair in the connections list
            if j != i[0]:
                dist = 0 #bond length. This line resets it for each iteration fo the loop
                dist = sqrt((float(atoms[int(i[0])-1][1])-float(atoms[int(j)-1][1]))**2+(float(atoms[int(i[0])-1][2])-float(atoms[int(j)-1][2]))**2+(float(atoms[int(i[0])-1][3])-float(atoms[int(j)-1][3]))**2) #measures bond length in Angstrom
                dist = round(dist, 3)   #rounds to nearest 0.001 of an angstrom. Is this good enough?
                bondlengths[key].append(dist) #appends it tot he list of bonds for atom pairs 'i j' 
    for pair in bondlengths:
        if len(bondlengths[pair]) > 1: #if there is recorded bond length information
            for i in bondlengths[pair]:
                if float(bondlengths[pair][0])-float(i) > 0.001: #if your bond lengths vary more than this, throw a fit.
                    logger.warning('ERROR! Tolerance exceeded for bonds between atoms of types %s.' % (str(pair)))
                    sys.exit()
            bondlengths[pair] = round(sum(bondlengths[pair])/len(bondlengths[pair]), 3) 
    return bondlengths #returns the bondlength dictionary of 'i j pair': length

#This writes a molecule file from the informaiton from dataextract and connectiontypes. It can overwrite your existing files though, beware!
def molwrite(species, atoms, connectiontypes, connections, sourcedir = './', forcefield = None):
    f = open("{0}{1}.mol" .format(sourcedir, species), 'w')
    f.write("### Music molecule construction information, generated from {0}{1}.pdb by the molfile generator \n\n" .format(sourcedir, species)) #Explains that this is autogenerated
    f.write("#Basic molecule information \n") #Required for music?
    f.write("Molecule_Name: {0} \n" .format(species)) #definitely required for music. And it hates tabs in your lines.
    f.write("Coord_Info: listed cartesian rigid\n") #here i'm telling the moelcule file your moecule has listed atoms, cartesian coordinates, and rigid relationships between them
    f.write("{0} \n" .format(len(atoms))) #this prints the number of atoms to your file
    if len(atoms) >1:
        for i in atoms:
            logger.debug(i)
            for j in i:
                f.write('{0} ' .format(str(j))) #prints the atom index, x, y, z, name
            f.write("{0} 0 0\n" .format('0' if forcefield == None else str(Forcefield.IntParams['{0}_{1}'.format(i[-1].split('_')[0], forcefield)][2]))) #prints q, ?, ?. change the '0' in format if you want to set a q to a variable
    else:
        f.write('1 0 0 0 {0} ' .format(str(atoms[-1]))) #prints the atom index, x, y, z, name
        f.write("{0} 0 0\n" .format('0' if forcefield == None else str(Forcefield.IntParams['{0}_{1}'.format(atoms[-1].split('_')[0], forcefield)][2]))) #prints q, ?, ?. change the '0' in format if you want to set a q to a variable
    if connectiontypes:
        f.write("Connect_Info: listed\n") #declares it's listing the number of bond dypes as well as the exact bond pairs
        f.write("{0} # number of unique bond types in your molecule\n" .format(len(connectiontypes)))
        for i in connectiontypes:
            x = str(i) #makes a string of you 'i j' pair
            x = x.strip("('" "')") #cuts out the quites from the string
            x = x.split("', '") #and the comma from between them
            logger.debug(x) #now it should just be 2 atom type names
            f.write("{0} {1} {2}\n" .format(x[0], x[1], connectiontypes[i])) #now it writes atom type 1, atom type 2, length
        for i in connections: #This loops over your CONECT reading list from dataextract to write the bond list
            for j in i:
                f.write("{0} " .format(str(j)))
            f.write('\n')
    f.write("# Degrees of freedom (optional)\n") #writes out the degrees of freedom, essentialy saying you have 3 for monoatomic, and 6 otherwise
    if len(atoms) == 1:
        dof = 3
    else:
        dof = 6
    f.write("Molecule_DOF: {0} #3 for monoatomic, 5 for rigid linear and 6 for  rigid nonlinear, if flexible then generally 3*number of atoms " .format(dof))        
    f.close()
    logger.debug("Molecule file {0}{1}.mol written!" .format(sourcedir, species))

#This function checks if you already have a directory with a certain name, and makes it if not.
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
    logger.debug("Directory {0} written!".format(dxout))

#This writes your atom-atom forcefield information based on the above dictionary at the start of the document
def AtmAtmMover(elements, forcefield, dxout = "./mapgen/test/"):
    tot_el_list = []
    tot_el_list.append(MOF_el_list) #takes the global list from the preamble to getyour mof molecule atom types
    tot_el_list.append(elements) #reads what elements are in yoru sorbent
    hicut = 18 #sets your hicut in A
    with open("{0}atom_atom_file" .format(dxout), "w") as file:
        file.write("""#This is an autogenereated interactions file for music, based on the python script Umbrella-v2. It's probably best to look this over before running your sims\n#Lennard-Jones interactions\n""")
        for i in range (0, len(MOF_el_list)):
            file.write("""{0} {0} LJ SIG@{1} EPS@{2} HICUT@{3}\n""".format(MOF_el_list[i], Forcefield.IntParams[MOF_el_list[i]][0], Forcefield.IntParams[MOF_el_list[i]][1], hicut)) #LJ self-parameters for the MOF-MOF interactions
            for j in range (i+1, len(MOF_el_list)):
                file.write("""{0} {1} LJ OFF\n""".format(MOF_el_list[i], MOF_el_list[j])) #LJ parameters for MOF 'i j' pairs, which are generally off
            for k in elements:
                file.write("""{0} {1} LJ SIG@{3} EPS@{4} HICUT@{2}\n""".format(MOF_el_list[i], k, hicut, str(round((Forcefield.IntParams[MOF_el_list[i]][0]+Forcefield.IntParams['{0}_{1}'.format(k.split('_')[0], forcefield)][0])/2, 3)), str(round(sqrt(Forcefield.IntParams[MOF_el_list[i]][1]*Forcefield.IntParams['{0}_{1}'.format(k.split('_')[0], forcefield)][1]), 3)))) #MOF-fluid forcefield paramters
            file.write("\n")
        for i in range (0, len(elements)):
            raw_el = str(elements[i].split('_')[0] + '_' + str(forcefield)) #This line says which interaction types your fluid will take from the dictionary at the start of the script
            file.write("""{0} {0} LJ SIG@{1} EPS@{2} HICUT@{3}\n""".format(elements[i], Forcefield.IntParams[raw_el][0], Forcefield.IntParams[raw_el][1], hicut)) #LJ self-parameters for the fluid-fluid interactions
            for j in range (i+1, len(elements)):
                logger.debug('i = '+ str(i)+ ', j= '+str(j) )
                file.write("""{0} {1} LJ SIG@{3} EPS@{4} HICUT@{2}\n""".format(elements[i], elements[j], hicut, str(round((Forcefield.IntParams['{0}_{1}'.format(elements[i].split('_')[0], forcefield)][0]+Forcefield.IntParams['{0}_{1}'.format(elements[j].split('_')[0], forcefield)][0])/2, 3)), str(round(sqrt(Forcefield.IntParams['{0}_{1}'.format(elements[i].split('_')[0], forcefield)][1]*Forcefield.IntParams['{0}_{1}'.format(elements[j].split('_')[0], forcefield)][1]),3)))) #LJ parameters for the fluid-fluid 'i j' interactions
            file.write("\n")
        file.write("#Coulombic region\n")
        for i in range (0, len(MOF_el_list)):
            file.write("""{0} {0} COUL OFF\n""".format(MOF_el_list[i])) #MOF-MOF self-coulombic interactions (off)
            for j in range (i+1, len(MOF_el_list)):
                file.write("""{0} {1} COUL OFF\n""".format(MOF_el_list[i], MOF_el_list[j])) #MOF-MOF interactiosn (also off)
            for k in elements:
                file.write("""{0} {1} COUL OFF\n""".format(MOF_el_list[i], k)) #MOF-fluid pairwise interactions (definitely off, use a map for this)
            file.write("\n")
        for i in range (0, len(elements)):
            file.write("""{0} {0} SUM FAST EWALD HICUT@{1}\n""".format(elements[i], hicut)) #fluid fluid self-interactions - use wolf cpoulombic. partial charge values are stored in the molecule files - see molfilewriter
            for j in range (i+1, len(elements)):
                file.write("""{0} {1} SUM FAST EWALD HICUT@{2}\n""".format(elements[i], elements[j], hicut)) #fluid fluid 'i j' interactions - use wolf cpoulombic. partial charge values are stored in the molecule files - see molfilewriter
            file.write("\n")
    logger.debug("Atom-Atom interaction file written!")

#This writes your sorb-sorb file telling music which intermolecular interactions are on
def SorbSorbWriter(species, elements,framework, dxout = "./mapgen/", pmap = True, emap = False):
    filename = "sorb_sorb_file"
    with open("{0}{1}" .format(dxout, filename), 'w') as file:
        file.write("{0} {0} NCOUL OFF\n{0} {0} COUL OFF\n\n" .format(framework)) #MOF-MOF interactions are off
        if isinstance(species, list) == True: #This checks if your species is a list or a string, which makes it work for multiple species at once
            if pmap == True: #using a MOF-fluid LJ potential map if true
                for i in species:            
                    file.write("{0} {1} NCOUL MAP@{1} FAST ".format(i, framework)) #declares it'll use a map
                    for j in elements:
                        file.write("{0}@PMAP@{1}.{0}.p " .format(j, framework)) #writes your map names. WARNING - it's a 18 character limit on map names, and no '-' allowed!
                    file.write("\n")
            else:
                for i in species:
                    file.write("{0} {1} NCOUL BASIC LJ FAST\n".format(i, framework)) #pairwise fluid-framework interactions
            if emap == True: #using a MOF-fluid coulomb map if true
                for i in species:
                    file.write("{0} {1} COUL MAP@{1} FAST ".format(i, framework)) #declares it'll use a map
                    for j in elements:
                        file.write("{0}@EMAP@{1}.{2}.e " .format(j, framework, 'Probe')) #writes your map names. WARNING - it's a 18 character limit on map names, and no '-' allowed!
                    file.write("\n")
            else:
                for i in species:
                    file.write("\n{0} {1} COUL OFF\n" .format(i, framework)) #pairwise coulomb interactionsdon't work. Trust me, this is better.
                logger.warning("Framework-fluid coulombic interactions are off") #warns you your coulomb fluid-framework stuff isn't happeneing
            file.write("\n{0} {0} NCOUL BASIC LJ FAST\n{0} {0} COUL SUM FAST EWALD KMAX@15 KAPPA@6.7\n\n" .format(i)) #I think this line shouldn't be here. Needs testing!
        elif isinstance(species, str) == True: 
            if pmap == True: #using a MOF-fluid LJ potential map if true
                file.write("{0} {1} NCOUL MAP@{1} FAST ".format(species, framework)) #declares it'll use a map
                for j in elements:
                    file.write("{0}@PMAP@{1}.{0}.p " .format(j, framework)) #writes your map names. WARNING - it's a 18 character limit on map names, and no '-' allowed!
                file.write("\n\n")
            else:
                file.write("{0} {1} NCOUL BASIC LJ FAST".format(species, framework)) #pairwise fluid-framework interactions
                file.write("\n\n")
            if emap == True: #using a MOF-fluid coulomb map if true
                file.write("{0} {1} COUL MAP@{1} FAST ".format(species, framework)) #declares it'll use a map
                for j in elements:
                    file.write("{0}@EMAP@{1}.{2}.e " .format(j, framework, 'Probe')) #writes your map names. WARNING - it's a 18 character limit on map names, and no '-' allowed!
                file.write("\n\n")
            else:
                file.write("{0} {1} COUL OFF\n" .format(species, framework)) #pairwise coulomb interactionsdon't work. Trust me, this is better.
                logger.warning("Framework-fluid coulombic interactions are off") #warns you your coulomb fluid-framework stuff isn't happeneing
        if isinstance(species, str) == True:
            file.write("\n{0} {0} NCOUL BASIC LJ FAST\n{0} {0} COUL SUM FAST EWALD KMAX@15 KAPPA@6.7\n\n" .format(species))  #fluid fluid interactions, with Wolf coulombic interactions 
        elif isinstance(species, list) == True: 
            for i in species:
                file.write("\n{0} {0} NCOUL BASIC LJ FAST\n{0} {0} COUL SUM FAST EWALD KMAX@15 KAPPA@6.7\n\n" .format(i))  #fluid fluid interactions, with Wolf coulombic interactions 
    logger.debug("Sorb-Sorb file written!")    

#Writes your intramolecular file
def IntraWriter(species, elements, framework, dxout = "./mapgen/"):
    with open("{0}intramolecular_file" .format(dxout), "w") as file:
        file.write("####This intramolecular interaction file is workable for making a map of {0} in a 2x2x2 IRMOF-1 framework\n" .format(species)) #outdated preamble
        if isinstance(species, list) == True: #Lets you handle multiple fluids at once
            for i in species:
                file.write("Intra: {0}\n" .format(i)) #no intramolecular interactions at all
        elif isinstance(species, str) == True: 
            file.write("Intra: {0}\n" .format(species))
        for i in list(elements): #lets you create a map of every element in your system. Having the extra lines in your production script don't matter
            file.write("Intra: {0}\n".format(i))
        file.write("Intra: {0}" .format(framework)) #No framework intramolecular interactions happening
    logger.debug("Intra file written!")

#bundles the above 3 functions into one for simplicity
def Intsetup(mol_name, elements, framework, forcefield, dxout, pmap = True, emap = False):
    AtmAtmMover(elements, forcefield, dxout)
    SorbSorbWriter(mol_name, elements, framework, dxout, pmap, emap)
    IntraWriter(mol_name, elements, framework, dxout)

#Writes a bashscript for making your maps, compatible with SLURM
def MapMakeRunWriter(species, elements, parentdir, dxout = "./mapgen/"):
    with open("{0}run.mapmaker" .format(dxout), "w") as file:
        file.write("""#!/usr/bin/env bash 
#SBATCH --job-name={0}.map
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jrhm21@bath.ac.uk #Change this to your email!
#SBATCH --account=free

# -- Set up the environment
module purge
module load group ce-molsim stack
module load music/std

# -- Find my directory
cd {2}{1}\n #Change this bit to get to your directory
export ATOMSDIR=../../atoms
export MOLSDIR=../../molecules
export PMAPDIR=../../maps/{0}
export EMAPDIR=../../maps
\n\n# --Run\n""" .format(species, dxout.split(".")[-1], parentdir))
        for count, i in enumerate(elements, 1):
            file.write("music_mapmaker makemap_{0}.ctr > logfile{1}\n" .format(i, count, parentdir))
    os.chmod("{0}run.mapmaker" .format(dxout), 0o777) #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    logger.debug("Runfile written!")   

#Writes a control file for making your maps
def mapctrlfilewriter(species, framework, dxout = "./mapgen/"):
    x = random.randint(0,99999) #makes the random seed 
    with open("{0}makemap_{1}.ctr" .format(dxout, species), "w") as file:
        file.write("""------ General Information ------------------------------------------
{0}_pmap map probe in {2}
1               # No. of iterations
1               # No. of steps between writes to output/log file
1               # No. of steps between writes to crash file
1               # No. of steps between writes to config. file
2                    # Start numbering simulations from .
{1}              # Iseeed
4                    # specifies contents of config file,
{0}.res           # Restart File to write to
{0}.con           # Configuration File
------ Atomic Types --------------------------------------------------
5                                  # number of atomic types

{0}
{0}.atm
           
Carbon                            # atom type
Carbon.atm                        # basic atom info file

Oxygen                             # atom type
Oxygen.atm                        # basic atom info file

Hydrogen                             # atom type
Hydrogen.atm                        # basic atom info file

Zinc                         # atom type
Zinc.atm                        # basic atom info file
------ Molecule Types -------------------------------------------------
2                       # number of sorbate types

{0}                   # sorbate 
{0}.mol               # sorbate coordinates file

{2}                    # sorbate 
{2}.mol                # sorbate coordinates file
------ Simulation Cell Information --------------------------------------
{2}                   # Fundamental cell type
2, 2, 2                 # No. of unit cells in x, y, z direction
1, 1, 1                 # (1 = Periodic) in x, y, z
------ Forcefield Information -------------------------------------------
BASIC
SPC
atom_atom_file     # atom-atom interaction file 
sorb_sorb_file     # sorbate-sorbate interaction file (optional)
intramolecular_file  # intramolecular interactions
------ Mapmaker Information --------------------------------------------
1              # Number of maps to make

{2}           # Sorbate to map
{0}       # Sorbate to probe map with
NCOUL LJ       # Interaction type to map
0.2            # Approximate grid spacing (Ang)
100.0          # High end potential cutoff (kJ/mol)
{2}.{0}.p           # Map filename or AUTO
------ Configuration Initialization -------------------------------------
{0}                            # Sorbate_Type  
Molecule NULL                              # Source Filename
{2}                            # Sorbate_Type
Fixed NULL                     # Source Filename""" .format(species, x, framework))
    logger.debug("Map control file written!")

#calculates your saturation pressure using the antoine equation parameters from the preamblem, so you can autogenerate pressures
def pSat(Species, T):#This function checks the Species is there and that you're in the right temperature range. 
    Chemical = Species.split("_")[0]
    logger.debug(Chemical)
    Press = 1 #pressure defaults to 1 kPa
    if Chemical in Antoine.Antoine:                              #It then reads out your saturation pressure in kPa for the user benefit and to get the rest of the program working
        if Antoine.Antoine[Chemical][3] <= T <= Antoine.Antoine[Chemical][4]:
            logger.info("I have Antoine parameters for {0} at {1} K." .format(Chemical, T))
        elif T < Antoine.Antoine[Chemical][3]:
            logger.warning("WARNING, that temperature is too low for my Antoine parameters for {0} (my limit is {1} K). I'll persevere anyway." .format(Chemical, Antoine.Antoine[Chemical][3]))
        elif Antoine.Antoine[Chemical][4] < T:
            logger.warning("WARNING, that temperature is too high for my Antoine parameters for {0} (my limit is {1} K). I'll persevere anyway." .format(Chemical, Antoine.Antoine[Chemical][4]))
        else:
            logger.warning("Something weird happened, I've probably got a bug. Oops!")
        Press = 100*(10**(Antoine.Antoine[Chemical][0]-(Antoine.Antoine[Chemical][1]/(Antoine.Antoine[Chemical][2]+T))))
    else:
        logger.warning("I'm sorry, I don't have that Species in my database.")
    return Press

#Creates a log-linear pressure isotherm using the saturation pressure from pSat and preamble defined number of points/predefined points 
def isothermcalculator(Pressure, n, minrelpress, maxrelpress):                                  #This function produces the isotherm as a data list, using pSat calculated by function pSat and the user defined isotherm length n
    m = n-len(FixedPoints)-1                                         #Pressure points re linearly distributed above 0.09 kPa to pSat, which may eb a bad idea. Who knows?
    AllPoints = [] #your list of total pressures
#    print(FixedPoints[-1])
#    print(Pressure-FixedPoints[-1])
    for point in range(1,n+1):
        relpress = minrelpress+(point*(maxrelpress-minrelpress)/n) #Calculates your partial pressure list between minrelpress and 0. Dividing the (point/float(n+1)) statement by a number lowers your max pressure to a fraction of the sat pressure
        press = Pressure*10**relpress #converts the above into absolute pressures
        press = round(press, 5) #Rounds the float to 3 decimal points, for simplicity
        AllPoints.append(press)
        AllPoints.sort() #sorts your points into ascending order
    isotherm = AllPoints
    logger.info("PSat = {0}KPa, Isotherm = {1}".format(Pressure, isotherm))
#    print(isotherm[2])
#    print(len(isotherm))
    return isotherm
    
#this function write s a .dat file based on the pressures you've calculated in isothermcalculator
def PressureFileWriter(Species, T, Pressure, isotherm, dirout):                 # This function writes your file to a specific directory, but currently doesn't say if your temp is out of range. Annoying!!#
#    print(isotherm)
    with open("{0}pressure.{1}.{2}.dat" .format(dirout, Species, T), "w") as file:
        file.write("{0!s} #PSat at {1!s} K is {2:1.3f} kPa. \n" .format(Species, T, Pressure))
        file.write(str(len(isotherm)) + "\n")
        file.write(', '.join(str(thing) for thing in isotherm))

#Writes a bashscript for taskfarming yuor simulations, compatible with SLURM
def TaskfarmRunWriter(Species, T, framework, parentdir, dirout):
    with open("{0}run.taskfarmer" .format(dirout), 'w') as file:
        file.write("""#!/usr/bin/env bash
#SBATCH --job-name={0}.{2}.{1}k
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --time=12:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jrhm21@bath.ac.uk #Change this to your email!
#SBATCH --account=re-ce1100

# -- Set up the environment
module purge
module load group ce-molsim stack
module load taskfarmer

#locate the experiment directory
cd {4} #Change this bit to get to your directory
cd {3}

# -- Run
mpirun -np 16 taskfarmer -f taskfarm

# -- python analyse the results
python ./isothermextractor.py""".format(Species, T, framework, dirout, parentdir)) #uses my isotherm extractor script to pull all your isotherms out and put them into .csv files in the directory abocve the taskfarm run file
    os.chmod("{0}run.taskfarmer" .format(dirout), 0o777)   #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    with open("{0}taskfarm" .format(dirout), 'w') as file1: #Writes the commands to taskfarm into a separate taskfram file
        for i in range(1,17):
            file1.write("bash ./{0:02d}/run.gcmc\n".format(i))
    with open("{0}taskfarm.backup" .format(dirout), 'w') as file2: #writes a backup you can restore the above to if there's a bug
        for i in range(1,17):
            file2.write("bash ./{0:02d}/run.gcmc\n".format(i))
    logger.debug("Taskfarm things file written!")

#Writes a python script to extract data from your experiments after taskfarming
def IsothermExtractMover(Species, T, framework, dirout, n = 20):
    with open("{0}isothermextractor.py" .format(dirout), 'w') as file:
        file.write("""from os import listdir
import numpy as np

isotherms = np.zeros((17, {3})) # this is a matrix of a 16x {3}-point isotherms with identical pressure values

def isothermextract(page): #this function scrapes subdirectories of name for isotherm data for a given isotherm file name
    z = 0
    z = int(page) #sets z to be your taskfarm run
#    print(z)
    drx = str("./{{0:02d}}/isotherm.{0}" .format(page)) # defines the isotherm
#    print(drx)
    f = open(drx, "r")
    x = []
    page = [] #redefines your input variable, but somehow still works
    for curline in f:
        if "#" not in curline:
            fish = curline.split()
            if len(fish) == 2:
                x.append(fish[0])
                page.append(fish[1])
    f.close()
    isotherms[0] = x #overwrites the first column of your np.zeros matrix to be your pressure oiubts
    isotherms[z] = page #overwrites column z in your np.zeros matrix with resutls from isotherm z
    return isotherms
#####################################################################################################################################################	
for i in range (1,17): #loops over your 16 isotherms
    SuccessCount = 0 #this variable preents issues downstream if some of your simulations fail
    filenames = listdir("./{{0:02d}}/" .format(i))
    for file in filenames:
        if file == "isotherm.{0}": #change the np.zeroes to be your results
            isothermextract(i)
            SuccessCount += 1
    if SuccessCount != 1: #Checks the abovle was successful
        print("WARNING! Something went wrong in simulation {{0:02d}}" .format(i)) #complains if not

x = isotherms[0] #This is now getting a mean and stdev of your results. Still works even if some of your simulations don't!
avg=[]
stdev=[]

avgisotherms = np.zeros((3, {3})) #now makes a 3x{3} matrix to hold your averaged data and stdevs
for i in range (0, {3}):
    spread = []
    for thing in range (2,17):
        if isotherms.item((thing, i)) != 0: #prevents your failed simulations from messing things up
            spread.append(isotherms.item((thing, i)))
    avg.append(np.mean(spread)) #gets your mean
    stdev.append(np.std(spread)) #gets your stdev
avgisotherms[0] = x
avgisotherms[1] = avg
avgisotherms[2] = stdev
np.savetxt("../{2}.{1}.{0}.results.csv", avgisotherms, delimiter=",") #output it to this file""" .format(Species, T, framework, n))
    logger.debug("Python isotherm extractor written!")

#Writes a .ctr file for your gcmc simulations (production and for getting your final configuarions at the end). This currently cannot cope with multiple sorbents at once, creates a 2x2x2 cell, and has bias insert atom softcoded in
def GcmcControlChanger(species, elements, T, n, framework, dirout, iterations = '750000', Restart = None, name = 'gcmc.ctr', pressure = 'file'):
    logger.debug(elements)
    x = random.randint(0,99999) #sets your random seed
    if isinstance(elements, list): #lets you use multi-element sorbent molecules
        n_species = len(MOF_el_list)+len(elements) #sets the total number of element types in your gcmc
    elif isinstance(elements, int): #lets you use single-element sorbents
        n_species = len(MOF_el_list)+elements
    with open("{0}{1}" .format(dirout, name), 'w') as file:
        file.write("""#This control file was written by the python scrpt umbrellav2-6 You probably ought to check me before use!
------ General Information ------------------------------------------
{0} molecule in {1} 
{4}              # No. of iterations, defaults to 750000
50000                # No. of steps between writes to output/log file
100000                # No. of steps between writes to crash file
2500                  # No. of steps between writes to config. file
{5}                   # Start numbering simulations fromhere. defaults to 1 for production run, 30 for config runs
{2}                #random seed
3                    # specifies contents of config file, outdated?
{1}.{0}.res         # Restart File to write to
{1}.{0}.con          # Configuration File
------ Atomic Types --------------------------------------------------
{3}                    # number of atomic types            \n\n""".format(species, framework, x, n_species, iterations, 1 if Restart == None else 30))
        for i in elements:
            file.write("""{0}   #atom type\n{0}.atm #atom file name\n\n""".format(i))
        for i in MOF_el_list:
            file.write("""{0}   #atom type\n{0}.atm #atom file name\n\n""".format(i))
        file.write("""------ Molecule Types -------------------------------------------------
2                    # number of sorbate types

{0}               # sorbate
{0}.mol           # sorbate coordinates file

{1}                # sorbate
{1}.mol             # sorbate coordinates file
------ Simulation Cell Information ------------------------------------
{1}                # Fundamental cell file
2, 2, 2              # No. of unit cells in x, y, z direction
1, 1, 1              # (1 = Periodic) in x, y, z
------ Forcefield Information -------------------------------------------
BASIC
SPC
atom_atom_file       # atom-atom interaction file
sorb_sorb_file       # sorbate-sorbate interaction file
intramolecular_file  # intramolecular interaction file/specification
------ Ideal Parameters -----------------------------------------------
Ideal                # Equation of State
1                    # no. of sorbates
{0}              # Sorbate Name
------ GCMC Information -----------------------------------------------
1                 # No. of iterations
{2}              # temperature
Ideal Parameters   # Tag for the equation of state (NULL = Ideal Gas)
{3}                  # No. of simulation points
5000                # Block size for statistics
1                  # no. of sorbates
          -------------------------
{0}            # Sorbate Name
{5}           #  pressure
Null               # sitemap filename (Null = no sitemap)
4                  # no of gcmc movetypes
1.0, 1.0, 1.0, 1.0      # move type weights
RINSERT                   # type of move.1
RDELETE                   # type of move.2
RTRANSLATE                # type of move.4
0.2, 1                    # Delta Translate, adjust delta option (0=NO, 1=YES)
RROTATE
0.2, 1
------ Configuration Initialization -------------------------------------
{0}             # Sorbate_Type\n""".format(species, framework, T, n, elements[-1], 'pressure.{0}.{1}.dat'.format(species, T) if pressure == 'file' else pressure)) #lets you define your pressure points manually, defaults to the values generated from the isothermgenerator
        if Restart == None: #for if you're starting a new isotherm
            file.write("""GCMC NULL
{0}              # Sorbate_Type
FIXED NULL
--------  Main Datafile Information --------
Energy, position, pair_energy  # contents of datafile""".format(framework))
        else: #lets you get configs from a main isotherm res file
            file.write("""{0}
{1}              # Sorbate_Type
FIXED NULL
--------  Main Datafile Information --------
Energy, position, pair_energy  # contents of datafile""".format(Restart, framework))
    logger.debug("GCMC control file written!")

#Writes a .ctr file for your postprocessing. There;s orobably a better way than using music_post, but I'm not there yet
def PostControlChanger(Species, n, framework, dirout, prefix, cutoff = '0'):
    with open("{0}/{1}.post.ctr" .format(dirout, 'full' if cutoff == '0' else 'trunc'), 'w') as file:
        file.write("""####This section is apparently required for working with any post code. Who knows why? not me!
#
#
------------------------------------------------------------
   ### Required section ######
-- Post Processor Information ------------
GCMC                            # Type of simulation GCMC, NVTMC , MD ....
./{1}.{0}.con                    # basename for config files
1, {2}                          # first and last file numbers
{5}.post.ctr.out                       # name for new ctrlfile that will regenerated
{3:02d}.{5}.postfile          # Base name for output files
{4}, 0                         # Percentages of data to skipped at start and end 


# The sections below are necessary only if you want the corresponding 
# analysis performed
# ---------------- ALL OF THEM ARE OPTIONAL ------------------------


####    This section is reqd for energy averages in your post code output files
####    as of now only total enrgies vs sim. step
------ Post : Energy Average Info -----------------------------------
100       # Number of blocks into which data should be divided for stats

####    This section is reqd for Loading averages in your post code outputfiles
####    as of now only species loading vs sim. step (for all species)
------ Post : Loading Average Info -----------------------------------
100       # Number of blocks into which data should be divided for stats""".format(Species, framework, n, prefix, cutoff, 'full' if cutoff == '0' else 'trunc'))
    logger.debug("Post control file written!")

#Writes a bashscript for actually running your simulations, compatible with SLURM
def GcmcRunWriter(Species, T, framework, parentdir, dirout, isotherm, iteration):
    with open("{0}run.gcmc".format(dirout), 'w') as file:
        file.write("""#!/usr/bin/env bash
#SBATCH --job-name={0}.{1}
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jrhm21@bath.ac.uk #Change this to your email!
#SBATCH --account=free

# -- Set up the environment
module purge
module load group ce-molsim stack
module load music/std

#locate the experiment directory
cd {3} #Change this bit to get to your directory
cd {2}

#create symbolic links to your interactions files and pressure file
ln -s ../../../../atoms atoms
ln -s ../../../../molecules molecules
ln -s ../../../../maps/{0} maps

ln -s ../atom_atom_file atom_atom_file
ln -s ../intramolecular_file intramolecular_file
ln -s ../sorb_sorb_file sorb_sorb_file
ln -s ../pressure.{0}.{1}.dat pressure.{0}.{1}.dat

#set the paths for atoms, molecules, pmaps, and emaps WARNING softcoded in umbrella, check it's right!
export ATOMSDIR=atoms
export MOLSDIR=molecules
export PMAPDIR=maps
export EMAPDIR=maps


# -- Run
music_gcmc gcmc.ctr > logfile.gcmc #runs your main simulation
music_post full.post.ctr > logfile.post.full #runs your postprocessing
cp {4:02d}.full.postoutput ../fullpostfiles/
music_post trunc.post.ctr > logfile.post.trunc #runs your postprocessing
cp {4:02d}.trunc.postoutput ../truncpostfiles/


""".format(Species, T, dirout, parentdir, iteration))
        for i, value in enumerate(isotherm): #for each isotherm point you're simulating
            file.write('music_gcmc {0}kpa_restart.ctr > {1}_restart.logfile\nmv finalconfig.xyz {3:02d}.{2}.{0}kpa.xyz\n'.format(value, int(i)+1, framework, int(iteration))) #set up a simulation to generate a new xyz file and move it to be names after the pressure point 
    os.chmod("{0}run.gcmc" .format(dirout), 0o777) #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    logger.debug("Run file written!")
