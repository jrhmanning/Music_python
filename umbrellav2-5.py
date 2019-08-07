#This is a rewrite of the earlier umbrella program for running music simulations given input parameters.
#The plan is to make this python script a one-stop shop for setting up/wunnign a pyython simulation, and leave the actual functions elsewhere.
#So this script will just take in all required inputs, and print template files for map generation (if required) and simulations.
#The workflow will go like this: 
#1. Set global variables for things to be solftcoded in
#2. Take in user input variables, (either 1 by 1 or as a list?)
#3. If the user inputs a pdb file, create a molfile from this
#4. write the appropraite pair interaction files (atom_atom, sorb_sorb, intra)
#5. Make an isotherm at your temperatuure of x points, equally spaced for a semilog graph
#6.Write the above out to the appropriate directories
#7. Report on what it's done.
##########
#0. preamble
import os
import re
import random
import numpy as np
from math import sqrt, log
#try:
#    input = raw_input("What is your sorbent species called?")
#except NameError:
#    pass
element_names = {'H': 'Hydrogen',
                 'C': 'Carbon',
                 'Cl': 'Chlorine',
                 'O': 'Oxygen',
                 'N': 'Nitrogen',
                 'F': 'Fluorine'
}
Antoine = {                                                 #This is your library of Antoine equatino parameters. The key is the name of the Species that will be written
    "CCL4":(4.02291, 1221.781, -45.739, 193, 350),          #Indices 0-2 are antoine parameters A, B, and C taken from NIST,
    "Chloroform":(4.20772, 1233.129, -40.953, 215, 334), #Indices 3 and 4 are the Antoine equation validity range, as quoted on the NIST website
    "DCM":(4.53691,1327.016,-20.474, 233, 313),
    "Chloromethane":(4.91858, 1427.529, 45.137, 303, 416),
    "Chloromethane2":(4.22507, 951.561, -23.468, 198, 278), 
    "Methane":(4.22061, 516.689, 11.223, 110, 190),
    "Methanol":(5.20409, 1581.341, -33.5, 288, 356.8),
    "THF":(4.12118, 1202.942, -46.818, 296.29, 372.8),
    "DMF":(3.93068, 1337.716, -82.648, 303, 363)
    }

IntParams = { #Sig (ang), Eps(K), Q (e)
"Zinc":(2.763, 62.34, 0), #From UFF
"Carbon":(3.898, 47.81, 0), #From Dreiding
"Hydrogen":(3.195, 7.642, 0), #Dreiding
"Oxygen":(3.3, 48.11, 0), #Dreiding
"Oxygen_clust":(3.66, 0.069, 0), #UFF
"Oxygen_ligand":(3.3, 48.11, 0), #Dreiding
"Oxygen_DMF":(3.3, 226, -0.5), #DMF parameters taken from 10.1039/C4CP05961A
"Me_DMF":(3.80, 69, 0.28),
"Nitrogen_DMF":(3.20, 144, -0.57),
"Carbon_DMF":(3.70, 47.3, 0.45),
"Hydrogen_DMF":(2.20, 7.18, 0.06),
"Carbon_Chloroform":(3.41, 68.94, -0.235), #Chlorfoorm parametrs taken from 10.1021/jp0535238
"Hydrogen_Chloroform":(2.81, 10.06, 0.355),
"Chlorine_Chloroform":(3.45, 138.58, -0.04)
}

MOF_el_list = ["Carbon", "Hydrogen", "Oxygen", "Zinc"]

FixedPoints = [0.01, 0.03, 0.05, 0.07, 0.09]              #These are the fixed low pressure points I'll include at the start fo your isotherm
#############
##1&2. Global variables:
##Framework - the name of the MOF to be used. Combined with version to get a molfile
framework = 'IRMOF1'
##Species - your sorbent. Only 1 for now
species = "Chloroform"
#species = input("What is your sorbent species called?")
#type(species)
##Interactions version (1. UFF, lj only; 2. UFF, coul + lj; 3. UFF/dreiding lj only; 4.UFF/dreiding lj+coul 5. UFF/dreiding, lj+coul, framework charges (CoRE MOF) )
interactions_version = 5 
#interactions_version = eval(input("What interactions version would you like to simulate?"))
#type(interactions_version)
##Temperature
#T = 298
T = eval(input("What temperature (in K) would you like to test?"))
type(T)
##number of isotherm points
n = 20
#n = eval(input("How many isotherm points do you want?"))
#type(n)
#############
##3. Check for (and make) a molecule file as required

#function to be exported to a repository. Checks if you already have a molecule file for your species
def filecheck(species, version, xt = 'mol', sourcedir = './', framework = 'IRMOF1'):
    filenames = os.listdir('{0}' .format(sourcedir))
    filelist = {"Test": 0}
    name = []
    exists = False
    correct_version = False
    for file in filenames:
        name = re.split("\W+|_", file)
        #print(name)
        if name[-1] == xt:
#            print(name[1][-1])
            if xt == 'mol':
                if len(name) == 3:
                    filelist[name[0]] = name[1][-1]
                elif 1 <= len(name) <= 2:
                    filelist[name[0]] = 0
                elif len(name) == 4:
                    filelist['{0}_{1}'.format(name[0], name[1])] = name[2][-1]
            if xt == 'pdb':
                filelist[name[0]] = 0
#    print(molfilelist)
    if species in filelist:
        print("File name {0} with type {1} found!".format(species, xt))
        exists = True
        print(filelist)
        if isinstance(filelist[species], list) == True:
            if str(filelist[species][-1]) == str(version):
                correct_version = True
        elif isinstance(filelist[species], str) == True:
            if str(filelist[species]) == str(version):
                correct_version = True
        elif isinstance(filelist[species], int) == True:
            if filelist[species] == version:
                correct_version = True
    return exists, correct_version

def mapcheck(species, version, xt = 'pmap', sourcedir = './', framework = 'IRMOF1'):
    filenames = os.listdir('{0}' .format(sourcedir))
    filelist = {"Test": 0}
    name = []
    exists = False
    correct_version = False
    for file in filenames:
        name = file.split(".")
#        print(name, len(name))
        if name[-1] == xt:
#            print(name[1])
#            print(name[1][-1])
            if len(name) == 3:
                sorb_name = name[1].split('_')
                #print(sorb_name)
                filelist['{0}_{1}'.format(sorb_name[0], sorb_name[1])] = [name[0], sorb_name[-1][-1]]
            else:
                print("Warning, pmap life name isn't framework_vX.atom_vX.pmap!")
                filelist[name[1]] = 0
    #print(filelist)
    if species in filelist:
        print("File type {0} {1} found!".format(species, xt))
        exists = True
        #print(filelist)
        if int(filelist[species][-1]) == version:
            correct_version = True
    return exists, correct_version

def dataextract(species, sourcedir = './'):
    atoms = []
    connections = []
    f = open('{0}{1}.pdb' .format(sourcedir, species), "r")
    for curline in f:
        fish = curline.split()
        #print(fish)
        if fish[0] == "HETATM":
            del(fish[8:10])
            del(fish[2:5])
            del(fish[0])
            atoms.append(fish)
        elif fish[0] == "CONECT":
            del(fish[0])
            connections.append(fish)
    for thing in atoms:
        if thing[4] in element_names:
            thing[4] = element_names[thing[4]]
    return atoms, connections

def connectiontypeswrite(atoms, connections):
    atomnames = []
    pairset = set()
    bondlengths = {}
    for i in atoms:
        atomnames.append(i[4])
    print(atomnames)
    #atomstypes = dict.fromkeys(atomnames)
    indices = []
    for i in connections:
        for j in i:
            if j != i[0]:
                pairset.add(tuple(sorted([atomnames[int(i[0])-1], atomnames[int(j)-1]])))
                for thing in pairset:
                    bondlengths[thing] = []
    for i in connections:
        for j in i:
            key = tuple(sorted([atomnames[int(i[0])-1], atomnames[int(j)-1]]))
            if j != i[0]:
                dist = 0
                dist = sqrt((float(atoms[int(i[0])-1][1])-float(atoms[int(j)-1][1]))**2+(float(atoms[int(i[0])-1][2])-float(atoms[int(j)-1][2]))**2+(float(atoms[int(i[0])-1][3])-float(atoms[int(j)-1][3]))**2)
                dist = round(dist, 3)   #rounds to nearest 0.001 of an angstrom. Is this good enough?
                bondlengths[key].append(dist)
    for pair in bondlengths:
        if len(bondlengths[pair]) > 1:
            for i in bondlengths[pair]:
                if float(bondlengths[pair][0])-float(i) > 0.001:
                    print('ERROR! Tolerance exceeded for bonds between atoms of types %s.' % (str(pair)))
                    sys.exit()
            bondlengths[pair] = round(sum(bondlengths[pair])/len(bondlengths[pair]), 3)
    return bondlengths

def molwrite(species, atoms, connectiontypes, connections, interactions_version, sourcedir = './'):
    f = open("{0}{1}_v{2}.mol" .format(sourcedir, species, interactions_version), "w")
    f.write("### Music molecule construction information, generated from {0}{1}.pdb by the molfile generator \n\n" .format(sourcedir, species))
    f.write("#Basic molecule information \n")
    f.write("Molecule_Name: {0}_v{1} \n" .format(species, interactions_version))
    f.write("#Atom coordinate \nCoord_Info: listed cartesian rigid\n") #here i'm telling the moelcule file your moecule has listed atoms, cartesian coordinates, and rigid relationships between them
    f.write("{0} \n" .format(len(atoms))) #this prints the number of atoms to your file
    if len(atoms) >1:
        for i in atoms:
            for j in i:
                f.write('{0} ' .format(str(j))) #prints the atoms
            f.write("{0} 0 0\n" .format('0'))
    else:
        f.write('1 0 0 0 {0} ' .format(str(atoms[-1])))
        f.write("{0} 0 0\n" .format('0'))
    if connectiontypes:
        f.write("Connect_Info: listed\n") #declares it's listing the number of bond dypes as well as the exact bond pairs
        f.write("{0} # number of unique bond types in your molecule\n" .format(len(connectiontypes)))
        for i in connectiontypes:
            x = str(i)
            x = x.strip("('" "')")
            x = x.split("', '")
            print(x)
            f.write("{0} {1} {2}\n" .format(x[0], x[1], connectiontypes[i]))
        for i in connections:
            for j in i:
                f.write("{0} " .format(str(j)))
            f.write('\n')
    f.write("# Degrees of freedom (optional)\n")
    if len(atoms) == 1:
        dof = 3
    else:
        dof = 6
    f.write("Molecule_DOF: {0} #3 for monoatomic, 5 for rigid linear and 6 for  rigid nonlinear, if flexible then generally 3*number of atoms " .format(dof))        
    f.close()
    print("Molecule file {0}{1}_v{2}.mol written!" .format(sourcedir, species, interactions_version))

def directorymaker(dxout = "./"):
    filename = "{0}test.txt" .format(dxout)
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    with open(filename, "w") as f:
        f.write("FOOBAR")
    print("Directory written!")
    
def AtmAtmMover(version, elements, dxout = "./mapgen/test/"):
    tot_el_list = []
    tot_el_list.append(MOF_el_list)
    tot_el_list.append(elements)
    hicut = 18
    with open("{0}atom_atom_file" .format(dxout), "w") as file:
        file.write("""#This is an autogenereated interactions file for music, based on the python script Umbrella-v2. It's probably best to look this over before running your sims\n#Lennard-Jones interactions\n""")
        for i in range (0, len(MOF_el_list)):
            file.write("""{0} {0} LJ SIG@{1} EPS@{2} HICUT@{3}\n""".format(MOF_el_list[i], IntParams[MOF_el_list[i]][0], IntParams[MOF_el_list[i]][1], hicut))
            for j in range (i+1, len(MOF_el_list)):
                file.write("""{0} {1} LJ OFF\n""".format(MOF_el_list[i], MOF_el_list[j]))
            for k in elements:
                file.write("""{0} {1} LJ SIG@LBMIX EPS@LBMIX HICUT@{2}\n""".format(MOF_el_list[i], k, hicut))
            file.write("\n")
        for i in range (0, len(elements)):
            raw_el = str(elements[i].split('_')[0] + '_' + species)
            file.write("""{0} {0} LJ SIG@{1} EPS@{2} HICUT@{3}\n""".format(elements[i], IntParams[raw_el][0], IntParams[raw_el][1], hicut))
            for j in range (i+1, len(elements)):
                file.write("""{0} {1} LJ SIG@LBMIX EPS@LBMIX HICUT@{2}\n""".format(elements[i], elements[j], hicut))
            file.write("\n")
        file.write("#Coulombic region\n")
        for i in range (0, len(MOF_el_list)):
            file.write("""{0} {0} COUL OFF\n""".format(MOF_el_list[i]))
            for j in range (i+1, len(MOF_el_list)):
                file.write("""{0} {1} COUL OFF\n""".format(MOF_el_list[i], MOF_el_list[j]))
            for k in elements:
                file.write("""{0} {1} COUL OFF\n""".format(MOF_el_list[i], k))
            file.write("\n")
        for i in range (0, len(elements)):
            file.write("""{0} {0} WFCOUL  HICUT@{1}\n""".format(elements[i], hicut))
            for j in range (i+1, len(elements)):
                file.write("""{0} {1} WFCOUL HICUT@{2}\n""".format(elements[i], elements[j], hicut))
            file.write("\n")
    print("Atom-Atom interaction file written!")


def SorbSorbWriter(species, elements, version, dxout = "./mapgen/", pmap = True, emap = False):
    file = "sorb_sorb_file"
    with open("{0}{1}" .format(dxout, file), 'w') as file:
        file.write("{0}_v{1} {0}_v{1} NCOUL OFF\n{0}_v{1} {0}_{1} COUL OFF\n\n" .format(framework, version)) #MOF-MOF interactions
        if isinstance(species, list) == True: 
            for i in species:            #swap these if and for statements?
                if pmap == True:
                    file.write("{0}_v{2} {1}_v{2} NCOUL MAP@{1}_v{2} FAST ".format(i, framework, version))
                    for j in elements:
                        file.write("{0}@PMAP@{2}.{0}_v{1}.pmap " .format(j, version, framework))
                    file.write("\n")
                else:
                    file.write("{0}_v{2} {1}_v{2} NCOUL BASIC LJ FAST\n".format(i, framework, version))
            for i in species:
                if emap == True:
                    file.write("{0}_v{2} {1}_v{2} COUL MAP@{1}_v{2} FAST ".format(i, framework, version))
                    for j in elements:
                        file.write("{0}@EMAP@{2}_v{1}.{0}_v{1}.emap " .format(j, version, framework))
                    file.write("\n")
                else:
                    file.write("\n{0}_v{2} {1}_v{1} COUL OFF\n" .format(i, framework, version))
                    print("Framework-fluid coulombic interactions are off")
            file.write("\n{0}_v{1} {0}_v{1} NCOUL BASIC LJ FAST\n{0}_v{1} {0}_v{1} COUL BASIC WFCOUL FAST\n\n" .format(i, version))
        elif isinstance(species, str) == True: 
            if pmap == True:
                file.write("{0}_v{2} {1}_v{2} NCOUL MAP@{1} FAST ".format(species, framework, version))
                for j in elements:
                    file.write("{0}@PMAP@{2}.{0}_v{1}.pmap " .format(j, version, framework))
                file.write("\n\n")
            else:
                file.write("{0}_v{2} {1}_v{2} NCOUL BASIC LJ FAST".format(species, framework, version))
                file.write("\n\n")
            if emap == True:
                file.write("{0}_v{2} {1}_v{2} COUL MAP@{1} FAST ".format(species, framework, version))
                for j in elements:
                    file.write("{0}@EMAP@{2}_v{1}.{0}_v{1}.emap " .format(j, version, framework))
                file.write("\n\n")
            else:
                file.write("{0}_v{2} {1}_v{2} COUL OFF\n" .format(i, framework, version))
                print("Framework-fluid coulombic interactions are off")
        file.write("\n{0}_v{1} {0}_v{1} NCOUL BASIC LJ FAST\n{0}_v{1} {0}_v{1} COUL BASIC WFCOUL FAST\n\n" .format(species, version))
    print("Sorb-Sorb file written!")    

def IntraWriter(species, elements, version, dxout = "./mapgen/"):
    with open("{0}intramolecular_file" .format(dxout), "w") as file:
        file.write("####This intramolecular interaction file is workable for making a map of {0} in a 2x2x2 IRMOF-1 framework\n" .format(species))
        if isinstance(species, list) == True: 
            for i in species:
                file.write("Intra: {0}_v{1}\n" .format(i, version))
        elif isinstance(species, str) == True: 
            file.write("Intra: {0}_v{1}\n" .format(species, version))
        for i in list(elements):
            file.write("Intra: {0}_v{1}\n".format(i, version))
        file.write("Intra: {0}" .format(framework))
    print("Intra file written!")

def Intsetup(mol_name, interactions_version, elements, dxout, pmap = True, emap = False):
    AtmAtmMover(interactions_version, elements, dxout)
    SorbSorbWriter(mol_name, elements, interactions_version, dxout, pmap, emap)
    IntraWriter(mol_name, elements, interactions_version, dxout)

def MapMakeRunWriter(species, elements, version, dxout = "./mapgen/"):
    with open("{0}run.mapmaker" .format(dxout), "w") as file:
        file.write("""#!/usr/bin/env bash
#SBATCH --job-name={0}.map
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jrhm21@bath.ac.uk
#SBATCH account=free

# -- Set up the environment
module purge
module load group ce-molsim stack
module load music/std

# -- Find my directory
cd ~/scratch
cd 01_music_IRMOF_solvent/{1}\n
export ATOMSDIR=../../atoms
export MOLSDIR=../../molecules
export PMAPDIR=../../maps
export EMAPDIR=../../maps
\n\n# --Run\n""" .format(species, dxout.split(".")[1]))
        for count, i in enumerate(elements, 1):
            file.write("music_mapmaker makemap_{0}_v{1}.ctr > logfile{2}\n" .format(i, version, count))
    os.chmod("{0}run.mapmaker" .format(dxout), 0o777) #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    print("Runfile written!")   

def mapctrlfilewriter(species, version, dxout = "./mapgen/"):
    x = random.randint(0,99999)
    with open("{0}makemap_{1}_v{2}.ctr" .format(dxout, species, version), "w") as file:
        file.write("""------ General Information ------------------------------------------
{0}_pmap_v{1} map probe in IRMOF1
1               # No. of iterations
1               # No. of steps between writes to output/log file
1               # No. of steps between writes to crash file
1               # No. of steps between writes to config. file
2                    # Start numbering simulations from .
{2}              # Iseeed
4                    # specifies contents of config file,
{0}_v{1}.res           # Restart File to write to
{0}_v{1}.con           # Configuration File
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

{0}_v{1}                   # sorbate 
{0}_v{1}.mol               # sorbate coordinates file

IRMOF1                    # sorbate 
IRMOF1.mol                # sorbate coordinates file
------ Simulation Cell Information --------------------------------------
IRMOF1                   # Fundamental cell type
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

IRMOF1           # Sorbate to map
{0}_v{1}       # Sorbate to probe map with
NCOUL LJ       # Interaction type to map
0.2            # Approximate grid spacing (Ang)
100.0          # High end potential cutoff (kJ/mol)
IRMOF1.{0}_v{1}.pmap           # Map filename or AUTO
------ Configuration Initialization -------------------------------------
{0}_v{1}                            # Sorbate_Type  
Molecule NULL                              # Source Filename
IRMOF1                            # Sorbate_Type
Fixed NULL                     # Source Filename""" .format(species, version, x))
    print("Map control file written!")
    
def pSat(Species, T):#This function checks the Species is there and that you're in the right temperature range. 
    Chemical = Species.split("_")[0]
    print(Chemical)
    Press = 1 #pressure defaults to 1 kPa
    if Chemical in Antoine:                              #It then reads out your saturation pressure in kPa for the user benefit and to get the rest of the program working
        if Antoine[Chemical][3] <= T <= Antoine[Chemical][4]:
            print("I have Antoine parameters for {0} at {1} K." .format(Chemical, T))
        elif T < Antoine[Chemical][3]:
            print("WARNING, that temperature is too low for my Antoine parameters for {0} (my limit is {1} K). I'll persevere anyway." .format(Chemical, Antoine[Chemical][3]))
        elif Antoine[Chemical][4] < T:
            print("WARNING, that temperature is too high for my Antoine parameters for {0} (my limit is {1} K). I'll persevere anyway." .format(Chemical, Antoine[Chemical][4]))
        else:
            print("Something weird happened, I've probably got a bug. Oops!")
        Press = 100*(10**(Antoine[Chemical][0]-(Antoine[Chemical][1]/(Antoine[Chemical][2]+T))))
    else:
        print("I'm sorry, I don't have that Species in my database.")
    return Press

def isothermcalculator(Pressure, n):                                  #This function produces the isotherm as a data list, using pSat calculated by function pSat and the user defined isotherm length n
    m = n-len(FixedPoints)-1                                         #Pressure points re linearly distributed above 0.09 kPa to pSat, which may eb a bad idea. Who knows?
    AllPoints = []
#    print(FixedPoints[-1])
#    print(Pressure-FixedPoints[-1])
    if FixedPoints[-1]/Pressure < 0.001:                  #This subroutine is used when the pSat isn't much higher than the fixed points. At this point it sets the minimum relative pressure to 10**-5 and makes a log linear range between this and 0 (+ a bit)
        minrelpress = -5
        for point in range(1,n+1):
            relpress = minrelpress*(1-(point/float(n+1)))
            press = Pressure*10**relpress
            press = round(press, 3)
            AllPoints.append(press)
            AllPoints.sort()
    else:
        minrelpress = log((FixedPoints[-1])/Pressure, 10)
        for thing in FixedPoints:                                 #This is the general loop for printing pressures the rest of the time
            AllPoints.append(thing)
        for point in range(1,m+2):
            relpress = minrelpress*(1-(point/float(m+1)))
            press = Pressure*10**relpress
            press = round(press, 3)
            AllPoints.append(press)
            AllPoints.sort()
    isotherm = AllPoints
    print("PSat = {0}KPa, Isotherm = {1}".format(Pressure, isotherm))
#    print(isotherm[2])
#    print(len(isotherm))
    return isotherm
    

def PressureFileWriter(Species, T, Pressure, isotherm, dirout, version):                 # This function writes your file to a specific directory, but currently doesn't say if your temp is out of range. Annoying!!#
#    print(isotherm)
    with open("{0}pressure.{1}.{2}.dat" .format(dirout, Species, T), "w") as file:
        file.write("{0!s}_v{3} #PSat at {1!s} K is {2:1.3f} kPa. \n" .format(Species, T, Pressure, version))
        file.write(str(len(isotherm)) + "\n")
        file.write(', '.join(str(thing) for thing in isotherm))
    
def TaskfarmRunWriter(Species, T, version, dirout):
    with open("{0}run.taskfarmer" .format(dirout), 'w') as file:
        file.write("""#!/usr/bin/env bash
#SBATCH --job-name={0}.IRMOF1.{1}k
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --time=06:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jrhm21@bath.ac.uk
#SBATCH account=free

# -- Set up the environment
module purge
module load group ce-molsim stack
module load taskfarmer

#locate the experiment directory
cd ~/scratch
cd 01_music_IRMOF_solvent/
cd {3}

# -- Run
mpirun -np 16 taskfarmer -f taskfarm

# -- python analyse the results
python ./isothermextractor.py""".format(Species, T, version, dirout))
    os.chmod("{0}run.taskfarmer" .format(dirout), 0o777)   #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    with open("{0}taskfarm" .format(dirout), 'w') as file1:
        for i in range(1,17):
            file1.write("bash ./{0:02d}/run.gcmc\n".format(i))
    with open("{0}taskfarm.backup" .format(dirout), 'w') as file2:
        for i in range(1,17):
            file2.write("bash ./{0:02d}/run.gcmc\n".format(i))
    print("Taskfarm things file written!")

def IsothermExtractMover(Species, T, version, dirout):
    with open("{0}isothermextractor.py" .format(dirout), 'w') as file:
        file.write("""############################################################################################################################################################################
####This python file takes a post-output and makes some xyz files for meta-analysis
####The 3 variables are molecule name ("your_species_here"), Temperature ("your_temperature_here"), and directory ("your_tree_here")
#######################################################################################################################################################

from os import listdir
import numpy as np

isotherms = np.zeros((17, 20)) # this is a matrix of a 16x 20-point isotherms with identical pressure values

def isothermextract(page): #this function scrapes subdirectories of name for isotherm data for a given isotherm file name
    z = 0
    z = int(page)
#    print(z)
    drx = str("./{{0:02d}}/isotherm.{0}" .format(page)) # defines the isotherm
#    print(drx)
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
for i in range (1,17):
    SuccessCount = 0
    filenames = listdir("./{{0:02d}}/" .format(i))
    for file in filenames:
        if file == "isotherm.{0}":
            isothermextract(i)
            SuccessCount += 1
    if SuccessCount != 1:
        print("WARNING! Something went wrong in simulation {{0:02d}}" .format(i))
x = isotherms[0]
avg=[]
stdev=[]

avgisotherms = np.zeros((3, 20))
for i in range (0, 20):
    spread = []
    for thing in range (2,17):
        if isotherms.item((thing, i)) != 0:
            spread.append(isotherms.item((thing, i)))
    avg.append(np.mean(spread))
    stdev.append(np.std(spread))
avgisotherms[0] = x
avgisotherms[1] = avg
avgisotherms[2] = stdev
np.savetxt("../{1}.{0}_v{2}.results.csv", avgisotherms, delimiter=",")""" .format(Species, T, version))
    print("Python isotherm extractor written!")

def GcmcControlChanger(species, elements,  T, n, version, dirout):
    x = random.randint(0,99999)
    n_species = 4+len(elements)
    with open("{0}gcmc.ctr" .format(dirout), 'w') as file:
        file.write("""#This control file was written by the python scrpt umbrellav2
------ General Information ------------------------------------------
{0}_v{1} molecule in IRMOF1_v{1}
500000              # No. of iterations
50000                # No. of steps between writes to output/log file
100000                # No. of steps between writes to crash file
2500                  # No. of steps between writes to config. file
1                   # Start numbering simulations from .
{2}                # Iseed
3                    # specifies contents of config file,
IRMOF1_v{1}.{0}_v{1}.res         # Restart File to write to
IRMOF1_v{1}.{0}_v{1}.con          # Configuration File
------ Atomic Types --------------------------------------------------
{3}                    # number of atomic types            \n\n""".format(species, version, x, n_species))
        for i in elements:
            file.write("""{0}   #atom type\n{0}.atm #atom file name\n\n""".format(i))
        file.write("""Carbon              # atom type
Carbon.atm          # basic atom info file

Oxygen               # atom type
Oxygen.atm           # basic atom info file

Hydrogen               # atom type
Hydrogen.atm           # basic atom info file

Zinc              # atom type
Zinc.atm           # basic atom info file

------ Molecule Types -------------------------------------------------
2                    # number of sorbate types

{0}_v{1}               # sorbate
{0}_v{1}.mol           # sorbate coordinates file

IRMOF1_v{1}                # sorbate
IRMOF1_v{1}.mol             # sorbate coordinates file
------ Simulation Cell Information ------------------------------------
IRMOF1_v{1}                # Fundamental cell file
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
{0}_v{1}              # Sorbate Name
------ GCMC Information -----------------------------------------------
1                 # No. of iterations
{2}              # temperature
Ideal Parameters   # Tag for the equation of state (NULL = Ideal Gas)
{3}                  # No. of simulation points
5000                # Block size for statistics
1                  # no. of sorbates
          -------------------------
{0}_v{1}            # Sorbate Name
pressure.{0}.{2}.dat           #  pressure
Null               # sitemap filename (Null = no sitemap)
4                  # no of gcmc movetypes
1.0, 1.0, 1.0, 1.0      # move type weights
BINSERT                   # type of move.1
IRMOF1.{4}_v{1}.pmap  # Bias Potential File, CAVITY-> Implies cavitybias
298.                       # Bias temperature for the bmap
BDELETE                   # type of move.2
RTRANSLATE                # type of move.4
0.2, 1                    # Delta Translate, adjust delta option (0=NO, 1=YES)
RROTATE
0.2, 1
------ Configuration Initialization -------------------------------------
{0}_v{1}             # Sorbate_Type
GCMC NULL
IRMOF1_v{1}              # Sorbate_Type
FIXED NULL
--------  Main Datafile Information --------
Energy, position, pair_energy  # contents of datafile""".format(species, version, T, n, elements[-1]))
    print("GCMC control file written!")

def PostControlChanger(Species, T, n, version, dirout):
    with open("{0}post.ctr" .format(dirout), 'w') as file:
        file.write("""####This section is apparently required for working with any post code. Who knows why? not me!
####This post control file is workable for chloroalkanes on a 2x2x2 IRMOF-1 only right now
####The 3 variables are molecule name ("your_species_here"), Temperature ("your_temperature_here"), and number of simulation points ("your_number_of_points_here")
------------------------------------------------------------
   ### Required section ######
-- Post Processor Information ------------
GCMC                            # Type of simulation GCMC, NVTMC , MD ....
./IRMOF1_v{1}.{0}_v{1}.con                    # basename for config files
1, {3}                          # first and last file numbers
post.ctr.out                       # name for new ctrlfile that will regenerated
IRMOF1_v{1}.{0}_v{1}.{2}K.post          # Base name for output files
40, 0                         # Percentages of data to skipped at start and end 


# The sections below are necessary only if you want the corresponding 
# analysis performed
# ---------------- ALL OF THEM ARE OPTIONAL ------------------------


####    This section is reqd for energy averages in your post code output files
####    as of now only total enrgies vs sim. step
------ Post : Energy Average Info -----------------------------------
20       # Number of blocks into which data should be divided for stats

####    This section is reqd for Loading averages in your post code outputfiles
####    as of now only species loading vs sim. step (for all species)
------ Post : Loading Average Info -----------------------------------
20       # Number of blocks into which data should be divided for stats""".format(species, version, T, n))
    print("Post control file written!")

def GcmcRunWriter(Species, T, version, dirout):
    with open("{0}run.gcmc".format(dirout), 'w') as file:
        file.write("""#########################################################################
####This control file is workable for chloroalkanes on a 2x2x2 IRMOF-1 only right now
####The 3 variables are molecule name ("your_species_here"), Temperature ("your_temperature_here"), and directory ("your_tree_here")
#############################################################################################################################################################################
#!/usr/bin/env bash
#SBATCH --job-name={0}_v{1}.{2}
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jrhm21@bath.ac.uk
#SBATCH account=free

# -- Set up the environment
module purge
module load group ce-molsim stack
module load music/std

#locate the experiment directory
cd ~/scratch
cd 01_music_IRMOF_solvent/
cd {3}

#set the paths for atoms, molecules, pmaps, and emaps
export ATOMSDIR=../../../../atoms
export MOLSDIR=../../../../molecules
export PMAPDIR=../../../../maps
export EMAPDIR=../../../../maps


# -- Run
music_gcmc gcmc.ctr > logfile.gcmc
music_post post.ctr > logfile.post


""".format(Species, version, T, dirout))
    os.chmod("{0}run.gcmc" .format(dirout), 0o777) #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    print("Run file written!")

    
##############################
##############################
##############################
#This bit looks for your molecule files and makes one in required
print("Checking mol files...")
directorymaker('./molecules')
if filecheck(species, interactions_version, 'mol', './molecules')[0] == True:
    if filecheck(species, interactions_version, 'mol', './molecules')[1] == True:
        print('{0} already exists as a molecule file with version {1}.' .format(species, interactions_version))
    elif filecheck(species, interactions_version, 'mol', './molecules')[1] == False:
        print("{0} already exists as a molecule file but with the wrong version. I'll need to make a new molecule file for {0}" .format(species))
        if filecheck(species, interactions_version, 'pdb')[0] ==True:
            atomlist, connlist = dataextract(species)
            bondinfo = connectiontypeswrite(atomlist, connlist)
            molwrite(species, atomlist, bondinfo, connlist, interactions_version, './molecules/')
        else:
            print("I'd need a .pdb file to get started though!")
else:
    print("{0} doesn't exist as a molecule file at all! I'll need to make a new molecule file for {0}" .format(species))
    print(filecheck(species, interactions_version, 'pdb')[0])
    if filecheck(species, interactions_version, 'pdb', './')[0] == True:
        atomlist, connlist = dataextract(species)
        bondinfo = connectiontypeswrite(atomlist, connlist)
        molwrite(species, atomlist, bondinfo, connlist, interactions_version, './molecules/')
    else:
        print("I'd need a .pdb file to get started though!")
        
#This bit now looks for all of your map files
print("Checking for maps of elements within {0}molecules/{1}_v{2}.mol" .format('./', species, interactions_version))
el_list = set()
if filecheck(species, interactions_version, 'mol', './molecules')[0] == True:
    with open("{0}molecules/{1}_v{2}.mol" .format('./', species, interactions_version)) as file:
        for line in file:
            if len(line.split()) == 8:
#                print(line.split()[4])
                el_list.add(line.split()[4])
    print("Atom list found is: " + str(el_list))
else:
    print("That's weird, I can't find your mol file. Something must have messed up!")

sorb_el_list = []
for i in el_list:
    sorb_el_list.append("{0}_sorb" .format(i.split("_")[0]))
print(sorb_el_list)
for i in sorb_el_list:
    #print(i)
    directorymaker('./maps')
    if mapcheck(i, interactions_version, 'pmap', './maps', framework)[0] == True:
        print("Interactions pmap found between fluid atom {0} and a framework {1}" .format(i, framework))
    else:
        print("{0} map on {1} doesn't exist with interactions version {2}! You'll have to make a new one!" .format(i, framework, interactions_version))
        if filecheck(i, interactions_version, 'mol', './molecules')[0] == True:
            if filecheck(species, interactions_version, 'mol', './molecules')[1] == True:
                print('{0} already exists as a molecule file with version {1}.' .format(species, interactions_version))
            elif filecheck(species, interactions_version, 'mol', './molecules')[1] == False:
                print("{0} already exists as a molecule file but with the wrong version. I'll need to make a new molecule file for {0}" .format(species))
                molwrite("{0}" .format(i), [i], None, None, interactions_version, "./molecules/")
        else:
            print("No map molecule found for {0}_v{1}, I'll have to write a new one!" .format(i, interactions_version, 'mol'))
            molwrite("{0}" .format(i), [i], None, None, interactions_version, "./molecules/")
        directorymaker('./mapgen/{0}_v{1}/' .format(species, interactions_version))
        Intsetup(sorb_el_list, interactions_version, sorb_el_list, './mapgen/{0}_v{1}/' .format(species, interactions_version), False)
        MapMakeRunWriter(i, sorb_el_list, interactions_version, './mapgen/{0}_v{1}/' .format(species, interactions_version))
        mapctrlfilewriter(i, interactions_version, './mapgen/{0}_v{1}/' .format(species, interactions_version))
        print("Your map is ready to be generated! Just go to ./mapgen/{0}_v{1}/ and run run.mapmaker" .format(species, interactions_version))
#This bit defines your pressure information
satP = pSat(species, T)
istm = isothermcalculator(satP, n)
#This bit now sets up your experiment. WARNING! It overwrites all rpevious experiment files in the directory chosen (but significantly doesn't overwrite any results files)        
directorymaker('./experiments/{0}_v{1}/{2}/'.format(species, interactions_version, T))
print("Working in the general directory:")
TaskfarmRunWriter(species, T, interactions_version, './experiments/{0}_v{1}/{2}/'.format(species, interactions_version, T))                                                          #Writes the taskfarmer run file
IsothermExtractMover(species, T, interactions_version, './experiments/{0}_v{1}/{2}/'.format(species, interactions_version, T))	                                                           #Writes isothermextractor.py into place
for directory in range (1,17):
    directorymaker('./experiments/{0}_v{1}/{2}/{3:02d}/'.format(species, interactions_version, T, directory))
    Intsetup(species, interactions_version, sorb_el_list, './experiments/{0}_v{1}/{2}/{3:02d}/'.format(species, interactions_version, T, directory), True, True)
    GcmcControlChanger(species, sorb_el_list, T, n, interactions_version, './experiments/{0}_v{1}/{2}/{3:02d}/'.format(species, interactions_version, T, directory))                            #puts in the gcmc control file
    PostControlChanger(species, T, n, interactions_version, './experiments/{0}_v{1}/{2}/{3:02d}/'.format(species, interactions_version, T, directory))                            #puts in the post control file
    GcmcRunWriter(species, T, interactions_version, './experiments/{0}_v{1}/{2}/{3:02d}/'.format(species, interactions_version, T, directory)) #puts in the runfile
    PressureFileWriter(species, T, satP, istm, './experiments/{0}_v{1}/{2}/{3:02d}/'.format(species, interactions_version, T, directory), interactions_version)
        
    
