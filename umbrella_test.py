#This python script is intended to be a one-stop shop for setting up/running a music simulation, so you can automate your processes in peace
################Specifically, this one works for mixed isotherms at constant fugacity###############################
#So this script will just take in all required inputs, and print template files for map generation (if required) and simulations.
#The workflow goes like this: 
#1. Set global variables for things to be softcoded in
#2. Take in user input variables
#3. If the user inputs a pdb file, create a molfile from this
#4. write the appropraite pair interaction files (atom_atom, sorb_sorb, intra)
#5. Make an isotherm at your temperatuure of x points, equally spaced for a semilog graph
#6. Write the above out to the appropriate directories
#7. Report on what it's done.
##########
#0. Preamble
import os
import re
import random
import numpy as np
from math import sqrt, log
from pathlib import Path
import logging
import datetime
now = datetime.datetime.now()
#try:
#    input = raw_input("What is your sorbent species called?")
#except NameError:
#    pass

##########
#1. Your required information
element_names = {'H': 'Hydrogen', #In my PDB files, they're as symbols, and in my music files, they're names. This dictionary allows for interconversion
                 'C': 'Carbon',
                 'Cl': 'Chlorine',
                 'O': 'Oxygen',
                 'N': 'Nitrogen',
                 'F': 'Fluorine',
                 'CH': 'Methynyl'
}
Antoine = {                                                 #This is your library of Antoine equation parameters. The key is the name of the Species that will be written
    "CCl4":(4.02291, 1221.781, -45.739, 193, 350),          #Indices 0-2 are antoine parameters A, B, and C taken from NIST,
    "Chloroform":(4.20772, 1233.129, -40.953, 215, 334), #Indices 3 and 4 are the Antoine equation validity range, as quoted on the NIST website
    "DCM":(4.53691,1327.016,-20.474, 233, 313),
    "Chloromethane":(4.91858, 1427.529, 45.137, 303, 416),
    "Chloromethane2":(4.22507, 951.561, -23.468, 198, 278), 
    "Methane":(4.22061, 516.689, 11.223, 110, 190),
    "Methanol":(5.20409, 1581.341, -33.5, 288, 356.8),
    "THF":(4.12118, 1202.942, -46.818, 296.29, 372.8),
    "DMF":(3.93068, 1337.716, -82.648, 303, 363),
    "Nitrogen":(3.7362, 264.651, -6.788, 63.14, 126)
    }

IntParams = {	 #Sig (ang), Eps(K), Q (e)   #Here are your atom-atom parameters for your forcefield. It's listed as atomtype_molecule or atomtype_forcefield. Check this before you do anything!
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
	"Carbon_Kam":(3.41, 68.94, -0.235), #Chlorfoorm parametrs taken from 10.1021/jp0535238
	"Hydrogen_Kam":(2.81, 10.06, 0.355),
	"Chlorine_Kam":(3.45, 138.58, -0.04),
	"Nitrogen_TRAPPE":(3.310, 36.0, -0.482),
	"COMN_TRAPPE":(0, 0, 0.964),
        "Methynyl_OPLS": (3.8, 40.26, 0.42),
        "Chlorine_OPLS": (3.47, 150.98, -0.14),
        "Hydrogen_CDP": (2.81, 10.06, -0.0551),
        "Carbon_CDP": (3.41, 68.94, 0.5609),
        "Chlorine_CDP": (3.45, 138.58, -0.1686),
	"Chlorine_UFF": (3.947, 114.128, 0),
        "Carbon_UFF": (3.851, 52.79, 0),
        "Hydrogen_UFF": (2.886, 22.12, 0)
        }

Components = {
        'DMF':('Me_DMF', 'Hydrogen_DMF', 'Oxygen_DMF', 'Nitrogen_DMF', 'Carbon_DMF'),
        'Chloroform':('Carbon_Kam', 'Chlorine_Kam', 'Hydrogen_Kam'),
        'DCM':('Carbon_Kam', 'Chlorine_Kam', 'Hydrogen_Kam'),
        'Chloromethane': ('Carbon_Kam', 'Chlorine_Kam', 'Hydrogen_Kam'),
        'CCl4': ('Carbon_Kam', 'Chlorine_Kam')
}

MOF_el_list = ["Carbon", "Hydrogen", "Oxygen", "Zinc"] #This is a list of elements in your MOF, needed for the control file atom types section
FixedPoints = [0.01, 0.03, 0.05, 0.07, 0.09]              #These are the fixed low pressure points I'll include at the start of your isotherm
errormessages = []

#############
##2. Input variables:
#Framework - the name of the MOF to be used. Combined with version to get a molfile
#framework = input("What is your framework species called?") #This line and the next denote the user will input it as part fot he script. Good for changing things quickly
#type(framework) #e.g. 'IRMOF1' (with apostrophes)
framework = 'IRMOF1'
##Species - your sorbent. Only 1 for now
species = "DMF"
species1 = "Chloroform"
#species = input("What is your sorbent species called?")
#type(species)
##Temperature
T = 298 
#T = eval(input("What temperature (in K) would you like to test?"))
#type(T)
##number of isotherm points
n = 20
#n = eval(input("How many isotherm points do you want?"))
#type(n)
##Forcefield name
forcefield = "DMF"
forcefield1 = 'Kam'
#n = input(input("What is yoru forcefield called?"))
#type(forcefield)
parentdir = Path('/home/r/jrhm21/scratch/05_music_mixed_sorbent_isotherms/')
xptpath = Path('./experiments/{0}.{2}/{1}/'.format(species, T, species1))
minrelpress = -3
maxrelpress = -0.5

#placeholder
filename = "{0}/test.txt" .format(xptpath) #Test file name
if not os.path.exists(os.path.dirname(filename)): #Checks if the test file exists
    try:
        os.makedirs(os.path.dirname(filename)) #Makes the file
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
    with open(filename, "w") as f:
        f.write("FOOBAR") #Writes something in the test file


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
            f.write("{0} 0 0\n" .format('0' if forcefield == None else str(IntParams['{0}_{1}'.format(i[-1].split('_')[0], forcefield)][2]))) #prints q, ?, ?. change the '0' in format if you want to set a q to a variable
    else:
        f.write('1 0 0 0 {0} ' .format(str(atoms[-1]))) #prints the atom index, x, y, z, name
        f.write("{0} 0 0\n" .format('0' if forcefield == None else str(IntParams['{0}_{1}'.format(atoms[-1].split('_')[0], forcefield)][2]))) #prints q, ?, ?. change the '0' in format if you want to set a q to a variable
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
def AtmAtmMover(elements, dxout = "./mapgen/test/"):
    hicut = 18 #sets your hicut in A
    with open("{0}atom_atom_file" .format(dxout), "w") as file:
        file.write("""#This is an autogenereated interactions file for music, based on the python script Umbrella-v3. It's probably best to look this over before running your sims\n#Lennard-Jones interactions\n""")
        for i in range (0, len(MOF_el_list)):
            file.write("""{0} {0} LJ SIG@{1} EPS@{2} HICUT@{3}\n""".format(MOF_el_list[i], IntParams[MOF_el_list[i]][0], IntParams[MOF_el_list[i]][1], hicut)) #LJ self-parameters for the MOF-MOF interactions
            for j in range (i+1, len(MOF_el_list)):
                file.write("""{0} {1} LJ OFF\n""".format(MOF_el_list[i], MOF_el_list[j])) #LJ parameters for MOF 'i j' pairs, which are generally off
            for k in elements:
                file.write("""{0} {1} LJ SIG@{3} EPS@{4} HICUT@{2}\n""".format(MOF_el_list[i], k, hicut, 
str(round((IntParams[MOF_el_list[i]][0]+IntParams[k][0])/2, 3)), str(round(sqrt(IntParams[MOF_el_list[i]][1]*IntParams[k][1]), 3)))) #MOF-fluid forcefield paramters
            file.write("\n")
        for i in range (0, len(elements)):
            file.write("""{0} {0} LJ SIG@{1} EPS@{2} HICUT@{3}\n""".format(elements[i], IntParams[elements[i]][0], IntParams[elements[i]][1], hicut)) #LJ self-parameters for the fluid-fluid interactions
            for j in range (i+1, len(elements)):
                logger.debug('i = '+ str(elements[i])+ ', j= '+str(elements[j]) )
                file.write("""{0} {1} LJ SIG@{3} EPS@{4} HICUT@{2}\n""".format(elements[i], elements[j], hicut,
 str(round((IntParams[elements[i]][0]+IntParams[elements[j]][0])/2, 3)), str(round(sqrt(IntParams[elements[i]][1]*IntParams[elements[j]][1]), 3)))) #LJ parameters for the fluid-fluid 'i j' interactions
            file.write("\n")
        file.write("#Coulombic region\n")
        for i in range (0, len(MOF_el_list)):
            file.write("""{0} {0} COUL OFF\n""".format(MOF_el_list[i])) #MOF-MOF self-coulombic interactions (off)
            for j in range (i+1, len(MOF_el_list)):
                file.write("""{0} {1} COUL OFF\n""".format(MOF_el_list[i], MOF_el_list[j])) #MOF-MOF interactiosn (also off)
            for k in elements:
                file.write("""{0} {1} COUL OFF\n""".format(MOF_el_list[i], k)) #MOF-fluid1 pairwise interactions (definitely off, use a map for this)
            file.write("\n")
        for i in range (0, len(elements)):
            file.write("""{0} {0} WFCOUL  HICUT@{1}\n""".format(elements[i], hicut)) #fluid fluid self-interactions - use wolf cpoulombic. partial charge values are stored in the molecule files - see molfilewriter
            for j in range (i+1, len(elements)):
                file.write("""{0} {1} WFCOUL HICUT@{2}\n""".format(elements[i], elements[j], hicut)) #fluid fluid 'i j' interactions - use wolf cpoulombic. partial charge values are stored in the molecule files - see molfilewriter
            file.write("\n")
    logger.debug("Atom-Atom interaction file written!")

#This writes your sorb-sorb file telling music which intermolecular interactions are on. It handles 3 cases:
#1 Map generation (species = list, elements = list, ****emap and pmap are both off****)
#2 simple X on Sorbent (species = string, elements  = list, emap and pmap can be on)
#3 multiple species on Sorbent (species  = list, elements = nested list (elements by species) emap and pmap can be on
def SorbSorbWriter(species, elements, framework, dxout = "./mapgen/", pmap = True, emap = False):
    filename = "sorb_sorb_file"
    with open("{0}{1}" .format(dxout, filename), 'w') as file:
        file.write("{0} {0} NCOUL OFF\n{0} {0} COUL OFF\n\n" .format(framework)) #MOF-MOF interactions are off
        if isinstance(species, list) == True: #This checks if your species is a list or a string, which makes it work for multiple species at once
            if pmap == True: #using a MOF-fluid LJ potential map if true make sure it's false for map genereation
                for i in species:            
                    file.write("{0} {1} NCOUL MAP@{1} FAST ".format(i, framework)) #declares it'll use a map
                    for j in Components[i]: #if species is a list, elements should be a nested list like [species0, species1] and species0 = [element0, element1] etc.
                        file.write("{0}@PMAP@{1}.{0}.p " .format(j, framework)) #writes your map names. WARNING - it's a 18 character limit on map names, and no '-' allowed!
                    file.write("\n")
            else:
                for i in species:
                    file.write("{0} {1} NCOUL BASIC LJ FAST\n".format(i, framework)) #pairwise fluid-framework interactions
            if emap == True: #using a MOF-fluid coulomb map if true
                for i in species:
                    file.write("{0} {1} COUL MAP@{1} FAST ".format(i, framework)) #declares it'll use a map
                    for j in Components[i]:
                        file.write("{0}@EMAP@{1}.{2}.e " .format(j, framework, 'Probe')) #writes your map names. WARNING - it's a 18 character limit on map names, and no '-' allowed!
                    file.write("\n")
            else:
                for i in species:
                    file.write("\n{0} {1} COUL OFF\n" .format(i, framework)) #pairwise coulomb interactionsdon't work. Trust me, this is better.
                logger.warning("Framework-fluid coulombic interactions are off") #warns you your coulomb fluid-framework stuff isn't happeneing
            for i in species:
                file.write("\n{0} {0} NCOUL BASIC LJ FAST\n{0} {0} COUL BASIC WFCOUL FAST\n\n" .format(species[0])) 
                for j in range(1, len(species)):
                    file.write("\n{1} {0} NCOUL BASIC LJ FAST\n{1} {0} COUL BASIC WFCOUL FAST\n\n" .format(i, species[j]))
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
            file.write("\n{0} {0} NCOUL BASIC LJ FAST\n{0} {0} COUL BASIC WFCOUL FAST\n\n" .format(species))
    logger.debug("Sorb-Sorb file written!")    

#Writes your intramolecular file
def IntraWriter(species, framework, dxout = "./mapgen/"):
    with open("{0}intramolecular_file" .format(dxout), "w") as file:
        file.write("####This intramolecular interaction file is workable for making a map of {0} in a 2x2x2 IRMOF-1 framework\n" .format(species)) #outdated preamble
        if isinstance(species, list) == True: #Lets you handle multiple fluids at once
            for i in species:
                file.write("Intra: {0}\n" .format(i)) #no intramolecular interactions at all
        elif isinstance(species, str) == True: 
            file.write("Intra: {0}\n" .format(species))
        file.write("Intra: {0}" .format(framework)) #No framework intramolecular interactions happening
    logger.debug("Intra file written!")

#bundles the above 3 functions into one for simplicity
def Intsetup(mol_name, elements, el_mix, framework, dxout, pmap = True, emap = False):
    AtmAtmMover(elements, dxout)
    SorbSorbWriter(mol_name, el_mix, framework, dxout, pmap, emap)
    IntraWriter(mol_name, framework, dxout)

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
#SBATCH account=free

# -- Set up the environment
module purge
module load group ce-molsim stack
module load music/std

# -- Find my directory
cd {2}{1}\n #Change this bit to get to your directory
export ATOMSDIR=../../atoms
export MOLSDIR=../../molecules
export PMAPDIR=../../maps/
export EMAPDIR=../../maps
\n\n# --Run\n""" .format(species, dxout.split(".")[-1], parentdir))
        for count, i in enumerate(elements, 1):
            file.write("music_mapmaker makemap_{0}.ctr > logfile{1}\n" .format(i, count))
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
    if Chemical in Antoine:                              #It then reads out your saturation pressure in kPa for the user benefit and to get the rest of the program working
        if Antoine[Chemical][3] <= T <= Antoine[Chemical][4]:
            logger.info("I have Antoine parameters for {0} at {1} K." .format(Chemical, T))
        elif T < Antoine[Chemical][3]:
            logger.warning("WARNING, that temperature is too low for my Antoine parameters for {0} (my limit is {1} K). I'll persevere anyway." .format(Chemical, Antoine[Chemical][3]))
        elif Antoine[Chemical][4] < T:
            logger.warning("WARNING, that temperature is too high for my Antoine parameters for {0} (my limit is {1} K). I'll persevere anyway." .format(Chemical, Antoine[Chemical][4]))
        else:
            logger.warning("Something weird happened, I've probably got a bug. Oops!")
        Press = 100*(10**(Antoine[Chemical][0]-(Antoine[Chemical][1]/(Antoine[Chemical][2]+T))))
    else:
        logger.warning("I'm sorry, I don't have that Species in my database.")
    return Press

#Creates a log-linear pressure isotherm using the saturation pressure from pSat and preamble defined number of points/predefined points 
def isothermcalculator(Pressure, n, minrelpress, maxrelpress):                                  #This function produces the isotherm as a data list, using pSat calculated by function pSat and the user defined isotherm length n
    m = n-len(FixedPoints)-1                                         #Pressure points re linearly distributed above 0.09 kPa to pSat, which may eb a bad idea. Who knows?
    AllPoints = [] #your list of total pressures
#    print(FixedPoints[-1])
#    print(Pressure-FixedPoints[-1])
    if FixedPoints[-1]/Pressure < 0.001:                  #This subroutine is used when the pSat isn't much higher than the fixed points. At this point it sets the minimum relative pressure to 10**-5 and makes a log linear range between this and 0 (+ a bit)
        #minrelpress = -2 #your minimum pressure points set as pSat*10^this
        for point in range(1,n+1):
            relpress = minrelpress+(point*maxrelpress/n) #Calculates your partial pressure list between minrelpress and 0. Dividing the (point/float(n+1)) statement by a number lowers your max pressure to a fraction of the sat pressure
            press = Pressure*10**relpress #converts the above into absolute pressures
            press = round(press, 3) #Rounds the float to 3 decimal points, for simplicity
            AllPoints.append(press)
            AllPoints.sort() #sorts your points into ascending order
    else:
        minrelpress = log((FixedPoints[-1])/Pressure, 10) #Now your minimum relative pressure is set as your top point of the fixed pressures list
        for thing in FixedPoints:                                 #This is the general loop for printing pressures the rest of the time
            AllPoints.append(thing)
        for point in range(1,m+2):
            relpress = minrelpress+(point*maxrelpress/n) 
            press = Pressure*10**relpress
            press = round(press, 3)
            AllPoints.append(press)
            AllPoints.sort()
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
#SBATCH --time=18:00:00
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

# -- python analyse the results\n""".format('.'.join(i for i in Species), T, framework, dirout, parentdir)) #uses my isotherm extractor script to pull all your isotherms out and put them into .csv files in the directory abocve the taskfarm run file
        for i in Species:
            file.write("""python isothermextractor.{0}.py\n""".format(i))
    os.chmod("{0}run.taskfarmer" .format(dirout), 0o777)   #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    with open("{0}taskfarm" .format(dirout), 'w') as file1: #Writes the commands to taskfarm into a separate taskfram file
        for i in range(1,17):
            file1.write("bash ./{0:02d}/run.gcmc > {0:02d}.log\n".format(i))
    with open("{0}taskfarm.backup" .format(dirout), 'w') as file2: #writes a backup you can restore the above to if there's a bug
        for i in range(1,17):
            file2.write("bash ./{0:02d}/run.gcmc > {0:02d}.log\n".format(i))
    logger.debug("Taskfarm things file written!")

#Writes a python script to extract data from your experiments after taskfarming
def IsothermExtractMover(Species, T, framework, dirout, n = 20):
    for i in Species:
        with open("{0}isothermextractor.{1}.py" .format(dirout, i), 'w') as file:
            file.write(
"""from os import listdir
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
np.savetxt("../{2}.{1}.{0}.{4}.csv", avgisotherms, delimiter=",") #output it to this file""".format(i, T, framework, n, 'To'.join(j for j in Species)))
    logger.debug("Python isotherm extractor written!")

#Writes a .ctr file for your gcmc simulations (production and for getting your final configuarions at the end). This currently cannot cope with multiple sorbents at once, creates a 2x2x2 cell, and has bias insert atom softcoded in
def GcmcControlChanger(species, elements, T, n, framework, dirout, iterations = '750000', Restart = None, name = 'gcmc.ctr', pressure = 'file'):
    logger.debug(elements)
    x = random.randint(0,99999) #sets your random seed
    if isinstance(elements, list): #lets you use multi-element sorbent molecules
        n_species = len(MOF_el_list)+len(elements) #sets the total number of element types in your gcmc
    else: #lets you use single-element sorbents
        n_species = len(MOF_el_list)+elements
    with open("{0}{1}" .format(dirout, name), 'w') as file:
        file.write("""#This control file was written by the python scirpt umbrellav3 You probably ought to check me before use!
------ General Information ------------------------------------------
{0} molecule in {1} 
{4}              # No. of iterations, defaults to 750000
50000                # No. of steps between writes to output/log file
100000                # No. of steps between writes to crash file
2500                  # No. of steps between writes to config. file
{5}                   # Start numbering simulations fromhere. defaults to 1 for production run, 30 for config runs (based on filename)
{2}                #random seed
3                    # specifies contents of config file, outdated?
{1}.{0}.res         # Restart File to write to
{1}.{0}.con          # Configuration File
------ Atomic Types --------------------------------------------------
{3}                    # number of atomic types            \n\n""".format('.'.join(i for i in species), framework, x, n_species, iterations, 1 if name == 'gcmc.ctr' else 30))
        for i in elements:
            file.write("""{0}   #atom type\n{0}.atm #atom file name\n\n""".format(i))
        for i in MOF_el_list:
            file.write("""{0}   #atom type\n{0}.atm #atom file name\n\n""".format(i))
        if isinstance(species, list):
            file.write(
"""------ Molecule Types -------------------------------------------------
{0}                    # number of sorbate types\n\n""".format(len(species)+1)) 
            for i in species:
                file.write("""{0}\n{0}.mol\n\n""".format(i))
        else:
            file.write(
"""------ Molecule Types -------------------------------------------------
2                    # number of sorbate types

{0}               # sorbate
{0}.mol           # sorbate coordinates file\n\n""".format(species))
        file.write(
"""{0}                # sorbate
{0}.mol             # sorbate coordinates file
------ Simulation Cell Information ------------------------------------
{0}                # Fundamental cell file
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
{1} #number of sorbates\n""".format(framework, len(species) if isinstance(species, list) else 1))
        if isinstance(species, list):
            for i in species:
                file.write(
"""{0} #sorbate name\n""".format(i))
        else:
            file.write("""{0} #sorbate name\n""".format(species))
        file.write(
"""------ GCMC Information -----------------------------------------------
1                 # No. of iterations
{0}              # temperature
Ideal Parameters   # Tag for the equation of state (NULL = Ideal Gas)
{1}                  # No. of simulation points
5000                # Block size for statistics
{2}                  # no. of sorbates\n""".format(T, n, len(list(species))))
        if isinstance(species, list):
            for i in species:
                file.write(
"""                 -------------------------
{0}            # Sorbate Name
{1}           #  pressure
Null               # sitemap filename (Null = no sitemap)
5                  # no of gcmc movetypes
1.0, 1.0, 1.0, 1.0, 1.0      # move type weights
RINSERT                   # type of move.1
RDELETE                   # type of move.2
RTRANSLATE                # type of move.4
0.2, 1                    # Delta Translate, adjust delta option (0=NO, 1=YES)
RROTATE
0.2, 1
IDFLIP
{2}\n""".format(i, 'pressure.dat' if pressure == 'file' else pressure, ', '.join(str(species.index(x)+1) for x in species if x != i)))
        else:
            file.write(
"""                 -------------------------
{0}            # Sorbate Name
{1}           #  pressure
Null               # sitemap filename (Null = no sitemap)
5                  # no of gcmc movetypes
1.0, 1.0, 1.0, 1.0, 1.0      # move type weights
RINSERT                   # type of move.1
RDELETE                   # type of move.2
RTRANSLATE                # type of move.4
0.2, 1                    # Delta Translate, adjust delta option (0=NO, 1=YES)
RROTATE
0.2, 1\n""".format(i, 'pressure.dat' if pressure == 'file' else pressure))
        if isinstance(species, list):
            file.write(
"""------ Configuration Initialization -------------------------------------
{0} #sorbate
{1}\n""".format(species[0], 'GCMC NULL' if Restart == None else 'RESTARTFILE {0}'.format(Restart)))
            for i in range(1,len(species)):
                file.write("{0}\n{1}\n".format(species[i], 'GCMC NULL' if name == 'gcmc.ctr' else Restart))
        else:
            file.write(
"""------ Configuration Initialization -------------------------------------
{0} #sorbate
{1}\n""".format(species, 'GCMC NULL' if Restart == None else 'RESTARTFILE {0}'.format(Restart)))
        file.write(
"""{0}              # Sorbate_Type
FIXED NULL
--------  Main Datafile Information --------
Energy, position, pair_energy  # contents of datafile

------ Movie Information ----------------------------------------------
movie.xyz            # Movie filename (output is in xyz format)
0, {1}            # Starting step, ending step
1000                    # Steps between frames
No                   # Include zeolite in movie: Yes or No
1, 1, 1             # Number of unit cell repeates to dump in x, y, z directions
-----------------------------------------------------------------------""".format(framework, n))
    logger.debug("GCMC control file written!")

#Writes a .ctr file for your postprocessing. There;s orobably a better way than using music_post, but I'm not there yet
def PostControlChanger(Species, n, framework, dirout):
    with open("{0}post.ctr" .format(dirout), 'w') as file:
        file.write("""####This section is apparently required for working with any post code. Who knows why? not me!
#
#
------------------------------------------------------------
   ### Required section ######
-- Post Processor Information ------------
GCMC                            # Type of simulation GCMC, NVTMC , MD ....
./{1}.{0}.con                    # basename for config files
1, {2}                          # first and last file numbers
post.ctr.out                       # name for new ctrlfile that will regenerated
postoutput          # Base name for output files
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
20       # Number of blocks into which data should be divided for stats""".format('.'.join(i for i in list(Species)), framework, n))
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
#SBATCH account=free

# -- Set up the environment
module purge
module load group ce-molsim stack
module load music/std

#locate the experiment directory
cd {3} #Change this bit to get to your directory
cd {2}

#set the paths for atoms, molecules, pmaps, and emaps WARNING softcoded in umbrella, check it's right!
export ATOMSDIR=../../../../atoms
export MOLSDIR=../../../../molecules
export PMAPDIR=../../../../maps
export EMAPDIR=../../../../maps

#create symbolic links to your interactions files and pressure file
ln -s ../atom_atom_file atom_atom_file
ln -s ../intramolecular_file intramolecular_file
ln -s ../sorb_sorb_file sorb_sorb_file
ln -s ../pressure.dat pressure.dat

# -- Run
music_gcmc gcmc.ctr > logfile.gcmc #runs your main simulation
music_post post.ctr > logfile.post #runs your postprocessing

""".format(Species, T, dirout, parentdir))
#        for i, value in enumerate(isotherm): #for each isotherm point you're simulating
#            file.write('music_gcmc {0}kpa_restart.ctr > {1}_restart.logfile\nmv finalconfig.xyz {3:02d}.{2}.{0}kpa.xyz\n'.format(value, int(i)+1, framework, int(iteration))) #set up a simulation to generate a new xyz file and move it to be names after the pressure point 
    os.chmod("{0}run.gcmc" .format(dirout), 0o777) #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    logger.debug("Run file written!")


def mixisotherm(species0, species1, p0, n=10, directory='./'):
    iso0 = []
    iso1 = []
    for i in range (0, n):
        iso0.append(round(float(p0)*float(1.-float(i)/(n-1)), 3))
        iso1.append(round(float(p0)*float(float(i)/(n-1)), 3))
    print('{0} isotherm: {1}\n{2} isotherm: {3}'.format(species0, iso0, species1, iso1))
    with open('{0}pressure.dat'.format(directory), 'w') as f:
        f.write('{0}\n'.format(species0))
        f.write('{0}\n'.format(n))
        f.write('{0}\n'.format(', '.join(str(thing) for thing in iso0)))
        f.write('{0}\n'.format(species1))
        f.write('{0}\n'.format(n))
        f.write('{0}'.format(', '.join(str(thing) for thing in iso1)))
    print('Exchange isotherm file written!')

    
##############################
##############################
##############################
#print('Warming up...')
#directorymaker('{0}/'.format(xptpath))
#for handler in logging.root.handlers[:]:
#    logging.root.removeHandler(handler)


#logger = logging.getLogger()
#fh = logging.FileHandler(filename='{0}/umbrella.log'.format(xptpath))
#fh.setLevel(logging.NOTSET)
#sh = logging.StreamHandler()
#sh.setLevel(logging.NOTSET)
#for handler in logging.root.handlers[:]:
#    print(handler)
#logger.addHandler(fh)
#logger.addHandler(sh)

#This bit looks for your molecule files and makes one if required. WARNING! It really likes to overwrite things. I've tried to change molwrite to 'x' not 'w', we'll see if that works

logger.info("Checking mol files...")
directorymaker('./molecules/')
speciespath = Path('./molecules/{0}.mol'.format(species))
if speciespath.exists(): #does your sorbent exist?
    logger.info('{0} already exists as a molecule file.' .format(species))
else:
    logger.warning("{0} doesn't exist as a molecule file at all! I'll need to make a new molecule file for {0}" .format(species))
    pdbpath = Path('./{0}.pdb'.format(species))
    if pdbpath.exists(): #do you have a .pdb of your sorbent in this directory?
        atomlist, connlist = dataextract(species) #read info from the pdb if so
        bondinfo = connectiontypeswrite(atomlist, connlist) #calculate bond information form that info
        logger.info(atomlist)
        molwrite(species, atomlist, bondinfo, connlist, './molecules/', forcefield) #write a new file if it doesn't exist already
    else:
        logger.warning("I'd need a .pdb file to get started though!")
        errormessages.append("""You don't have a sorbent molecule and I couldn't make one! You need to make one before running.""")


speciespath1 = Path('./molecules/{0}.mol'.format(species1))
if speciespath1.exists(): #does your sorbent exist?
    logger.info('{0} already exists as a molecule file.' .format(species1))
else:
    logger.warning("{0} doesn't exist as a molecule file at all! I'll need to make a new molecule file for {0}" .format(species1))
    pdbpath = Path('./{0}.pdb'.format(species1))
    if pdbpath.exists(): #do you have a .pdb of your sorbent in this directory?
        atomlist, connlist = dataextract(species1) #read info from the pdb if so
        bondinfo = connectiontypeswrite(atomlist, connlist) #calculate bond information form that info
        logger.info(atomlist)
        molwrite(species1, atomlist, bondinfo, connlist, './molecules/', forcefield1) #write a new file if it doesn't exist already
    else:
        logger.warning("I'd need a .pdb file to get started though!")
        errormessages.append("""You don't have a sorbent molecule and I couldn't make one! You need to make one before running.""")

#Now I'm doing the same for your sorbate
logger.info("Checking mol files...")
directorymaker('./molecules/')
if Path('./molecules/{0}.mol'.format(framework)).exists():
    logger.info('{0} already exists as a molecule file.' .format(framework))
else:
    logger.info("{0} doesn't exist as a molecule file at all! You'll need to make a new molecule file for {0}" .format(framework)) #But I'm not clever enough to make a molfile of your framework, sorry!
    errormessages.append("""You don't have a framework molecule! You need to make one before running.""")
    #should I add a sys.exit statement to this?

logger.info('------------------------------Finished looking for molecules--------------------')
        
#This bit now looks for all of your map files
logger.info("Checking for maps of elements within {0}" .format(speciespath))
el_list = set() #creates a list of elements in your sorbent, excluding multiples
if speciespath.exists():
    with speciespath.open() as file:
        for line in file:
            if len(line.split()) == 8: #should be 8 (index x y z name charge ? ?), but you might have added comments
                logger.info(line.split()[4])
                el_list.add(line.split()[4])
    logger.info("Atom list found is: " + str(el_list))
else:
    logger.warning("That's weird, I can't find your mol file. Something must have messed up!")
    errormessages.append("""I lost your moelcule file while looking for your maps. I have no idea what happened, but this must be bad.""")
    #should I add a sys.exit statement to this?

el_list1 = set()
if speciespath1.exists():
    with speciespath1.open() as file:
        for line in file:
            if len(line.split()) == 8: #should be 8 (index x y z name charge ? ?), but you might have added comments
                logger.info(line.split()[4])
                el_list1.add(line.split()[4])
    logger.info("Atom list found is: " + str(el_list))
else:
    logger.warning("That's weird, I can't find your mol file. Something must have messed up!")
    errormessages.append("""I lost your moelcule file while looking for your maps. I have no idea what happened, but this must be bad.""")
    #should I add a sys.exit statement to this?

sorb_el_list = [] #i want to differentiate between framework elements and sorbent elements here, so i rename them as X_sorb
for i in el_list:
    sorb_el_list.append("{0}_{1}" .format(i.split("_")[0], forcefield)) #Sometimes I'll call them X_something esle to denote their forcefield for instance. This hould cut it back to just X
for i in el_list1:
    sorb_el_list.append("{0}_{1}" .format(i.split("_")[0], forcefield1)) #Sometimes I'll call them X_something esle to denote their forcefield for instance. This hould cut it back to just X
logger.info(sorb_el_list)

for i in sorb_el_list: #Now I have a list, lets look for those maps!
    logger.info(i)
    directorymaker('./maps/')
    mappath = Path('./maps/{0}.{1}.p'.format(framework, i))
    if mappath.exists():#(if you have a map there, and it's the right framework)
        logger.info("Interactions pmap found between fluid atom {0} and a framework {1}" .format(i, framework))
    else:
        errormessages.append("""I couldn't find the pmaps you probably needed, so I tried to make some setup files for you.""")
        logger.warning("{0} map on {1} doesn't exist! You'll have to make a new one!" .format(i, framework))
        if Path('./molecules/{0}.mol'.format(i)).exists(): #do you have a map molecule file?
            logger.warning("{0} already exists as a molecule file." .format(species))
        else:
            logger.warning("No map molecule found for {0}, I'll have to write a new one!" .format(i, 'mol'))
            molwrite("{0}".format(i), [i], None, None, "./molecules/") #makes a new one if needed
            errormessages.append("""I had to write new molecule files for your map generation.""")
        directorymaker('./mapgen/{0}.{1}/' .format(species, species1)) #makes a map generation directory for your framework
        Intsetup(sorb_el_list, sorb_el_list, sorb_el_list, framework, './mapgen/{0}.{1}/' .format(species, species1), False, False) #sets up you simulation interactions
        MapMakeRunWriter(i, sorb_el_list, parentdir, './mapgen/{0}.{1}/' .format(species, species1)) #makes a run file
        mapctrlfilewriter(i, framework, './mapgen/{0}.{1}/' .format(species, species1)) 
        logger.info("Your map is ready to be generated! Just go to ./mapgen/{0}.{1}/ and run run.mapmaker" .format(species, species1))

#This bit defines your pressure information
logger.info('Finished looking for maps--------------------')
logger.info("""Now I'll make an isotherm for you""")
satP = pSat(species, T)
satP1 = pSat(species1, T)
satPuse = min(satP, satP1)
logger.info('Your sorbents have the following saturation pressures:\n {0}: {1}\n{2} {3}\nTherefore I\'ll use {4} as my total fugacity.'.format(species, satP, species1, satP1, satPuse))
#istm = isothermcalculator(satP, n, minrelpress, maxrelpress)


logger.info('Finished calculating your isotherm--------------------')
logger.info("""Now I'll set up your simulation documents for you""")
#This bit now sets up your experiment. WARNING! It overwrites all previous experiment files in the directory chosen (but significantly doesn't overwrite any results files)        
directorymaker('{0}/'.format(xptpath)) #makes your general experiment directory for this isotherm
directorymaker('./experiments/finalconfigs/{0}.{2}/{1}/'.format(species, T, species1)) #makes your 16 subdirectories for taskfarming
logger.info("Working in directory {0}:".format(xptpath))
mol_mix = [species, species1]
spec0 = []
spec0.append(i for i in el_list)
spec1 = []
spec1.append(i for i in el_list1)
el_mix = []
el_mix.append(spec0)
el_mix.append(spec1)
IsothermExtractMover(mol_mix, T, framework, '{0}/'.format(xptpath), n)
TaskfarmRunWriter(mol_mix, T, framework, parentdir, '{0}/'.format(xptpath))                                                          #Writes the taskfarmer run file
Intsetup(mol_mix, sorb_el_list, el_mix, framework, '{0}/'.format(xptpath), True, True)
mixisotherm(species, species1, satPuse, n, '{0}/'.format(xptpath))
#PressureFileWriter(species, T, satP, istm, '{0}/'.format(xptpath, directory))  #writes your isotherm
for directory in range (1,17):
    directorymaker('{0}/{1:02d}/'.format(xptpath, directory))
    logger.info("Working in directory {0}:".format('.{0}/{1:02d}/'.format(xptpath, directory)))
    GcmcControlChanger(mol_mix, sorb_el_list, T, n, framework, '{0}/{1:02d}/'.format(xptpath, directory), '1000000')                            #puts in the gcmc control file
    PostControlChanger(mol_mix, n, framework, '{0}/{1:02d}/'.format(xptpath, directory))                            #puts in the post control file
    GcmcRunWriter(species, T, framework, parentdir, '{0}/{1:02d}/'.format(xptpath, directory), None, directory) #puts in the runfile

logger.info('Finished writing your simulation files--------------------')
logger.warning("""Now I'll wrap up.""")
if len(errormessages) == 0:
    logger.warning('This script ran without generating any error messages')
else:
    logger.warning('The following errors were generated during the script running.')
    for i in errormessages:
        logger.warning(i + '\n')
logger.warning('Your simulations are ready to go in directory {0}!'.format(xptpath))
