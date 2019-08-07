
######This python script takes in pdb files and spits out a .mol file for music
#It requires the following variables: 
# species - the molecule you're considering (also used to read [species].pdb and create [species].mol)
# sourcedir - the directory you're reading [species].pdb from
# (optional) - forcefield to load coulobic information from a forcefield dictionary directly into your .mol file
######
#While working, it creates the following variables:
#atoms - a list of atoms in your species
#bonds - a list of bonds in your species
#bondlengths - a dictionary of the bonds in your molecule, with their respective lengths
import sys
from math import sqrt, log
import logging
import datetime
now = datetime.datetime.now()

element_names = {'H': 'Hydrogen', #In my PDB files, they're as symbols, and in my music files, they're names. This dictionary allows for interconversion
                 'C': 'Carbon',
                 'Cl': 'Chlorine',
                 'O': 'Oxygen',
                 'N': 'Nitrogen',
                 'F': 'Fluorine',
                 'CH': 'CH',
                 'ArC':'ArC',
                 'CH2': 'CH2',
                 'Me': 'Me',
                 'CO': 'Carbonyl',
                 'HO': 'AldH',
                 'X': 'COMN1',
                 'Y': 'COMN2'
}

####This section sorts out your messages from this script to the console and a log file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler('{0}/molfilemaker.log'.format('.'))
handler.setLevel(logging.DEBUG)
handler2 = logging.StreamHandler()
handler2.setLevel(logging.INFO)
logger.addHandler(handler)
logger.addHandler(handler2)
logger.debug('-----------------------------------')
logger.debug('Generating new files on {0}-{1}-{2} at {3}:{4}'.format(now.year, now.month, now.day, now.hour, now.minute))
logger.debug('-----------------------------------')
######

species = 'EG_TRAPPE'

######
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
            thing[4] = '{0}_{1}'.format(element_names[thing[4]],species) #replaces the atom symbol from above with the element name
    return atoms, connections

#This sorts the connections detected from apdb into a list for a molfile
def bondlengthswrite(atoms, connections):
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
            print(pair)
            print(bondlengths[pair])
            for i in bondlengths[pair]:
                if float(bondlengths[pair][0])-float(i) > 0.01: #if your bond lengths vary more than this, throw a fit.
                    logger.warning('ERROR! Tolerance exceeded for bonds between atoms of types %s.' % (str(pair)))
                    sys.exit()
            bondlengths[pair] = round(sum(bondlengths[pair])/len(bondlengths[pair]), 3) 
    return bondlengths #returns the bondlength dictionary of 'i j pair': length

#This writes a molecule file from the informaiton from dataextract and bondlengths. It can overwrite your existing files though, beware!
def molwrite(species, atoms, bondlengths, connections, sourcedir = './', forcefield = None):
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
    if bondlengths:
        f.write("Connect_Info: listed\n") #declares it's listing the number of bond dypes as well as the exact bond pairs
        f.write("{0} # number of unique bond types in your molecule\n" .format(len(bondlengths)))
        for i in bondlengths:
            x = str(i) #makes a string of you 'i j' pair
            x = x.strip("('" "')") #cuts out the quites from the string
            x = x.split("', '") #and the comma from between them
            logger.debug(x) #now it should just be 2 atom type names
            f.write("{0} {1} {2}\n" .format(x[0], x[1], bondlengths[i])) #now it writes atom type 1, atom type 2, length
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
######

atoms, connections = dataextract(species)
bondlengths = bondlengthswrite(atoms, connections)
molwrite(species, atoms, bondlengths, connections)
