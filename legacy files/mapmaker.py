#This is a python script that tales a molecule structure from a .pdb file and generates a map for each unique atom type
import os
import errno
import random
import sys
from math import sqrt

element_names = {'H': 'Hydrogen',
                 'C': 'Carbon',
                 'Cl': 'Chlorine',
                 'O': 'Oxygen',
                 'N': 'Nitrogen',
                 'F': 'Fluorine'
}
    

def dataextract(pdb = "Chloroform.pdb", dxin = "./"):    #I stole this from the molfile generator and added the last few lines creating an element set, so i could potentially combine them eventually?
    molecule = pdb.split(".")[0]
    atoms = []
    connections = []
    f = open("{0}{1}" .format(dxin, pdb), "r")
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
    elements = set()                   #not included in molfile generator
    for thing in atoms:
        elements.add(thing[4])
    print("pdb file read!")
    return atoms, connections, elements, molecule

def MapMolWrite(elements, version, dxout = "./molecules/"):
    for i in elements:
        with open('./experiments/repository/Element_pmap_uncharge.mol.template', 'r') as f:
            filedata = f.read()
        filedata = filedata.replace('*your_element_name*', i)
        filedata = filedata.replace('*your_version_number*', version)
        with open("{0}{1}_v{2}.mol" .format(dxout, i, version), "w") as f:
            f.write(filedata)
        print("Molecule file {0} written!" .format(i))

def AtmAtmMover(version, dxout = "./mapgen/test/"):
    with open('./experiments/repository/atom_atom_file_v{0}.template' .format(version), 'r') as f:
        filedata = f.read()
    with open("{0}atom_atom_file" .format(dxout), "w") as f:
            f.write(filedata)
    print("Atom-Atom interaction file (version{0}) written!" .format(version))

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


def SorbSorbWriter(Species, elements, version, dxout = "./mapgen/"):
    file = "sorb_sorb_file"
    with open("{0}{1}" .format(dxout, file), 'w') as file:
        file.write("IRMOF1 IRMOF1 NCOUL OFF\nIRMOF1 IRMOF1 COUL OFF\n\n") #MOF-MOF interactions
        file.write("{0}_v{1} IRMOF1 NCOUL MAP@IRMOF1 FAST ".format(Species, version))
        for i in elements:
            file.write("{0}_sorb@PMAP@IRMOF1.{0}_v{1}.pmap " .format(i, version))
        file.write("\n{0}_v{1} IRMOF1 COUL OFF\n\n" .format(Species, version))
        file.write("{0}_v{1} {0}_v{1} NCOUL BASIC LJ FAST\n{0}_v{1} {0}_v{1} COUL OFF\n\n" .format(Species, version))
	for i in elements:
            file.write("{0}_v{1} IRMOF1 NCOUL BASIC LJ FAST\n{0}_v{1} IRMOF1 COUL OFF\n\n" .format(i, version))
            file.write("{0}_v{1} {0}_v{1} NCOUL BASIC LJ FAST\n{0}_v{1} {0}_v{1} COUL OFF\n\n" .format(i, version))
    print("Sorb-Sorb file written!")    

def IntraWriter(species, elements, version, dxout = "./mapgen/"):
    with open("{0}intramolecular_file" .format(dxout), "w") as file:
        file.write("####This intramolecular interaction file is workable for amking a map of {0} in a 2x2x2 IRMOF-1 framework\n" .format(species))
        for i in elements:
            file.write("Intra: {0}_v{1}\n" .format(i, version))
        file.write("Intra: IRMOF1")
    print("Intra file written!")

def MapMakeRunWriter(Species, elements, dxout = "./mapgen/"):
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
cd 01_music_IRMOF_solvent{1}\n
export ATOMSDIR=../../atoms
export MOLSDIR=../../molecules
export PMAPDIR=../../maps
export EMAPDIR=../../maps
\n\n# --Run\n""" .format(Species, dxout.split(".")[1]))
        for count, i in enumerate(elements, 1):
            file.write("music_mapmaker makemap_{0}_v{1}.ctr > logfile{2}\n" .format(i, version, count))
        print("Runfile written!")
    
def controlfilewriter(species, elements, version, dxout = "./mapgen/"):
    x = random.randint(0,99999)
    for i in elements:
        with open('./experiments/repository/makemap.ctr.template', 'r') as f:
            filedata = f.read()  
        filedata = filedata.replace("*your_Species_name*", species)
        filedata = filedata.replace("*your_version_here*", version)
        filedata = filedata.replace("*your_seed_here*", "{0:05d}" .format(x))
        filedata = filedata.replace("*your_element_name*", i)
        with open("{0}makemap_{1}_v{2}.ctr" .format(dxout, i, version), "w") as file:
            file.write(filedata)
    print("Control files written!")

def connectiontypeswrite(atoms, connections):
    atomnames = []
    pairset = set()
    pairlengths = {}
    for i in atoms:
        atomnames.append(i[4])
    print(atomnames)
    #atomstypes = dict.fromkeys(atomnames)
    indices = []
    for i in connections:
        for j in i:
            if j != i[0]:
                pairset.add(tuple(sorted([atomnames[int(i[0])-1], atomnames[int(j)-1]])))
		pairdict = dict()
		for k in pairset:
			pairdict[k] = []
    for i in connections:
        for j in i:
            key = tuple(sorted([atomnames[int(i[0])-1], atomnames[int(j)-1]]))
            if j != i[0]:
                dist = 0
                dist = sqrt((float(atoms[int(i[0])-1][1])-float(atoms[int(j)-1][1]))**2+(float(atoms[int(i[0])-1][2])-float(atoms[int(j)-1][2]))**2+(float(atoms[int(i[0])-1][3])-float(atoms[int(j)-1][3]))**2)
                dist = round(dist, 3)   #rounds to nearest 0.001 of an angstrom. Is this good enough?
                pairdict[key].append(dist)
    for pair in pairdict:
        if len(pairdict[pair]) > 1:
            for i in pairdict[pair]:
                if float(pairdict[pair][0])-float(i) > 0.002:
                    print('ERROR! Tolerance exceeded for bonds between atoms of types {0}. Bond discrepancy measured at {1} Angstrom' .format(str(pair), str(round(float(pairdict[pair][0])-float(i), 4))))
                    #sys.exit()
            pairdict[pair] = round(sum(pairdict[pair])/len(pairdict[pair]), 3)
    return pairdict
    

def molwrite(species, version, atoms, connectiontypes, connections, dxout = "./molecules/"):
    f = open("{0}{1}_v{2}.mol" .format(dxout, species, version), "w")
    f.write("### Music molecule construction information, generated from {0} by the molfile generator \n\n" .format(species))
    f.write("#Basic molecule information \n")
    f.write("Molecule_Name: {0}_v{1} \n" .format(species, version))
    f.write("#Atom coordinate \nCoord_Info: listed cartesian rigid\n") #here i'm telling the moelcule file your moecule has listed atoms, cartesian coordinates, and rigid relationships between them
    f.write("{0} \n" .format(len(atoms))) #this prints the number of atoms to your file
    print(atoms)
    for i in atoms:
        f.write('{0} {1} {2} {3} {4}_sorb ' .format(i[0], i[1], i[2], i[3], i[4])) #prints the atoms
        f.write("0 0 0\n")# .format(charges[int(i[0])-1]))
    f.write("Connect_Info: listed\n") #declares it's listing the number of bond dypes as well as the exact bond pairs
    f.write("{0} # number of unique bond types in your molecule\n" .format(len(connectiontypes)))
    for i in connectiontypes:
        x = str(i)
        x = x.strip("('" "')")
        x = x.split("', '")
        f.write("{0}_sorb {1}_sorb {2}\n" .format(x[0], x[1], connectiontypes[i]))
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
    print("Molecule file {0} written!" .format(species))
	
###############################
random.seed()
atoms, connections, elements, mol_name = dataextract("THF.pdb")
version = "3"
dxout = "./mapgen/{0}_v{1}/" .format(mol_name, version)
directorymaker(dxout)
directorymaker("./molecules/")

bondtypes = connectiontypeswrite(atoms, connections)
#charges = [0, 0, 0, 0, 0]
molwrite(mol_name, version, atoms, bondtypes, connections)
MapMolWrite(elements, version)
AtmAtmMover(version, dxout)
SorbSorbWriter(mol_name, elements, version, dxout)
IntraWriter(mol_name, elements, version, dxout)
MapMakeRunWriter(mol_name, elements, dxout)
controlfilewriter(mol_name, elements, version, dxout)


