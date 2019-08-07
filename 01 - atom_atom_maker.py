######This python script takes in a list of your elements and some forcefield information, then spits out an atom_atom_file for music
#It requires the following variables: 
#elements - a list of the elemetns involved in your sorbent species
#MOF_el_list - a list of your framework atom types
#IntParams - a dictionary of your IntParams
#forcefield - the name of the forcefield, which is a suffix to your species atoms in IntParams
####
# If you're wanting to generate your element list from [species].mol, you'll be able to use this too.

#################
from math import sqrt, log
import logging
from pathlib import Path
import datetime
now = datetime.datetime.now()
#################

forcefield = 'DMFAA' #for example


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
    "Carbon_Kamath":(3.41, 68.94, -0.235), #Chlorfoorm parametrs taken from 10.1021/jp0535238
    "Hydrogen_Kamath":(2.81, 10.06, 0.355),
    "Chlorine_Kamath":(3.45, 138.58, -0.04),
    "Nitrogen_TRAPPE":(3.310, 36.0, -0.482),
    "COMN_TRAPPE":(0, 0, 0.964),
    "Methynyl_OPLS": (3.8, 40.26, 0.42),
    "Chlorine_OPLS": (3.47, 150.98, -0.14),
    "Hydrogen_CDP": (2.81, 10.06, -0.0551),
    "Carbon_CDP": (3.41, 68.94, 0.5609),
    "Chlorine_CDP": (3.45, 138.58, -0.1686),
    "Chlorine_UFF": (3.947, 114.128, 0),
    "Carbon_UFF": (3.851, 52.79, 0),
    "Hydrogen_UFF": (2.886, 22.12, 0),
    
    "Carbon_DMFAA":(3.50, 33.20,-0.110), #from 10.1016/j.molliq.2015.03.004
    "Hydrogen_DMFAA":(2.50, 15.16,0.060),
    "Nitrogen_DMFAA":(3.25, 85.52,0.040),
    "Oxygen_DMFAA":(2.96, 105.61,-0.680),
    "Carbonyl_DMFAA":(3.75, 52.80,0.5),
    "AldH_DMFAA":(2.5, 15.16,0.00),
    
    "Carbon_DMFYang":(3.50,33.20,-0.110), ##from 10.1016/j.molliq.2015.03.004
    "Hydrogen_DMFYang":(2.5,15.16,0.060),
    "Nitrogen_DMFYang":(3.25,85.52,-0.14),
    "Oxygen_DMFYang":(2.96,105.61,-0.5),
    "Carbonyl_DMFYang":(3.75,52.80,0.5),
    "AldH_DMFYang":(2.5,15.16,0),
    
    "Carbon_DMFCaleman":(3.5,33.20,-0.110), #from 10.1016/j.molliq.2015.03.004
    "Hydrogen_DMFCaleman":(2.5,15.16,0.06),
    "Nitrogen_DMFCaleman":(3.25,85.52,-0.14),
    "Oxygen_DMFCaleman":(2.96,105.61,-0.5),
    "Carbonyl_DMFCaleman":(3.75,52.80,0.5),
    "AldH_DMFCaleman":(2.5,15.16,0)
    }

MOF_el_list = ["Carbon", "Hydrogen", "Oxygen", "Zinc"] #This is a list of elements in your MOF, needed for the control file atom types section

####This section sorts out your messages from this script to the console and a log file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler('{0}/atom_atom_maker.log'.format('.'))
handler.setLevel(logging.DEBUG)
handler2 = logging.StreamHandler()
handler2.setLevel(logging.INFO)
logger.addHandler(handler)
logger.addHandler(handler2)
logger.debug('-----------------------------------')
logger.debug('Generating new files on {0}-{1}-{2} at {3}:{4}'.format(now.year, now.month, now.day, now.hour, now.minute))
logger.debug('-----------------------------------')

################ This section creates a list of your sorbent atom types from a .mol file, 
speciespath = Path('./DMFAA.mol')
el_list = set()
if speciespath.exists():
    with speciespath.open() as file:
        for line in file:
            if len(line.split()) == 8: #should be 8 (index x y z name charge ? ?), but you might have added comments
                logger.info(line.split()[4])
                el_list.add(line.split()[4])
    logger.info("Atom list found is: " + str(el_list))

sorb_el_list = [] #i want to differentiate between framework elements and sorbent elements here, so I rename them as X_sorb
for i in el_list:
    sorb_el_list.append("{0}_{1}" .format(i.split("_")[0], forcefield)) #Sometimes I'll call them X_something esle to denote their forcefield for instance. This hould cut it back to just X
logger.info(sorb_el_list)


################

#This writes your atom-atom forcefield information based on the above dictionary at the start of the document
def AtmAtmMover(elements, forcefield, MOF_el_list = MOF_el_list, IntParams = IntParams, dxout = "./mapgen/test/"):
    tot_el_list = []
    tot_el_list.append(MOF_el_list) #takes the global list from the preamble to getyour mof molecule atom types
    tot_el_list.append(elements) #reads what elements are in yoru sorbent
    hicut = 18 #sets your hicut in A
    with open("{0}atom_atom_file" .format(dxout), "w") as file:
        file.write("""#This is an autogenereated interactions file for music, based on the python script Umbrella-v2. It's probably best to look this over before running your sims\n#Lennard-Jones interactions\n""")
        for i in range (0, len(MOF_el_list)):
            file.write("""{0} {0} LJ SIG@{1} EPS@{2} HICUT@{3}\n""".format(MOF_el_list[i], IntParams[MOF_el_list[i]][0], IntParams[MOF_el_list[i]][1], hicut)) #LJ self-parameters for the MOF-MOF interactions
            for j in range (i+1, len(MOF_el_list)):
                file.write("""{0} {1} LJ OFF\n""".format(MOF_el_list[i], MOF_el_list[j])) #LJ parameters for MOF 'i j' pairs, which are generally off
            for k in elements:
                file.write("""{0} {1} LJ SIG@{3} EPS@{4} HICUT@{2}\n""".format(MOF_el_list[i], k, hicut, str((IntParams[MOF_el_list[i]][0]+IntParams['{0}_{1}'.format(k.split('_')[0], forcefield)][0])/2), str(sqrt(IntParams[MOF_el_list[i]][1]*IntParams['{0}_{1}'.format(k.split('_')[0], forcefield)][1])))) #MOF-fluid forcefield paramters
            file.write("\n")
        for i in range (0, len(elements)):
            raw_el = str(elements[i].split('_')[0] + '_' + str(forcefield)) #This line says which interaction types your fluid will take from the dictionary at the start of the script
            file.write("""{0} {0} LJ SIG@{1} EPS@{2} HICUT@{3}\n""".format(elements[i], IntParams[raw_el][0], IntParams[raw_el][1], hicut)) #LJ self-parameters for the fluid-fluid interactions
            for j in range (i+1, len(elements)):
                logger.debug('i = '+ str(i)+ ', j= '+str(j) )
                file.write("""{0} {1} LJ SIG@{3} EPS@{4} HICUT@{2}\n""".format(elements[i], elements[j], hicut, str((IntParams['{0}_{1}'.format(elements[i].split('_')[0], forcefield)][0]+IntParams['{0}_{1}'.format(elements[j].split('_')[0], forcefield)][0])/2), str(sqrt(IntParams['{0}_{1}'.format(elements[i].split('_')[0], forcefield)][1]*IntParams['{0}_{1}'.format(elements[j].split('_')[0], forcefield)][1])))) #LJ parameters for the fluid-fluid 'i j' interactions
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
            file.write("""{0} {0} WFCOUL  HICUT@{1}\n""".format(elements[i], hicut)) #fluid fluid self-interactions - use wolf cpoulombic. partial charge values are stored in the molecule files - see molfilewriter
            for j in range (i+1, len(elements)):
                file.write("""{0} {1} WFCOUL HICUT@{2}\n""".format(elements[i], elements[j], hicut)) #fluid fluid 'i j' interactions - use wolf cpoulombic. partial charge values are stored in the molecule files - see molfilewriter
            file.write("\n")
    logger.debug("Atom-Atom interaction file written!")
###################

AtmAtmMover(sorb_el_list, forcefield, MOF_el_list, IntParams, "./")
