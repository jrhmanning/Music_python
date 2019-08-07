######This python script take in a list of your elements and some forcefield information, then spits out a lot of control files to make maps of them all
#It requires the following variables: 
#species - the name of your sorbent species. Since we're making a map, this is actually one of your sorbent elements
#framework - the name of your MOF (or other porous material
#################
from math import sqrt, log
import logging
from pathlib import Path
import datetime
now = datetime.datetime.now()
import random
#################
####This section sorts out your messages from this script to the console and a log file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler('{0}/sorb_sorb_maker.log'.format('.'))
handler.setLevel(logging.DEBUG)
handler2 = logging.StreamHandler()
handler2.setLevel(logging.INFO)
logger.addHandler(handler)
logger.addHandler(handler2)
logger.debug('-----------------------------------')
logger.debug('Generating new files on {0}-{1}-{2} at {3}:{4}'.format(now.year, now.month, now.day, now.hour, now.minute))
logger.debug('-----------------------------------')
##################
species = ['XXXX']
sorb_el_list = ['''### your species list here###''']
framework = 'YYYY'
##################
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
###########
for i in sorb_el_list:
	mapctrlfilewriter(i, framework, "./mapgen/{0}".format(species[-1]))