######This python script take in a list of your elements and some forcefield information, then spits out a control file to manage your gcmc simulations
#It requires the following variables: 
#species - the name of your sorbent species
#elements - to tell your control file how many and what elements to expect
#n - your number of smiulation points
#T - your simulaiton temperature
#framework - the name of your MOF (or other porous material
#iterations = the number of iterations you're using
#Restart = if you want to restart your simulation from a certian point, or to get a snapshot
#name - useful for autogenerating control files for snapshots
#Pressure - useful if you want to make a snapshot
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
species = ['DMFYang']
sorb_el_list = ['AldH_s', 'O_s', 'Carb_s', 'C_s', 'N_s', 'H_s']
MOF_el_list = ['Carbon', 'Oxygen', 'Zinc', 'Hydrogen']
framework = 'IRMOF1'
T = 298
n = 9
##################
#Writes a .ctr file for your gcmc simulations (production and for getting your final configuarions at the end). This currently cannot cope with multiple sorbents at once, creates a 2x2x2 cell, and has bias insert atom softcoded in
def GcmcControlChanger(species, elements, T, n, framework, dirout, iterations = '750000', Restart = None, name = 'gcmc.ctr', pressure = 'file'): # where the function definition already gives a value, this means it's optional and the function defaults to the one here
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
{3}                    # number of atomic types            \n\n""".format(species, framework, x, n_species, iterations, 1 if Restart == None else 21))
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
          ------------------------- #Repeat this section for each sorbate
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
#######################
for directory in range(1, 17):
    GcmcControlChanger(species[-1], sorb_el_list, T, n, framework, './{0:02d}/'.format(directory), '1000000', 'RESTARTFILE {0}.{1}.res.{2}'.format(framework, species[-1], 20)) #for production runs
#for i, value in enumerate(istm): #now i make the extra .ctr files for yoru subdirectories
#    GcmcControlChanger(species[-1], sorb_el_list, T, 1, framework, '{0}/{1:02d}/'.format(xptpath, directory), '1', 'RESTARTFILE {0}.{1}.res.{2}'.format(framework, species[-1], i+1), '{0}kpa_restart.ctr'.format(value), value) #makes your control files for all of your pressure points

