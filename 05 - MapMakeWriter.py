######This python script take in a list of your species and elements, then spits out a bash runscript to make maps for them all in conjunction with control files
#It requires the following variables: 
#species - the name of your sorbent species
#elements - a list of the elemetns involved in your sorbent species, for maps
#framework - the name of your MOF (or other porous material
#parentdir = a path to navigate you to your testing directory
#################
import os
from math import sqrt, log
import logging
from pathlib import Path
import datetime
now = datetime.datetime.now()
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
framework = 'YYYY'
sorb_el_list = ['''### your species list here###''']
parentdir = Path('/home/r/jrhm21/scratch/03_music_chloroform_forcefield/')
#################
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
export PMAPDIR=../../maps/{0}
export EMAPDIR=../../maps
\n\n# --Run\n""" .format(species, dxout.split(".")[-1], parentdir))
        for count, i in enumerate(elements, 1):
            file.write("music_mapmaker makemap_{0}.ctr > logfile{1}\n" .format(i, count, parentdir))
    os.chmod("{0}run.mapmaker" .format(dxout), 0o777) #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    logger.debug("Runfile written!")   
###########
MapMakeRunWriter(species[-1], sorb_el_list, '{0}/'.format(parentdir), "./mapgen/{0}".format(species[-1]))