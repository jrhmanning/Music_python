######This python script take in a list of your species and temperature, then spits out a bash runscript to run your GCMC simulations
#It requires the following variables: 
#species - the name of your sorbent species
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
T = 298
framework = 'YYYY'
parentdir = Path('/home/r/jrhm21/scratch/03_music_chloroform_forcefield/')
istm = ['''### your pressure points here###''']
##################
#Writes a bashscript for actually running your simulations, compatible with SLURM
def GcmcRunWriter(Species, T, framework, parentdir, dirout, isotherm):
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

#set the paths for atoms, molecules, pmaps, and emaps WARNING softcoded in umbrella, check it's right!
export ATOMSDIR=../../../../atoms
export MOLSDIR=../../../../molecules
export PMAPDIR=../../../../maps/{0}
export EMAPDIR=../../../../maps


# -- Run
music_gcmc gcmc.ctr > logfile.gcmc #runs your main simulation
music_post post.ctr > logfile.post #runs your postprocessing

""".format(Species, T, dirout, parentdir))
        for i, value in enumerate(isotherm): #for each isotherm point you're simulating
            file.write('music_gcmc {0}kpa_restart.ctr > {1}_restart.logfile\nmv finalconfig.xyz {2}.{0}kpa.xyz\n'.format(value, int(i)+1, framework)) #set up a simulation to generate a new xyz file and move it to be names after the pressure point 
    os.chmod("{0}run.gcmc" .format(dirout), 0o777) #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    logger.debug("Run file written!")
####################
GcmcRunWriter(species[-1], T, framework, parentdir, './experiments/{0}/{1}/'.format(species[-1], T), istm)