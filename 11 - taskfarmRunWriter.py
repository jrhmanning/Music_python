######This python script take in a list of your species and temperature, then spits out a bash runscript to taskfarm 16 simulations at once, for more efficient computation
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
##################
#Writes a bashscript for taskfarming yuor simulations, compatible with SLURM
def TaskfarmRunWriter(Species, T, framework, parentdir, dirout):
    with open("{0}run.taskfarmer" .format(dirout), 'w') as file:
        file.write("""#!/usr/bin/env bash
#SBATCH --job-name={0}.{2}.{1}k
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --time=06:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jrhm21@bath.ac.uk #Change this to your email!
#SBATCH account=free

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
############
TaskfarmRunWriter(species[-1], T, framework, parentdir, './experiments/{0}/{1}/'.format(species[-1], T))