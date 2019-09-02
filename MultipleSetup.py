######This python script prepares a simulaiton for you assuming you've already got maps, molecules, atoms, and interactions files for everything in your simulation
###### It writes the appropriate fiels for yout o 
#It requires the following variables: 
#species - the name of your sorbent species. Since we're making a map, this is actually one of your sorbent elements
#T - your simulation temperature
#n - the number of isotherm points you're using
#minrelpress - the minimum pressure point you want to simulate, relative to the saturation pressure These are good for zeroing in on your pressure region
#maxrelpress - the maximum pressure point you want to simulate, relative to the saturation pressure
#################
import sys
import os
import re
import random
import numpy as np
from math import sqrt, log
from pathlib import Path
import logging
import datetime
import argparse

import musicpy.setup as setup
import musicpy.Antoine as Antoine
import musicpy.Forcefield as Forcefield
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
######your experiment variables
#### Directory stuff
username = 'jrhm21' #your username, for writing paths and slurm outputs
parentdir = Path('/home/r/{0}/scratch/03_music_chloroform_forcefield/'.format(username)) #The directory this python script is in
targetdir = Path('./experiments/{0}/{1}/'.format(species[-1], T)) #The directory you want to output to

#### simulation basics
species = ['DMF', 'MeOH_TRAPPE', 'EtOH_TRAPPE' 'EG_TRAPPE', 'Nitrogen_TRAPPE'] #Name of your sorbents, as a list
sorb_el_list = [
['AldH_s', 'O_s', 'Carb_s', 'C_s', 'N_s', 'H_s'],
['Me_s', 'O_s', 'H_s'],
['Me_s', 'CH2_s', 'O_s', 'H_s'],
['CH2_s', 'O_s', 'H_s'],
['N_s', 'COM_s']] #different elements in yoru sorbent, as a list
framework = 'IRMOF1' #name of your framework
T = 298 #simulation temperature

##### forcefield information
forcefield = ['DMFAA', 'TRAPPE', 'TRAPPE', 'TRAPPE', 'TRAPPE'] # lookup code for you interaction parameters in Forcefield.py 
hicut = 18 #atom-atom interaciton cutoff, in angstrom
coultype = 'Ewald' #fluid-fluid interactions calcualtion type
pmap = True #do you have fluid-framework pmaps
emap = True #do you have fluid-framework emaps?

##### adsorption information
iso_length = 9 #number of isotherm points
minrelpress = -1.5 #Minimum pressure to be considered, as a fration of saturation pressure (pMin = pSat*10^minrelpress) 
maxrelpress = -1 #maximum pressure to be considered, accoring to the above equation


#simulation information
n_iterations = '1000000' #number of MC steps int he main simulation
restart = None #if you want to write a restart control file, change this from None to 'RESTARTFILE {0}.{1}.res.{2}'.format(framework, species[-1], 20), where 20 is the restartfile you're going from
ctrl_file_name = 'gcmc.ctr' #name for your control file. Main one is gcmc.ctr, other ones can be called as you like
pressure = 'file' #Useful for restarts if you want the pressure to eb asingle value e.g. '10' (kPa)

######## Taskfarming information
Taskfarm_total = 15 
n_per_species = 3
n_species = len(species)

if Taskfarm_total != n_per_species*n_species:
    logger.warning('''###########################
	###########################
	You've got a different number of taskfarm jobs ({0}) to your total number of species* replicates per species ({1})!
	Better go back and fix that!
	###########################
	###########################'''.format(Taskfarm_total, n_per_species*n_species))
	sys.exit()

#################
###########################
##### Initial housekeeping - make the directory your experiment will be in, and subdirectories for each replicate:

setup.directorymaker(logger, targetdir)
for directory in range(1, Taskfarm_total+1):
    setup.directorymaker(logger, '{0}{1:02d}/'.format(targetdir, directory))
setup.directorymaker(logger, '{0}{1}/'.format(targetdir, 'fullpostfiles'))
setup.directorymaker(logger, '{0}{1}/'.format(targetdir, 'truncpostfiles'))

#####First we'll calculate your isotherm and write a pressure.dat file in your target directory
workingdir = 1
for i in range(1, n_species+1):
	for j in range(workingdir, workingdir+n_per_species):
	    satP = setup.pSat(logger, species[i], T)
        istm = setup.isothermcalculator(logger, satP, iso_length, minrelpress, maxrelpress)
        setup.PressureFileWriter(logger, species[i], T, satP, istm, '{0}{1:02d}/'.format(targetdir, j)) #goes in the directory above your individual one, gets symbolic linked later
    workingdir += n_per_species

##### Now to make some control files for you
workingdir = 1
for i in range(1, n_species+1):
	for directory in range(workingdir, workingdir+n_per_species):
        setup.GcmcControlChanger(logger, species[i], sorb_el_list, T, n, framework, '{0}{1:02d}/'.format(targetdir, directory), n_iterations, restart, ctrl_file_name, pressure)
        for j, value in enumerate(istm): #now I make the extra .ctr files to give you final .xyz files for each simulaiton
            setup.GcmcControlChanger(logger, species[i], sorb_el_list, T, '1', framework, '{0}{1:02d}/'.format(targetdir, directory), '1', 'RESTARTFILE {0}.{1}.res.{2}'.format(framework, species[-1], j+1), '{0}kpa_restart.ctr'.format(value), value) #makes your control files for all of your pressure points
        setup.PostControlChanger(logger, species[i], T, n, framework, '{0}{1:02d}/'.format(targetdir, directory), '0')
        setup.PostControlChanger(logger, species[i], T, n, framework, '{0}{1:02d}/'.format(targetdir, directory), '60')
        setup.GcmcRunWriter(logger, species[i], T, framework, parentdir, '{0}{1:02d}/'.format(targetdir, directory), istm, '')#the last variable is the relative location of your interactions files
        setup.AtmAtmMover(logger, sorb_el_list, forcefield, coultype, hicut, frameowrk, '{0}{1:02d}/'.format(targetdir, directory))
        setup.SorbSorbWriter(logger, species[i], sorb_el_list,framework,'{0}{1:02d}/'.format(targetdir, directory), pmap, emap)
        setup.IntraWriter(logger, species[i], sorb_el_list, framework, '{0}{1:02d}/'.format(targetdir, directory))

setup.TaskfarmRunWriter(logger, 'multiple', T, framework, parentdir, targetdir, Taskfarm_total)


logger.info('''
##################################################################
So I've done the following:
Created the directory {0} and {1} subdirectories for taskfarming.
Calculated and written isotherms for all {2} species you've got, {3} times over each.
Created a gcmc.ctr file, 2 post control files, and a run.gcmc file for each subdirectory.
Created a taskfarm document and a backup, and finally run.taskfarmer in {0}. It might complain about not having a python script, but don't worry about that.
Finally I wrote an atom_atom_file, sorb_sorb_file, and intramolecular_file into each subdirectory.
You're good to go!
##################################################################
'''.format(targetdir, Taskfarm_total, n_species, n_per_species))