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

#### simulation basics
species = ['MeOH'] #Name of your sorbent, as a list
sorb_el_list = ['Me_s', 'O_s', 'H_s'] #different elements in yoru sorbent, as a list
framework = 'IRMOF1step0' #name of your framework
T = 298 #simulation temperature

#### Directory stuff
username = 'jrhm21' #your username, for writing paths and slurm outputs
parentdir = Path('/home/r/{0}/scratch/02_music_squished_IRMOF1/'.format(username)) #The directory this python script is in
targetdir = Path('./mapgen/tension/{0}/{1}'.format(species[-1], framework)) #The directory you want to output to


##### forcefield information
forcefield = 'MeOH' # lookup code for you interaction parameters in Forcefield.py 
hicut = 18 #atom-atom interaciton cutoff, in angstrom
coultype = None #fluid-fluid interactions calcualtion type
pmap = False #do you have fluid-framework pmaps
emap = False #do you have fluid-framework emaps?
FluFra = True #do you have explicit fluid-framework lennard-jones interactions?


#################
###########################
##### Initial housekeeping - make the directory your experiment will be in, and subdirectories for each replicate:

setup.directorymaker(logger, targetdir)
##### Now to make some control files for you
setup.MapMakeRunWriter(logger, species[-1], sorb_el_list, '{0}/'.format(parentdir), targetdir)
for i in sorb_el_list:
    setup.mapctrlfilewriter(logger, i, framework, targetdir)
setup.AtmAtmMover(logger, sorb_el_list, forcefield, coultype, hicut, framework, targetdir, FluFra)
setup.SorbSorbWriter(logger, sorb_el_list, sorb_el_list,framework,targetdir, pmap, emap)
setup.IntraWriter(logger, species[-1], sorb_el_list, framework, targetdir)

logger.info('''
##################################################################
So I've done the following:
First I made directory {0}.
Then I produced a map making control file and runscript for you to make it.
Finally I wrote an atom_atom_file, sorb_sorb_file, and intramolecular_file into {0}.
You're good to go!
##################################################################
'''.format(targetdir))
