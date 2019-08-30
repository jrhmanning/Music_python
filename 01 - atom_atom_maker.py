######This python script takes in a list of your elements and some forcefield information, then spits out an atom_atom_file for music
#It requires the following variables: 
#elements - a list of the elemetns involved in your sorbent species
#MOF_el_list - a list of your framework atom types
#IntParams - a dictionary of your IntParams
#forcefield - the name of the forcefield, which is a suffix to your species atoms in IntParams
####
# If you're wanting to generate your element list from [species].mol, you'll be able to use this too.

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
handler = logging.FileHandler('{0}/atom_atom_maker.log'.format('.'))
handler.setLevel(logging.DEBUG)
handler2 = logging.StreamHandler()
handler2.setLevel(logging.INFO)  ##If the messages are getting annoying, change this to logging.DEBUG
logger.addHandler(handler)
logger.addHandler(handler2)
logger.debug('-----------------------------------')
logger.debug('Generating new files on {0}-{1}-{2} at {3}:{4}'.format(now.year, now.month, now.day, now.hour, now.minute))
logger.debug('-----------------------------------')
##################
####This section lists your user-defined input parameters:

species = 'DMFAA'
forcefield = 'DMFAA' #for example
framework = 'IRMOF1'
coultype = 'Ewald' # 'Ewald', 'Wolf', or 'None' for different fluid-fluid coulombic interactions)
cutoff=18
speciespath = Path('./{0}.mol'.format(species))
targetdirectory = './'
#Or if you're going to define tyour atom list manually, this list will replace the below section
el_list=None


logger.info('''These are my input parameters:
I'm making an Atom-Atom interaction file in directory {0} for a simulation of sorbent {1} on sorbate {2}!'''.format(targetdirectory, species, framework))
if el_list:
    logger.info('You\'ve told me the names of your sorbent atoms already, they are: {0}'.format(el_list))
else:
    logger.info('You\'ve not told me the names of your sorbent atoms already, so I\'ll work it out now.')
    logger.info('To do this, I\'ll look for {0}.mol in directory {1} and read it for atom names'.format(species, speciespath))
logger.info('''Then I\'ll write a atom-atom interaction file using LB mixing rules, \
Ewald summation between sorbent molecules, a cutoff of {1} Angstrom, and finally\
print it out to the directory {0}'''.format(targetdirectory, cutoff))
logger.debug('NB the above info is not auto-updated if you change the Ewald summation and LB mixing, so it might be incorrect.')

################
####This section creates a list of your sorbent atom types from a .mol file, 
if el_list:
    logger.info('I received a user-defined sorbent list, so I\'m not going to go looking for a new one.')
else:
    el_list = set()
    if speciespath.exists():
        with speciespath.open() as file:
            for line in file:
                if len(line.split()) == 8: #should be 8 (index x y z name charge ? ?), but you might have added comments
                    logger.info(line.split()[4])
                    el_list.add(line.split()[4])
        logger.info("Atom list found is: " + str(el_list))
    else:
        logger.warning('I didn\'t find a file named {0}.mol, so I don\'t have any sorbent atoms to write!'.format(species))
sorb_el_list = [] #i want to differentiate between framework elements and sorbent elements here, so I rename them as X_s
for i in el_list:
    sorb_el_list.append("{0}_{1}" .format(i.split("_")[0], 's')) #Sometimes I'll call them X_something esle to denote their forcefield for instance. This hould cut it back to just X_s
logger.info('The list of sorbent atoms I\'ve found are: {0}'.format(sorb_el_list))
################

###############


setup.AtmAtmMover(
    logger, #this variablwe logs out to the one made by this file
    sorb_el_list, #Your elements to be considered in the fluid
    forcefield, #Your forcefield tag for these elements
    coultype, #Your fluid-fluid coulombic interaction model
    cutoff, #Your nonbonded interaction cutoff distance
    framework, #Your framework name, to look up in musicpy.Forcefield.MOF_el_list
    targetdir #Where the atom_atom_interaction file will be spit out to
)
