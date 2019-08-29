
######This python script takes in pdb files and spits out a .mol file for music
#It requires the following variables: 
# species - the molecule you're considering (also used to read [species].pdb and create [species].mol)
# sourcedir - the directory you're reading [species].pdb from
# (optional) - forcefield to load coulobic information from a forcefield dictionary directly into your .mol file
######
#While working, it creates the following variables:
#atoms - a list of atoms in your species
#bonds - a list of bonds in your species
#bondlengths - a dictionary of the bonds in your molecule, with their respective lengths

import os
import re
import random
import numpy as np
import sys
from math import sqrt, log
from pathlib import Path
import logging
import datetime
import musicpy.setup as setup

now = datetime.datetime.now()


####This section sorts out your messages from this script to the console and a log file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler('{0}/molfilemaker.log'.format('.'))
handler.setLevel(logging.DEBUG)
handler2 = logging.StreamHandler()
handler2.setLevel(logging.INFO)
logger.addHandler(handler)
logger.addHandler(handler2)
logger.debug('-----------------------------------')
logger.debug('Generating new files on {0}-{1}-{2} at {3}:{4}'.format(now.year, now.month, now.day, now.hour, now.minute))
logger.debug('-----------------------------------')
######

species = 'EG_TRAPPE'


atoms, connections = setup.dataextract(logger, species)
bondlengths = setup.connectiontypeswrite(logger, atoms, connections)
setup.molwrite(logger, species, atoms, bondlengths, connections)
