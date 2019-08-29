######This python script take in a list of your elements and some forcefield information, then spits out a lot of control files to make maps of them all
#It requires the following variables: 
#species - the name of your sorbent species. Since we're making a map, this is actually one of your sorbent elements
#framework - the name of your MOF (or other porous material
#################
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
species = ['XXXX']
sorb_el_list = ['''### your sorbent element list here###''']
framework = 'YYYY'
targetdirectory="./mapgen/{0}".format(species[-1])
###########
for i in sorb_el_list:
	setup.mapctrlfilewriter(logger, i, framework, targetdirectory)
