######This python script takes in your species, frameowkr, and the elements you have, then spits out a sorb_sorb_file for music
#It requires the following variables: 
#species - the name of your sorbent species. Here it needs to be a list, otherwise things go badly.
#elements - a list of the elemetns involved in your sorbent species, for maps
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
framework = 'YYYY'
sorb_el_list = ['''### your species list here###''']
targetdir = Path('./experiments/{0}/298/'.format(species))
#########
setup.SorbSorbwriter(logger, species, elements, framework, targetdir, True, True)
