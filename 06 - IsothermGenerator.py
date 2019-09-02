######This python script take in a your species name and a dictionary of information about its boiling point, and produces a pressure file for your simulations
#It requires the following variables: 
#species - the name of your sorbent species. Since we're making a map, this is actually one of your sorbent elements
#T - your simulation temperature
#n - the number of isotherm points you're using
#minrelpress - the minimum pressure point you want to simulate, relative to the saturation pressure These are good for zeroing in on your pressure region
#maxrelpress - the maximum pressure point you want to simulate, relative to the saturation pressure
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
species = ['DMF'] #NAme of your sorbent, is used to look up anointe equation parameters in Antoine.py
T = 298 #Temperature for antoine equation
iso_length = 9 #numebr of isotherm points
minrelpress = -1.5 #Minimum pressure to be considered, as a fration of saturation pressure (pMin = pSat*10^minrelpress) 
maxrelpress = -1 #maximum pressure to be considered, accoring to the above equation
targetdir = Path('./experiments/{0}/{1}/'.format(species[-1], T)) #directory you're putting your experiment files into
#################
###########################
satP = setup.pSat(logger, species[-1], T)
istm = setup.isothermcalculator(logger, satP, iso_length, minrelpress, maxrelpress)
#print(istm)
setup.PressureFileWriter(logger, species[-1], T, satP, istm, targetdir)
