######This python script take in a list of your elements and some forcefield information, then spits out a control file to manage your gcmc simulations
#It requires the following variables: 
#species - the name of your sorbent species
#elements - to tell your control file how many and what elements to expect
#n - your number of smiulation points
#T - your simulaiton temperature
#framework - the name of your MOF (or other porous material
#iterations = the number of iterations you're using
#Restart = if you want to restart your simulation from a certian point, or to get a snapshot
#name - useful for autogenerating control files for snapshots
#Pressure - useful if you want to make a snapshot
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
species = ['DMFYang']
sorb_el_list = ['AldH_s', 'O_s', 'Carb_s', 'C_s', 'N_s', 'H_s']
framework = 'IRMOF1'
T = 298
n = 9
targetdir = Path('{0}/{1:02d}/'.format(xptpath, directory))
##################
setup.GcmcControlChanger(logger, species[-1], sorb_el_list, T, n, framework, tagretdir, '1000000') #, 'RESTARTFILE {0}.{1}.res.{2}'.format(framework, species[-1], 20)) #for production runs, these 3 optional extra arguments are for running restarts.
#for i, value in enumerate(istm): #now i make the extra .ctr files for yoru subdirectories
#    setup.GcmcControlChanger(logger, species[-1], sorb_el_list, T, 1, framework, targetdir, '1', 'RESTARTFILE {0}.{1}.res.{2}'.format(framework, species[-1], i+1), '{0}kpa_restart.ctr'.format(value), value) #makes your control files for all of your pressure points

