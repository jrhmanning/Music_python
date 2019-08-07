######This python script takes in your species, elements and framework, then spits out an intramolecular file for music
#It requires the following variables: 
#species - the name of your sorbent species. Here it needs to be a list, otherwise things go badly.
#elements - a list of the elemetns involved in your sorbent species, for maps
#framework - the name of your MOF (or other porous material
#################
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
framework = 'YYYY'
sorb_el_list = ['''### your species list here###''']
#################
#Writes your intramolecular file
def IntraWriter(species, elements, framework, dxout = "./mapgen/"):
    with open("{0}intramolecular_file" .format(dxout), "w") as file:
        for i in species:
            file.write("Intra: {0}\n" .format(i)) #no intramolecular interactions at all
        for i in list(elements): #lets you create a map of every element in your system. Having the extra lines in your production script don't matter
            file.write("Intra: {0}\n".format(i))
        file.write("Intra: {0}" .format(framework)) #No framework intramolecular interactions happening
    logger.debug("Intra file written!")
##################
ntraWriter(species, sorb_el_list, framework, './experiments/{0}/298/'.format(species))