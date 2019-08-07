######This python script takes in your species, frameowkr, and the elements you have, then spits out a sorb_sorb_file for music
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

##################
#This writes your sorb-sorb file telling music which intermolecular interactions are on
def SorbSorbWriter(species, elements,framework, dxout = "./mapgen/", pmap = True, emap = False):
    filename = "sorb_sorb_file"
    with open("{0}{1}" .format(dxout, filename), 'w') as file:
        file.write("{0} {0} NCOUL OFF\n{0} {0} COUL OFF\n\n" .format(framework)) #MOF-MOF interactions are off
        if pmap == True: #using a MOF-fluid LJ potential map if true
            for i in species:            
                file.write("{0} {1} NCOUL MAP@{1} FAST ".format(i, framework)) #declares it'll use a map
                for j in elements:
                    file.write("{0}@PMAP@{1}.{0}.p " .format(j, framework)) #writes your map names. WARNING - it's a 18 character limit on map names, and no '-' allowed!
                file.write("\n")
        else:
            for i in species:
                file.write("{0} {1} NCOUL BASIC LJ FAST\n".format(i, framework)) #pairwise fluid-framework interactions
        if emap == True: #using a MOF-fluid coulomb map if true
            for i in species:
                file.write("{0} {1} COUL MAP@{1} FAST ".format(i, framework)) #declares it'll use a map
                for j in elements:
                    file.write("{0}@EMAP@{1}.{2}.e " .format(j, framework, 'Probe')) #writes your map names. WARNING - it's a 18 character limit on map names, and no '-' allowed!
                file.write("\n")
        else:
            for i in species:
                file.write("\n{0} {1} COUL OFF\n" .format(i, framework)) #pairwise coulomb interactionsdon't work. Trust me, this is better.
            logger.warning("Framework-fluid coulombic interactions are off") #warns you your coulomb fluid-framework stuff isn't happeneing
        for i in species:
            file.write("\n{0} {0} NCOUL BASIC LJ FAST\n{0} {0} COUL BASIC WFCOUL FAST\n\n" .format(i))  #fluid fluid interactions, with Wolf coulombic interactions 
    logger.debug("Sorb-Sorb file written!")    
#########
SorbSorbwriter(species, elements, framework, './experiments/{0}/298/'.format(species), True, True)