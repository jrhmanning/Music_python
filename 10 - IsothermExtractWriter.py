######This python script take in a list of your species and temperature, then spits out a python script which will analyse the data files produced by music post to get an averaged isotherm across 16 replicates
#It requires the following variables: 
#species - the name of your sorbent species
#elements - a list of the elemetns involved in your sorbent species, for maps
#framework - the name of your MOF (or other porous material
#parentdir = a path to navigate you to your testing directory
#################
import os
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
T = 298
n  20
framework = 'YYYY'
##################
#Writes a python script to extract data from your experiments after taskfarming
def IsothermExtractMover(Species, T, framework, dirout, n = 20):
    with open("{0}isothermextractor.py" .format(dirout), 'w') as file:
        file.write("""from os import listdir
import numpy as np

isotherms = np.zeros((17, {3})) # this is a matrix of a 16x {3}-point isotherms with identical pressure values

def isothermextract(page): #this function scrapes subdirectories of name for isotherm data for a given isotherm file name
    z = 0
    z = int(page) #sets z to be your taskfarm run
#    print(z)
    drx = str("./{{0:02d}}/isotherm.{0}" .format(page)) # defines the isotherm
#    print(drx)
    f = open(drx, "r")
    x = []
    page = [] #redefines your input variable, but somehow still works
    for curline in f:
        if "#" not in curline:
            fish = curline.split()
            if len(fish) == 2:
                x.append(fish[0])
                page.append(fish[1])
    f.close()
    isotherms[0] = x #overwrites the first column of your np.zeros matrix to be your pressure oiubts
    isotherms[z] = page #overwrites column z in your np.zeros matrix with resutls from isotherm z
    return isotherms
#####################################################################################################################################################	
for i in range (1,17): #loops over your 16 isotherms
    SuccessCount = 0 #this variable preents issues downstream if some of your simulations fail
    filenames = listdir("./{{0:02d}}/" .format(i))
    for file in filenames:
        if file == "isotherm.{0}": #change the np.zeroes to be your results
            isothermextract(i)
            SuccessCount += 1
    if SuccessCount != 1: #Checks the abovle was successful
        print("WARNING! Something went wrong in simulation {{0:02d}}" .format(i)) #complains if not

x = isotherms[0] #This is now getting a mean and stdev of your results. Still works even if some of your simulations don't!
avg=[]
stdev=[]

avgisotherms = np.zeros((3, {3})) #now makes a 3x{3} matrix to hold your averaged data and stdevs
for i in range (0, {3}):
    spread = []
    for thing in range (2,17):
        if isotherms.item((thing, i)) != 0: #prevents your failed simulations from messing things up
            spread.append(isotherms.item((thing, i)))
    avg.append(np.mean(spread)) #gets your mean
    stdev.append(np.std(spread)) #gets your stdev
avgisotherms[0] = x
avgisotherms[1] = avg
avgisotherms[2] = stdev
np.savetxt("../{2}.{1}.{0}.results.csv", avgisotherms, delimiter=",") #output it to this file""" .format(Species, T, framework, n))
    logger.debug("Python isotherm extractor written!")
##############
IsothermExtractMover(species[-1], T, framework, './experiments/{0}/{1}/'.format(species[-1], T), n)