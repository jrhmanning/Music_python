######This python script take in a your species name and a dictionary of information about its boiling point, and produces a pressure file for your simulations
#It requires the following variables: 
#species - the name of your sorbent species. Since we're making a map, this is actually one of your sorbent elements
#T - your simulation temperature
#n - the number of isotherm points you're using
#minrelpress - the minimum pressure point you want to simulate, relative to the saturation pressure These are good for zeroing in on your pressure region
#maxrelpress - the maximum pressure point you want to simulate, relative to the saturation pressure
#################
from math import sqrt, log
import logging
from pathlib import Path
import datetime
now = datetime.datetime.now()
#import numpy as np
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
species = ['DMF']
FixedPoints= []
T = 298
n = 9
minrelpress = -1.5
maxrelpress = -1
Antoine = {                                                 #This is your library of Antoine equation parameters. The key is the name of the Species that will be written
    "CCl4":(4.02291, 1221.781, -45.739, 193, 350),          #Indices 0-2 are antoine parameters A, B, and C taken from NIST,
    "Chloroform":(4.20772, 1233.129, -40.953, 215, 334), #Indices 3 and 4 are the Antoine equation validity range, as quoted on the NIST website
    "DCM":(4.53691,1327.016,-20.474, 233, 313),
    "Chloromethane":(4.91858, 1427.529, 45.137, 303, 416),
    "Chloromethane2":(4.22507, 951.561, -23.468, 198, 278), 
    "Methane":(4.22061, 516.689, 11.223, 110, 190),
    "Methanol":(5.20409, 1581.341, -33.5, 288, 356.8),
    "THF":(4.12118, 1202.942, -46.818, 296.29, 372.8),
    "DMF":(3.93068, 1337.716, -82.648, 303, 363),
    "Nitrogen":(3.7362, 264.651, -6.788, 63.14, 126)
    }

#################
#calculates your saturation pressure using the antoine equation parameters from the preamblem, so you can autogenerate pressures
def pSat(Species, T):#This function checks the Species is there and that you're in the right temperature range. 
    Chemical = Species.split("_")[0]
    logger.debug(Chemical)
    Press = 1 #pressure defaults to 1 kPa
    if Chemical in Antoine:                              #It then reads out your saturation pressure in kPa for the user benefit and to get the rest of the program working
        if Antoine[Chemical][3] <= T <= Antoine[Chemical][4]:
            logger.info("I have Antoine parameters for {0} at {1} K." .format(Chemical, T))
        elif T < Antoine[Chemical][3]:
            logger.warning("WARNING, that temperature is too low for my Antoine parameters for {0} (my limit is {1} K). I'll persevere anyway." .format(Chemical, Antoine[Chemical][3]))
        elif Antoine[Chemical][4] < T:
            logger.warning("WARNING, that temperature is too high for my Antoine parameters for {0} (my limit is {1} K). I'll persevere anyway." .format(Chemical, Antoine[Chemical][4]))
        else:
            logger.warning("Something weird happened, I've probably got a bug. Oops!")
        Press = 100*(10**(Antoine[Chemical][0]-(Antoine[Chemical][1]/(Antoine[Chemical][2]+T))))
    else:
        logger.warning("I'm sorry, I don't have that Species in my database.")
    return Press

#Creates a log-linear pressure isotherm using the saturation pressure from pSat and preamble defined number of points/predefined points 
def isothermcalculator(Pressure, n, minrelpress, maxrelpress):                                  #This function produces the isotherm as a data list, using pSat calculated by function pSat and the user defined isotherm length n
    m = n-len(FixedPoints)-1                                         #Pressure points re linearly distributed above 0.09 kPa to pSat, which may eb a bad idea. Who knows?
    AllPoints = [] #your list of total pressures
    for point in range(1,n+1):
        relpress = minrelpress+point*abs((minrelpress-maxrelpress)/(n)) #Calculates your partial pressure list between minrelpress and 0. Dividing the (point/float(n+1)) statement by a number lowers your max pressure to a fraction of the sat pressure
        print(relpress)
        press = Pressure*10**relpress #converts the above into absolute pressures
        print(press)
        press = round(press, 5) #Rounds the float to 3 decimal points, for simplicity
        AllPoints.append(press)
        AllPoints.sort() #sorts your points into ascending order
    isotherm = AllPoints
    logger.info("PSat = {0}KPa, Isotherm = {1}".format(Pressure, isotherm))
    return isotherm
    
#this function write s a .dat file based on the pressures you've calculated in isothermcalculator
def PressureFileWriter(Species, T, Pressure, isotherm, dirout):                 # This function writes your file to a specific directory, but currently doesn't say if your temp is out of range. Annoying!!#
#    print(isotherm)
    with open("{0}pressure.{1}.{2}.dat" .format(dirout, Species, T), "w") as file:
        file.write("{0!s} #PSat at {1!s} K is {2:1.3f} kPa. \n" .format(Species, T, Pressure))
        file.write(str(len(isotherm)) + "\n")
        file.write(', '.join(str(thing) for thing in isotherm))
###########################
satP = pSat(species[-1], T)
istm = isothermcalculator(satP, n, minrelpress, maxrelpress)
print(istm)
PressureFileWriter(species[-1], T, satP, istm, './'.format(species[-1], T))
