#Hi there! I'm a python script who's job it is to automatically generate isotherms for your Music simulations
#I eat Antoine equation parameters (hardcoded in), sorbent chemical names, temperatures in Kelvin, and desired number of isotherm points
#Using these, I determine the saturation pressure of the chemical Species at your temperature, create an isotherm with 5 set low-pressure points and (n-5) points equidistant between this and your saturation pressure
#Finally, I create a data file called pressure.*Species*.*temperature*.dat with all of this information in
#Written by Joe Manning with the help of Gaee Donval and Mat Tolladay, March 2019 in the university of Bath, UK

from math import log

Antoine = {                                                 #This is your library of Antoine equatino parameters. The key is the name of the Species that will be written
    "CCL4":(4.02291, 1221.781, -45.739, 193, 350),          #Indices 0-2 are antoine parameters A, B, and C taken from NIST,
    "Chloroform":(4.20772, 1233.129, -40.953, 215, 334), #Indices 3 and 4 are the Antoine equation validity range, as quoted on the NIST website
    "DCM":(4.53691,1327.016,-20.474, 233, 313),
    "Chloromethane":(4.91858, 1427.529, 45.137, 303, 416),
    "Chloromethane2":(4.22507, 951.561, -23.468, 198, 278), 
    "Methane":(4.22061, 516.689, 11.223, 110, 190),
    "Methanol":(5.20409, 1581.341, -33.5, 288, 356.8),
    "THF":(4.12118, 1202.942, -46.818, 296.29, 372.8)
}

FixedPoints = [0.01, 0.03, 0.05, 0.07, 0.09]              #These are the fixed low pressure points I'll include atthe start fo your isotherm


def pSat(Species, T):#This function checks the Species is there and that you're in the right temperature range. 
    Chemical = Species.split("_")[0]
    print(Chemical)
    Press = 1 #pressure defaults to 1 kPa
    if Chemical in Antoine:                              #It then reads out your saturation pressure in kPa for the user benefit and to get the rest of the program working
        if Antoine[Chemical][3] <= T <= Antoine[Chemical][4]:
            print("I have Antoine parameters for {0} at {1} K." .format(Chemical, T))
        elif T < Antoine[Chemical][3]:
            print("WARNING, that temperature is too low for my Antoine parameters for {0} (my limit is {1} K). I'll persevere anyway." .format(Chemical, Antoine[Chemical][3]))
        elif Antoine[Chemical][4] < T:
            print("WARNING, that temperature is too high for my Antoine parameters for {0} (my limit is {1} K). I'll persevere anyway." .format(Chemical, Antoine[Chemical][4]))
        else:
            print("Something weird happened, I've probably got a bug. Oops!")
        Press = 100*(10**(Antoine[Chemical][0]-(Antoine[Chemical][1]/(Antoine[Chemical][2]+T))))
    else:
        print("I'm sorry, I don't have that Species in my database.")
    return Press

def isothermcalculator(Pressure, n):                                  #This function produces the isotherm as a data list, using pSat calculated by function pSat and the user defined isotherm length n
    m = n-len(FixedPoints)-1                                         #Pressure points re linearly distributed above 0.09 kPa to pSat, which may eb a bad idea. Who knows?
    AllPoints = []
#    print(FixedPoints[-1])
#    print(Pressure-FixedPoints[-1])
    if FixedPoints[-1]/Pressure < 0.001:                  #This subroutine is used when the pSat isn't much higher than the fixed points. At this point it sets the minimum relative pressure to 10**-5 and makes a log linear range between this and 0 (+ a bit)
        print('Defaulting to range from 10e-5 to 10e0')
	minrelpress = -5
        for point in range(1,n+1):
            relpress = minrelpress*(1-(point/(n+1)))
            press = Pressure*10**relpress
            press = round(press, 3)
            AllPoints.append(press)
            AllPoints.sort()
    else:
	print('Working between my fixed points up to {0} kPa and {1}' .format(FixedPoints[-1], round(Pressure)))
        minrelpress = log((FixedPoints[-1])/Pressure, 10)
	#print(minrelpress)
        for thing in FixedPoints:                                 #This is the general loop for printing pressures the rest of the time
            AllPoints.append(thing)
	#print(AllPoints)
        for point in range(1,m+2):
	    #print(point)
            #print(float(1-(float(point)/(m+1))))
	    relpress = minrelpress*(1-(float(point)/(m+1)))
            #print(relpress)
	    press = Pressure*10**relpress
           # print(press)
	    AllPoints.append(press)
        AllPoints.sort()
    isotherm = AllPoints
#    print(isotherm[2])
#    print(isotherm)
    return isotherm
    

def PressureFileWriter(Species, T, Pressure, isotherm, directory, dirout):                 # This function writes your file to a specific directory, but currently doesn't say if your temp is out of range. Annoying!!#
#    print(isotherm)
    f = open("%s/%02d/pressure.%s.%s.dat" % (dirout, directory, Species, T), "w")
    f.write("{0!s} #PSat at {1!s} K is {2:1.3f} kPa. \n" .format(Species, T, Pressure))
    f.write(str(len(isotherm)) + "\n")
    f.write(', '.join(str(thing) for thing in isotherm))
    f.close
################################################################################################################################################################
##User input section
#Species = input("What is your sorbent Species called?")
#type(Species)
#T = eval(input("What temperature (in K) would you like to test?"))
#type(T)
#n = eval(input("How many isotherm points do you want?"))
#type(n)
##Or you could hard code any part of it:
#Species = str("Chloroform_v1")
#T = 298
#n = 20
###############################################################################################################################################################
##This bit does the actual work          
#Pressure = pSat(Species, T)
#print("{0:1.3f}" .format(Pressure))
#isotherm = isothermcalculator(Pressure, n)
#print(isotherm)
#print(len(isotherm))
#dirout = "./"
#directory = 0
#PressureFileWriter(Species, T, Pressure, isotherm)
