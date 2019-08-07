#Hi there! I'm a python script who's job it is to automatically generate isotherms for your Music simulations
#I eat Antoine equation parameters (hardcoded in), sorbent chemical names, temperatures in Kelvin, and desired number of isotherm points
#Using these, I determine the saturation pressure of the chemical Species at your temperature, create an isotherm with 5 set low-pressure points and (n-5) points equidistant between this and your saturation pressure
#Finally, I create a data file called pressure.*Species*.*temperature*.dat with all of this information in
#Written by Joe Manning with the help of Gaee Donval and Mat Tolladay, March 2019 in the university of Bath, UK

Antoine = {                                                 #This is your library of Antoine equatino parameters. The key is the name of the Species that will be written
    "CCL4_v1":(4.02291, 1221.781, -45.739, 193, 350),          #Indices 0-2 are antoine parameters A, B, and C taken from NIST,
    "Chloroform_v1":(4.20772, 1233.129, -40.953, 215, 334), #Indices 3 and 4 are the Antoine equation validity range, as quoted on the NIST website
    "Chloroform_v2":(4.20772, 1233.129, -40.953, 215, 334),
    "DCM_v1":(4.53691, 1327.016, -20.474, 233, 313),
    "Chloromethane_v1":(4.91858, 1427.529, 45.137, 303, 416),
    "Chloromethane2_v1":(4.22507, 951.561, -23.468, 198, 278), 
    "Methane_v1":(4.22061, 516.689, 11.223, 110, 190),
    "Methanol_v1":(5.20409, 1581.341, -33.5, 288, 356.8),
    "THF_v1":(4.12118, 1202.942, -46.818, 296.29, 372.8)
}

FixedPoints = [0.01, 0.03, 0.05, 0.07, 0.09]              #These are the fixed low pressure points I'll include atthe start fo your isotherm


def pSat(Species, T):                                   #This function checks the Species is there and that you're in the right temperature range. 
    Press = 0
    if Species in Antoine:                              #It then reads out your saturation pressure in kPa for the user benefit and to get the rest of the program working
        if Antoine[Species][3] <= T <= Antoine[Species][4]:
            print("I have Antoine parameters for %s at %s K." % (Species, T))
            Press = 100*(10**(Antoine[Species][0]-(Antoine[Species][1]/(Antoine[Species][2]+T))))
            print("At this temperature, the saturation pressure is {0:1.3f} kPa." .format(Press))
            rnge = 1
        elif T < Antoine[Species][3]:
            print("WARNING, that temperature is too low for my Antoine parameters for %s (my limit is %s K). I'll persevere anyway." % (Species, Antoine[Species][3]))
            Press = 100*(10**(Antoine[Species][0]-(Antoine[Species][1]/(Antoine[Species][2]+T))))
            print("At this temperature, the saturation pressure is possibly {0:1.3f} kPa. But who knows, you range-busting maniac." .format(Press))
            rnge = 0
        elif Antoine[Species][4] < T:
            print("WARNING, that temperature is too high for my Antoine parameters for %s (my limit is %s K). I'll persevere anyway" % (Species, Antoine[Species][4]))
            Press = 100*(10**(Antoine[Species][0]-(Antoine[Species][1]/(Antoine[Species][2]+T))))
            print("At this temperature, the saturation pressure is possibly {0:1.3f} kPa. But who knows, you range-busting maniac." .format(Press))
            rnge = 2
        else:
            print("Something weird happened, I've probably got a bug. Oops!")
            exit()
    else:
        print("I'm sorry, I don't have that Species in my database.")
        exit()
    return Press

def isothermcalculator(Pressure, n):                                  #This function produces the isotherm as a data list, using pSat calculated by function pSat and the user defined isotherm length n
    m = n-len(FixedPoints)-1                                         #Pressure points re linearly distributed above 0.09 kPa to pSat, which may eb a bad idea. Who knows?
    AllPoints = []
#    print(FixedPoints[-1])
#    print(Pressure-FixedPoints[-1])
    if Pressure-FixedPoints[-1] < FixedPoints[-1]:                  #This subroutine is used when the pSat isn't much higher than the fixed points. At this point it ignores them and generates more points to compensate
        for point in range(1,m+2+len(FixedPoints)):
            press = point*(Pressure)/(m+len(FixedPoints))
            press = round(press, 3)
            AllPoints.append(press)
            AllPoints.sort()
    else:
        for thing in FixedPoints:                                 #This is the general loop for printing pressures the rest of the time
            AllPoints.append(thing)
        for point in range(1,m+2):
            press = point*(Pressure-FixedPoints[-1])/m
            press = round(press, 3)
            AllPoints.append(press)
            AllPoints.sort()
    isotherm = AllPoints
#    print(isotherm[2])
#    print(len(isotherm))
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
##Species = str("CCL4")
##T = 298K
##n = 20
###############################################################################################################################################################
##This bit does the actual work          
#Pressure = pSat(Species, T)
##print("{0:1.3f}" .format(pSat("CCL4", 298)))
#isotherm = isothermcalculator(Pressure, n)
##print(isotherm)
#dirout = "./"
#directory = 0
#PressureFileWriter(Species, T, Pressure, isotherm)
