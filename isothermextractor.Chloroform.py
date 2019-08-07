from os import listdir
import numpy as np

isotherms = np.zeros((17, 20)) # this is a matrix of a 16x 20-point isotherms with identical pressure values

def isothermextract(page): #this function scrapes subdirectories of name for isotherm data for a given isotherm file name
    z = 0
    z = int(page) #sets z to be your taskfarm run
#    print(z)
    drx = str("./{0:02d}/isotherm.Chloroform" .format(page)) # defines the isotherm
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
    filenames = listdir("./{0:02d}/" .format(i))
    for file in filenames:
        if file == "isotherm.Chloroform": #change the np.zeroes to be your results
            isothermextract(i)
            SuccessCount += 1
    if SuccessCount != 1: #Checks the abovle was successful
        print("WARNING! Something went wrong in simulation {0:02d}" .format(i)) #complains if not

x = isotherms[0] #This is now getting a mean and stdev of your results. Still works even if some of your simulations don't!
avg=[]
stdev=[]

avgisotherms = np.zeros((3, 20)) #now makes a 3x20 matrix to hold your averaged data and stdevs
for i in range (0, 20):
    spread = []
    for thing in range (2,17):
        if isotherms.item((thing, i)) != 0: #prevents your failed simulations from messing things up
            spread.append(isotherms.item((thing, i)))
    avg.append(np.mean(spread)) #gets your mean
    stdev.append(np.std(spread)) #gets your stdev
avgisotherms[0] = x
avgisotherms[1] = avg
avgisotherms[2] = stdev
np.savetxt("../IRMOF1.358.Chloroform.DMFToChloroform.csv", avgisotherms, delimiter=",") #output it to this file