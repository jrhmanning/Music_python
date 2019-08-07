############################################################################################################################################################################
####This python file takes a post-output and makes some xyz files for meta-analysis
####The 3 variables are molecule name ("your_species_here"), Temperature ("your_temperature_here"), and directory ("your_tree_here")
#######################################################################################################################################################

import numpy as np
isotherms = np.zeros((17, 20)) # this is a matrix of a 16x 20-point isotherms with identical pressure values
def isothermextract(page): #this function scrapes subdirectories of name for isotherm data for a given isotherm file name
    z = 0
    z = int(page)
    print(z)
    drx = str("./%02d/isotherm.CCL4_v1" % page) # defines the isotherm
    print(drx)
    f = open(drx, "r")
    x = []
    page = []
    for curline in f:
        if "#" not in curline:
            fish = curline.split()
            if len(fish) == 2:
                x.append(fish[0])
                page.append(fish[1])
    f.close()
    isotherms[0] = x
    isotherms[z] = page
    return isotherms
#####################################################################################################################################################	
for i in range (1,17):
    isothermextract(i)
x = isotherms[0]
avg=[]
stdev=[]

avgisotherms = np.zeros((3, 20))
for i in range (0, 20):
	print(i)
	spread = []
	
	for thing in range (2,17):

		spread.append(isotherms.item((thing, i)))
	avg.append(np.mean(spread))
	stdev.append(np.std(spread))
avgisotherms[0] = x
avgisotherms[1] = avg
avgisotherms[2] = stdev
np.savetxt("../298.CCL4_v1.results.csv", avgisotherms, delimiter=",")
