
import numpy as np


def isothermextract(page, species): #this function scrapes subdirectories of name for isotherm data for a given isotherm file name
    isotherms = np.zeros((17, 20)) # this is a matrix of a 16x 20-point isotherms with identical pressure values
    z = 0
    z = int(page)
    print(z)
    drx = str("./{0:02d}/isotherm.{1}".format(page,species)) # defines the isotherm
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

