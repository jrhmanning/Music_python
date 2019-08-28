import numpy as np
import matplotlib.pyplot as plt

isotherms = np.zeros((17, 20)) # this is a matrix of a 16x 20-point isotherms with identical pressure values
def isotherm(page): #this function scrapes subdirectories of name for isotherm data for a given isotherm file name  
    z = 0
    z = int(page)
    print(z)
    drx = str("./%02d/isotherm.Chloroform_v1" % page) # defines the isotherm 
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
for i in range (1,17):
    isotherm(i)
    plt.semilogx(isotherms[0],isotherms[i], 'ro')

#print(isotherms)
plt.title('Chloroform in IRMOF-1 at 298K')
plt.xlabel('Pressure (kPa), Psat = 26 kPa')
plt.ylabel('Quantity adsorbed (molecules/unit cell)')
plt.show()
