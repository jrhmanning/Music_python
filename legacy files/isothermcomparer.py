import numpy as np
import matplotlib.pyplot as plt

isotherms = np.zeros((17, 20)) # this is a matrix of a 16x 20-point isotherms with identical pressure values
def isotherm(page): #this function scrapes subdirectories of name for isotherm data for a given isotherm file name  
    z = 0
    z = int(page)
    print(z)
    drx = str("./298k_repeat/%02d/isotherm.Chloroform_v1" % page) # defines the isotherm 
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

plt.semilogx(0.001, 0, 'ro', label='no coulombic')

drx2 = str("../chloroform_v2/298k/isotherm.Chloroform_v2") # defines the isotherm
print(drx2)
f2 = open(drx2, "r")
a = []
b = []
for curline in f2:
    if "#" not in curline:
        fish = curline.split()
        if len(fish) == 2:
            a.append(fish[0])
            b.append(fish[1])
f2.close()
plt.semilogx(a, b, 'bo', label='fluid-fluid coul')


drx3 = str("./ccl3hxptl.txt") # defines the isotherm
print(drx3)
f3 = open(drx3, "r")
c = []
d = []
for curline in f3:
    if "#" not in curline:
        fish = curline.split()
        if len(fish) == 2:
            c.append(fish[0])
            d.append(fish[1])
f3.close()
plt.semilogx(c, d, 'g--', label='xptl')

isotherms2 = np.zeros((17, 20)) # this is a matrix of a 16x 20-point isotherms with identical pressure values
def isotherm2(page4): #this function scrapes subdirectories of name for isotherm data for a given isotherm file name
    z2 = 0
    z2 = int(page4)
    print(z2)
    drx4 = str("./313k/%02d/isotherm.Chloroform_v1" % page4) # defines the isotherm
    print(drx4)
    f4 = open(drx4, "r")
    x4 = []
    page4 = []
    for curline in f4:
        if "#" not in curline:
            fish = curline.split()
            if len(fish) == 2:
                x4.append(fish[0])
                page4.append(fish[1])
    f4.close()
    isotherms2[0] = x4
    isotherms2[z2] = page4
for i in range (1,17):
    isotherm2(i)
    plt.semilogx(isotherms2[0],isotherms2[i], 'rx')

plt.semilogx(0.001, 0, 'rx', label='313k')



#print(isotherms)
plt.title('Chloroform in IRMOF-1 at 298K')
plt.xlabel('Pressure (kPa), Psat = 26 kPa')
plt.ylabel('Quantity adsorbed (molecules/unit cell)')
plt.legend(loc='lower right', shadow=True, ncol=1 )
plt.show()
