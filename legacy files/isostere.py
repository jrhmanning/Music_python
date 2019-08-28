import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
#import maplotlib as plt
filen = 'OccupancyTotal'
T = [298, 313, 328, 343, 358, 373, 388]
R = 8.314 #J/mol.K
pfile = 'pressure.*.dat'
Occupancy = []

def getOccupancy(temp, filename):
    with open('./{0}/{1}'.format(temp, filename), 'r') as f:
        f.readline()
        linedata = []
        for line in f:
            rawdata = line.strip()
            #print(rawdata)
            linedata = rawdata.split(', ')
            #print(linedata)
            for i in linedata[1:]:
                Occupancy.append(i)
for i in T:
    getOccupancy(i, filen)
jeff = np.array(Occupancy).astype(np.float)
#print(Occupancy)
hist, bin_edges = np.histogram(jeff, 20)

print(hist)
#print(sum(hist))
#print(len(Occupancy))
print('---------------')
print(sorted(bin_edges))

histo = {}
for i in bin_edges:
    histo[i] = []
#print(histo)

Q_st = []
n = []
quant = []

def fillbins(temp, filename1, filename2):
    fugacity = []
    n = []
    for file_name in os.listdir('./{0}/'.format(temp)):
        if fnmatch.fnmatch(file_name, filename2):
            print(file_name)
            break
    #print(file_name)
    with open('./{0}/{1}'.format(temp, file_name), 'r') as f:
        #print(f.readline())
        #f.readline()
        for line in f:
            rawdata = line.strip()
            #print(rawdata)
            fugacity = rawdata.split(', ')
    #print(fugacity)
    #######################
    with open('./{0}/{1}'.format(temp, filename1), 'r') as f:
        f.readline()
        linedata = []
        for count, line in enumerate(f):
            rawdata = line.strip()
            #print(rawdata)
            linedata = rawdata.split(', ')
            #print(count)
            #print(line)
            if count<len(fugacity):
                fug = fugacity[count]
            for i in linedata[1:]:
                for key in histo:
                    if float(i)<float(key):
                        histo[key].append([fug,temp])
                quant.append(i)
                Q_st.append(round(R*temp**2*(np.log(float(fug)*1000)/temp), 3))
for i in T:
    fillbins(i, filen, pfile)
#print(histo)
print(len(quant))
print(len(Q_st))
with open('Q_st_test.csv', 'w') as f:
    for i in range(len(quant)):
        f.write('{0}, {1}\n'.format(quant[i], Q_st[i]))
plt.plot(quant, Q_st)
print('gate 1')
plt.savefig('test.jpg')
print('gate 2')
#plt.show()
