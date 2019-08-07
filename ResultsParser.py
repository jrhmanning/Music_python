#organise by MC step (denoted by ./IRMOF1.DCM.Chloroform.con.x)
#Then by type (total energy, interaction strength, loading* no.species)
# each of those 3 is separated by empty lines, and loading* no species by commented lines
import os
import numpy as np
import math
#from scipy import stats
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from pymbar import timeseries #Written by Michael Shrits and Paul chodera - cite http://dx.doi.org/10.1063/1.2978177
species = 'DMF_v4'
Template = 'IRMOF1.{0}.con'.format(species)
directory = './postfiles'
#print(alldata)
#print('\n'.join(sorted(confiles)))
def GetEnergy(line):
    data = {}
    linedata = []
    while True:
        #print(alldata[line])
        linedata = alldata[line].split()
        #print(linedata)
        data[linedata[0]] = linedata[2]
        line +=1
        if len(alldata[line]) == 0:
            break
    
    datalist = np.zeros(len(data))
    for i, key in enumerate(data):
        datalist[i] = data[key]
#    [t0, g, Neff_max] = timeseries.detectEquilibration(datalist)
#    ################print('Energy atuocorrelation analysis:\nt0 = {0}, g = {1}, N_max = {2}'.format(t0, g, Neff_max))
#    aveData = np.mean(datalist[t0:])
#    stdData = np.std(datalist[t0:])
    finval = float(datalist[-1])
    accept = []
    N_accept = 1
    accept.append(finval)
    for i in datalist[-2::-1]:
        #print((float(i)-finval)/finval)
        #print(0.05*finval)
        if abs((float(i)-finval)) < abs(0.01*finval):
            accept.append(i)
            N_accept +=1
        else:
            break
    ###########print(N_accept)
    CrudeAvg=np.mean(accept)
    CrudeStd = np.std(accept)
    ############print('autocorrelation -> {0:.0f} pm {1:.0f}, crude -> {2:.0f} pm {3:.0f}'.format(aveData, stdData, CrudeAvg, CrudeStd))
    #points = []
    #for i in datalist[t0::math.ceil(g)]:
    #    points.append(i)
    #print('points = {0}'.format(len(points)))
    #print('The average energy is {0} kJ, with a standard deviation of {1} kJ'.format(aveData, stdData))
    return N_accept, 100, 100, CrudeAvg, CrudeStd, line
def GetInteraction(line):
    Coulombic= {}
    NonCoul = {}
    linedata = []
    while True:
        #print(alldata[line])
        linedata = alldata[line].split()
        #print(linedata)
        if len(linedata) < 1:
            break
        elif linedata[2] != ':':
            #Coulombic[(linedata[0].split('--')[0], linedata[0].split('--')[-1])] = 'None'
            line +=1
        else:
            Coulombic[(linedata[0].split('--')[0], linedata[0].split('--')[-1], linedata[1])] = linedata[5]
            line +=1
        linedata = alldata[line].split()
        #print(linedata)
        line +=1
        if len(linedata) < 1:
            break
        elif linedata[2] != ':':
            #NonCoul[(linedata[0].split('--')[0], linedata[0].split('--')[-1])] = 'None'
            pass
        else:
            NonCoul[(linedata[0].split('--')[0], linedata[0].split('--')[-1], linedata[1])] = linedata[5]
    #print(Coulombic)
    #print(NonCoul)
    return Coulombic, NonCoul, line 

def GetLoading(line):
    linedata = []
    data = {}
    while True:
        linedata = alldata[line].split()
        #print(linedata)
        if linedata[0] == '##':
            #print('I\'m done now!')
            break
        data[linedata[0]] = linedata[2]
        line +=1
    datalist = np.zeros(len(data))
    for i, key in enumerate(data):
        datalist[i] = data[key]
#    [t0, g, Neff_max] = timeseries.detectEquilibration(datalist)
#    ##########print('Loading autocorrelation analysis:\nt0 = {0}, g = {1}, N_max = {2}'.format(t0, g, Neff_max))
#    #print('t0 = {0}, g = {1}, N_max = {2}'.format(t0, g, Neff_max))
#    aveData = np.mean(datalist[t0::math.ceil(g)])
#    stdData = np.std(datalist[t0::math.ceil(g)])
    finval = float(datalist[-1])
    accept = []
    N_accept = 1
    accept.append(finval)
    for i in datalist[-2::-1]:
        #print((float(i)-finval)/finval)
        #print(0.05*finval)
        if abs((float(i)-finval)) < abs(0.01*finval):
            accept.append(i)
            N_accept +=1
        else:
            break
    ##########print(N_accept)
    CrudeAvg=np.mean(accept)
    CrudeStd = np.std(accept)
    ###########print('autocorrelation -> {0:.0f} pm {1:.0f}, crude -> {2:.0f} pm {3:.0f}'.format(aveData, stdData, CrudeAvg, CrudeStd))
    #points = []
    #for i in datalist[t0::math.ceil(g)]:
    #    points.append(i)
    #print('points = {0}'.format(len(points)))
    return CrudeAvg, CrudeStd, line

def ParseStep(line):
    t0, g, Neff, ave1, std1, line = GetEnergy(line)
    line += 4
    Coulombic, NonCoul, line = GetInteraction(line)
    while True:
        if alldata[line] == '':
            break
        line +=1
    line +=5
    aveA, stdA, line = GetLoading(line)
    return ave1, std1, Coulombic, NonCoul, aveA, stdA, t0, g, Neff

def FindStep(thing):
    #print(thing)
    instances = []
    for i, entry in enumerate(alldata):
        if thing in entry:
            i +=2
            instances.append(i)
    #print('Confile {0} appears on the lines: {1}'.format(thing, instances))
    return instances[0]
################################################################################################

EnergiesAve = {}
EnergiesStd = {}
Coul = {}
Coul_FluFlu = {}
Coul_FluFra = {}
NC_FluFlu = {}
NC_FluFra = {}
NonCoul = {}
Loadingave = {}
Loadingstd = {}
burnin = {}
efficiency = {}
NEff = {}

EnergiesAve[0] = []
EnergiesStd[0] = []
Coul[0] = []
NonCoul[0] = []
Loadingave[0] = []
Loadingstd[0] = []
burnin[0] = []
efficiency[0] = []
NEff[0] = []
Coul_FluFlu[0] = []
Coul_FluFra[0] = []
NC_FluFlu[0] = []
NC_FluFra[0] = []


for filename in os.listdir(directory):
    print(filename)
    startlines = []
    alldata = []
    confiles = set()
    if filename.endswith('postfile'):
        print('your filename is ', filename)
        with open('{0}/{1}'.format(directory, filename), 'r') as f:
        #alldata = f.readlines()
            for line in f:
                alldata.append(line.strip())
                if Template in line:
                    confiles.add(line.split(' ')[-1].split('.')[-1])


    for i in range (1, len(confiles)+1):
        #print(i)
        #print('Data for pressure point {0} begins at line {1}'.format(i, FindStep('{1}.{0}'.format(i, Template))))
        startlines.append(FindStep('{1}.{0}'.format(i, Template)))
    print(startlines)

    for n, i in enumerate(startlines, 1):
        if n not in EnergiesAve:
            EnergiesAve[n] = []
        if n not in EnergiesStd:
            EnergiesStd[n] = []
        if n not in Coul:
            Coul[n] = []
        if n not in NonCoul:
            NonCoul[n] = []
        if n not in Loadingave:
            Loadingave[n] = []
        if n not in Loadingstd:
            Loadingstd[n] = []
        if n not in burnin:
            burnin[n] = []
        if n not in efficiency:
            efficiency[n] = []
        if n not in NEff:
            NEff[n] = []
        if n not in Coul_FluFlu:
            Coul_FluFlu[n] = []
        if n not in Coul_FluFra:
            Coul_FluFra[n] = []
        if n not in NC_FluFlu:
            NC_FluFlu[n] = []
        if n not in NC_FluFra:
            NC_FluFra[n] = []
    EnergiesAve[0].append('{0} Energy ave (KJ)'.format(filename.split('.')[0]))
    EnergiesStd[0].append('{0} Energy stdev (KJ)'.format(filename.split('.')[0]))
    Coul[0].append('{0} Coulombic interactions'.format(filename.split('.')[0]))
    NonCoul[0].append('{0} Noncoulombic interactions'.format(filename.split('.')[0]))
    Loadingave[0].append('{0} ave {1} (mol/uc)'.format(filename.split('.')[0], species))
    Loadingstd[0].append('{0} stdev {1} (mol/uc)'.format(filename.split('.')[0], species))
    burnin[0].append('{0} burnin time'.format(filename.split('.')[0]))
    efficiency[0].append('{0} statistical inefficiency'.format(filename.split('.')[0]))
    NEff[0].append('{0} max samples'.format(filename.split('.')[0]))
    Coul_FluFlu[0].append('{0}'.format(filename.split('.')[0]))
    Coul_FluFra[0].append('{0}'.format(filename.split('.')[0]))
    NC_FluFlu[0].append('{0}'.format(filename.split('.')[0]))
    NC_FluFra[0].append('{0}'.format(filename.split('.')[0]))

    
    for n, i in enumerate(startlines, 1):
        #print ('Simulation step {0}'.format(n))
        a, b, c, d, e, f, i, j, k = ParseStep(i)
        EnergiesAve[n].append(round(a,3))
        EnergiesStd[n].append(round(b,3))
        Coul_FluFlu[n].append(c[('{0}'.format(species), '{0}'.format(species), 'Coulombic')])
        Coul_FluFra[n].append(c[('{0}'.format(species), 'IRMOF1', 'Coulombic')])
        NC_FluFlu[n].append(d[('{0}'.format(species), '{0}'.format(species), 'NonCoulom')])
        NC_FluFra[n].append(d[('{0}'.format(species), 'IRMOF1', 'NonCoulom')])
        #Coul[n].append(c)
        #NonCoul[n].append(d)
        print(e)
        Loadingave[n].append(round(e,3))
        Loadingstd[n].append(round(f,3))
        burnin[n].append(round(i,3))
        efficiency[n].append(round(j,3))
        NEff[n].append(round(k,3))

#print(Energies)
###################################################   
with open('EnergyTotal', 'w') as file:
    #file.write('Step, ')
    #file.write(', '.join(i for i in Energies['0']))
    #file.write('\n')
    for i, key in enumerate(EnergiesAve, 0):
        file.write('{0}, '.format(i))
        file.write(', '.join(str(j) for j in EnergiesAve[key]))
        file.write('\n')
with open('OccupancyTotal', 'w') as file:
    for i, key in enumerate(Loadingave, 0):
        file.write('{0}, '.format(i))
        file.write(', '.join(str(j) for j in Loadingave[key]))
        file.write('\n')
    file.write('\n')   
#print(Loading)
with open('NvsE', 'w') as file:
    for key in EnergiesAve:
        file.write('simulation step {0}\n'.format(key))
        file.write('{0}, {1}, {2}, {3}\n'.format(
                Loadingave[0][0], Loadingstd[0][0], EnergiesAve[0][0], EnergiesStd[0][0]
            ))
        for j in range(0,len(EnergiesAve[1])):
            file.write('{0}, {1}, {2}, {3}\n'.format(
                Loadingave[key][j], Loadingstd[key][j], EnergiesAve[key][j], EnergiesStd[key][j]
            ))
        file.write('\n')
with open('Energy stats', 'w') as file:
    for i, key in enumerate(burnin, 0):
        file.write('{0}, '.format(i))
        file.write(', '.join(str(j) for j in burnin[key]))
        file.write('\n')
    file.write('\n')
    for i, key in enumerate(efficiency, 0):
        file.write('{0}, '.format(i))
        file.write(', '.join(str(j) for j in efficiency[key]))
        file.write('\n')
    file.write('\n')   
    for i, key in enumerate(NEff, 0):
        file.write('{0}, '.format(i))
        file.write(', '.join(str(j) for j in NEff[key]))
        file.write('\n')
    file.write('\n')  

with open('specific interactions', 'w') as file:
    file.write('Fluid-fluid interactions\nNon-Coulombic\n')
    for i, key in enumerate(NC_FluFlu):
        file.write('{0}, '.format(i))
        file.write(', '.join(str(j) for j in NC_FluFlu[key]))
        file.write('\n')
    file.write('Coulombic\n')
    for i, key in enumerate(Coul_FluFlu):
        file.write('{0}, '.format(i))
        file.write(', '.join(str(j) for j in Coul_FluFlu[key]))
        file.write('\n')
    file.write('\n')
    file.write('Fluid-framework interactions\nNon-Coulombic\n')
    for i, key in enumerate(NC_FluFra):
        file.write('{0}, '.format(i))
        file.write(', '.join(str(j) for j in NC_FluFra[key]))
        file.write('\n')
    file.write('Coulombic\n')
    for i, key in enumerate(Coul_FluFra):
        file.write('{0}, '.format(i))
        file.write(', '.join(str(j) for j in Coul_FluFra[key]))
        file.write('\n')

NC_FluFlu_ave=[]
NC_FluFlu_std=[]
for i in range (1, len(confiles)+1):
    #print(NC_FluFlu[i])
    L = [float(n) for n in NC_FluFlu[i]]
    print(round(np.mean(L), 3))
    print(round(np.std(L), 3))
    NC_FluFlu_ave.append(round(np.mean(L), 3))
    NC_FluFlu_std.append(round(np.std(L), 3))
Coul_FluFlu_ave=[]
Coul_FluFlu_std=[]
for i in range (1, len(confiles)+1):
    #print(NC_FluFlu[i])
    L = [float(n) for n in Coul_FluFlu[i]]
    #print(round(np.mean(L), 3))
    #print(round(np.std(L), 3))
    Coul_FluFlu_ave.append(round(np.mean(L), 3))
    Coul_FluFlu_std.append(round(np.std(L), 3))
NC_FluFra_ave=[]
NC_FluFra_std=[]
for i in range (1, len(confiles)+1):
    #print(NC_FluFlu[i])
    L = [float(n) for n in NC_FluFra[i]]
    #print(round(np.mean(L), 3))
    #print(round(np.std(L), 3))
    NC_FluFra_ave.append(round(np.mean(L), 3))
    NC_FluFra_std.append(round(np.std(L), 3))
Coul_FluFra_ave=[]
Coul_FluFra_std=[]
for i in range (1, len(confiles)+1):
    #print(NC_FluFlu[i])
    L = [float(n) for n in Coul_FluFra[i]]
    #print(round(np.mean(L), 3))
    #print(round(np.std(L), 3))
    Coul_FluFra_ave.append(round(np.mean(L), 3))
    Coul_FluFra_std.append(round(np.std(L), 3))
print('--------------------------')
N_ave = []
N_std =[]
for i in range (1, len(confiles)+1):
    print(Loadingave[i])
    L = [float(n) for n in Loadingave[i]]
    print(L)
    #print(round(np.mean(L), 3))
    #print(round(np.std(L), 3))
    N_ave.append(round(np.mean(L), 3))
    N_std.append(round(np.std(L), 3))
    
    
with open('NvsInt', 'w') as file:
    file.write('\
Average loading (mol/uc), Stdev loading,\
average NC Fluid-fluid intereaction (KJ/mol), stdev NC lfuid-fluid interaction,\
average coul Fluid-fluid intereaction (KJ/mol), stdev coul fluid-fluid interaction,\
average NC Fluid-framework interaction (KJ/mol), stdev NC fluid-framework interaction,\
average coul Fluid-framework intereaction (KJ/mol), stdev coul fluid-framework interaction, sum NC\n')
    for i in range(len(confiles)):
        file.write('{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}\n'.format(
N_ave[i], N_std[i], NC_FluFlu_ave[i], NC_FluFlu_std[i], Coul_FluFlu_ave[i], Coul_FluFlu_std[i], NC_FluFra_ave[i], NC_FluFra_std[i], Coul_FluFra_ave[i], Coul_FluFra_std[i], NC_FluFlu_ave[i] +  NC_FluFra_ave[i]
    ))
    #file.write('Step, ')
    #file.write(', '.join(i for i in Energies['0']))
    #file.write('\n')
    
#    for i, key in enumerate(Energies, 0):
#        file.write('{0}, '.format(i))
#        file.write(', '.join(str(j) for j in Energies[key]))
#        file.write('\n')
#    file.write('\n')
#    for i, key in enumerate(DCMLoading, 0):
#        file.write('{0}, '.format(i))
#        file.write(', '.join(str(j) for j in DCMLoading[key]))
#        file.write('\n')
#    file.write('\n')   
#    for i, key in enumerate(ChlLoading, 0):
#        file.write('{0}, '.format(i))
#        file.write(', '.join(str(j) for j in ChlLoading[key]))
#        file.write('\n')
#    file.write('\n')  
#print(Energies)
#
#print(ChlLoading)
#for i in Coul[1:]:
#    pass
#pprint.pprint(NC_FluFlu)
#pprint.pprint(NC_FluFra)

