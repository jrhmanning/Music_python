species0 = 'DMF'
species1 = 'Chloroform'
n = 10
p0 = 0.523

def mixisotherm(species0, species1, p0, n=10):
    iso0 = []
    iso1 = []
    for i in range (0, n):
        iso0.append(round(float(p0)*float(1.-float(i)/(n-1)), 3))
        iso1.append(round(float(p0)*float(float(i)/(n-1)), 3))
    print('{0} isotherm: {1}\n{2} isotherm: {3}'.format(species0, iso0, species1, iso1))
    with open('pressure.dat', 'w') as f:
        f.write('{0}\n'.format(species0))
        f.write('{0}\n'.format(n))
        f.write('{0}\n'.format(', '.join(str(thing) for thing in iso0)))
        f.write('{0}\n'.format(species1))
        f.write('{0}\n'.format(n))
        f.write('{0}'.format(', '.join(str(thing) for thing in iso1)))
    print('Exchange isotherm file written!')

mixisotherm(species0, species1, p0, n)
