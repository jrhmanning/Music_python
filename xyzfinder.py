import shutil

path = "../../../../../03_music_chloroform_forcefield/experiments/DCM/298"
for i in range (1,17):
    shutil.copyfile('{0}/{1:02d}/{1:02d}.IRMOF1.25.729kpa.xyz'.format(path, i), 'DCM.298.{0:02d}.25.729.xyz'.format(i))
    print('{0:02d}'.format(i))
