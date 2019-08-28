######################################################################################################################################################
##Hello! I'm a program who's job it is to make producing multiple music isotherms easier!
##I eat a set of working music files from a repository (the directory has to be called "./experiments/repository/") and your desired sorbent, temperature and number of points
##I copy these files into a new directory (./experiments/*your sorbent here*/*your temperature here*/) set up to run 16 parallel simulations of the same system using taskfarmer
##I won't run it for you (right now) because my author is a scaredypants, but We'll get there eventuallly, I'm sure
######################################################################################################################################################
##The directory tree I'll expect will go something like this:
### "Music"
###  | -- "atoms"
###  | -- "molecules"
###  | -- "maps"
###  | -- "mapgen"
###  |
###  | -- results
###  \ -- "experiments"
###          | -- Me!
###          | -- "repository"
###          |      | -- run.taskfarmer
###          |      | -- run.gcmc
###          |      | -- gcmc.ctr
###          |      | -- post.ctr
###          |      | -- intramolecular_file
###          |      | -- atom_atom_file
###          |      \ -- sorb_sorb_file
###          \ -- *your sorbent name*
###                 \ *your simulation temperature*
###						| -- run.taskfarmer
###						| -- taskfarm.list
###						\ -- 01-16
###							   | -- As repository, plus:
###							   | -- pressure.*your sorbent*.*temperature*.dat
###							   \ -- setpath
######################################################################################################################################################
import os
import errno
import random
import logisothermgenerator as ig  ## I wrote a separate isotherm writer script, containing its own Antoine parameters library and the functions "pSat", "isothermcalculator", amd "Pressurefilewriter"
Pressure = 0
isotherm = []
######################################################################################################################################################
####The file copying function definitions
######################################################################################################################################################
def directorymaker(directory, dirout):
    filename = "%s/%02d/test.txt" % (dirout, directory)
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    with open(filename, "w") as f:
        f.write("FOOBAR")
    print("Directory written!")
		
def GcmcControlChanger(Species, T, n, directory, dirout):
    # Read in the file
    with open("./repository/gcmc.ctr.template", 'r') as file :
        filedata = file.read()
    x = random.randint(0,99999)
    # Replace the target string
    filedata = filedata.replace('*your_temperature_here*', str(T))
    filedata = filedata.replace('*your_number_of_points_here*', str(n))
    filedata = filedata.replace('*your_Species_here*', Species)
    filedata = filedata.replace('*your_random_seed*', "%05d" % x)
    filedata = filedata.replace('*your_pressure_file_here*', "pressure.%s.%s.dat" % (Species, T))
        # Write the file out again
    fileout = ("%s/%02d/gcmc.ctr" % (dirout, directory))
    with open(fileout, 'w') as file:
        file.write(filedata)
    print("Control file written!")

def PostControlChanger(Species, T, n, directory, dirout):
    # Read in the file
    with open("./repository/post.ctr.template", 'r') as file :
        filedata = file.read()
    # Replace the target string
    filedata = filedata.replace('*your_temperature_here*', str(T))
    filedata = filedata.replace('*your_number_of_points_here*', str(n))
    filedata = filedata.replace('*your_Species_here*', Species)

        # Write the file out again
    fileout = ("%s/%02d/post.ctr" % (dirout, directory))
    with open(fileout, 'w') as file:
        file.write(filedata)
    print("Post control file written!")

def GcmcRunWriter(Species, T, directory, dirout):
    # Read in the file
    with open("./repository/run.gcmc.template", 'r') as file :
        filedata = file.read()
    # Replace the target string
    filedata = filedata.replace('*your_temperature_here*', str(T))
    filedata = filedata.replace('*your_Species_here*', Species)
    filedata = filedata.replace('*your_tree_here*', "%s/%02d/" % (dirout, directory))
    # write the file out again
    with open("%s/%02d/run.gcmc" % (dirout, directory), 'w') as file:
        file.write(filedata)
	os.chmod("%s/%02d/run.gcmc" % (dirout, directory), 0o777) #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    print("Run file written!")

def TaskfarmRunWriter(Species, T, dirout):
    # Read in the file
    with open("./repository/run.taskfarmer.template", 'r') as file :
        filedata = file.read()
    # Replace the target string
    filedata = filedata.replace('*your_temperature_here*', str(T))
    filedata = filedata.replace('*your_Species_here*', Species)
        # Write the file out again
    filedata = filedata.replace('*your_tree_here*', dirout)
    with open("%s/run.taskfarmer" %(dirout), 'w') as file:
        file.write(filedata)
	os.chmod("%s/run.taskfarmer" %(dirout), 0o777)   #chmod so it's fully rwx for everyone (0o denotes the following number is in octal)
    print("Taskfarmer file written!")

def InteractionMover(Species, directory, dirout, version = '1'):
	file1 = "atom_atom_file"
	file2 = "sorb_sorb_file"
	file3 = "intramolecular_file"
	with open("./repository/{0}_v{1}.template" .format(file1, version), 'r') as file :
 	       filedata1 = file.read()
	with open("./repository/{0}_v{1}.template" .format(file2, version), 'r') as file :
		filedata2 = file.read()
	with open("./repository/%s.template" %(file3), 'r') as file :
	        filedata3 = file.read()
	filedata2 = filedata2.replace('*your_Species_here*', Species)
	filedata3 = filedata3.replace('*your_Species_here*', Species)
	with open("%s/%02d/%s" %(dirout, directory, file1), 'w') as file:
        	file.write(filedata1)
	with open("%s/%02d/%s" %(dirout, directory, file2), 'w') as file:
        	file.write(filedata2)
	with open("%s/%02d/%s" %(dirout, directory, file3), 'w') as file:
        	file.write(filedata3)
	print("Interactions files written!")


def TaskfarmMover(dirout):
	file1 = "taskfarm"
	with open("./repository/%s.template" %(file1), 'r') as file :
        	filedata = file.read()
	with open("%s/%s" %(dirout, file1), 'w') as file:
        	file.write(filedata)
	with open("%s/%s.backup" %(dirout, file1), 'w') as file:
        	file.write(filedata)	
	print("Taskfarmer joblist written!")

def IsothermExtractMover(Species, T, dirout):
	file1 = "isothermextractor.py"
	with open("./repository/%s.template" %(file1), 'r') as file :
	        filedata = file.read()
	filedata = filedata.replace('*your_temperature_here*', str(T))
	filedata = filedata.replace('*your_species_here*', Species)
	with open("%s/%s" %(dirout, file1), 'w') as file:
	        file.write(filedata)	
	print("Python isotherm extractor written!")

######################################################################################################################################################
####The program
######################################################################################################################################################


random.seed()                                             #This is the random seed for putting a number in the control file

#####User input section
Species = input("What is your sorbent Species called?")
Species = str(Species)
T = eval(input("What temperature (in K) would you like to test?"))
T = int(T)
n = eval(input("How many isotherm points do you want?"))
n= int(n)
#Or you could hard code any part of it:
#Species = str("CCL4")
#T = 298K
#n = 20

version = Species.split('_v')[1]
print(version)

#####Pressure value generation
satP = ig.pSat(Species, T)
istm = ig.isothermcalculator(satP, n)
print(istm)

#####file writing section
dirout = "./%s/%s" % (Species, T)                      #Sets the directory you're aiming at
for directory in range(1,17):
	print("Working in directory %02d:" % (directory))
	directorymaker(directory, dirout)                                               #makes the overall experiment directory, and 16 repetition subdirectories
	GcmcControlChanger(Species, T, n, directory, dirout)                            #puts in the gcmc control file
	PostControlChanger(Species, T, n, directory, dirout)                            #puts in the post control file
	GcmcRunWriter(Species, T, directory, dirout)							                    #puts in the runfile
	InteractionMover(Species, directory, dirout, version)                                             #puts in the 3 interaction files
	print("Pressure check %s" % (satP))
	ig.PressureFileWriter(Species, T, satP, istm, directory, dirout)         #Writes the new pressure file into this directory based on the output of isothermgenerator.isothermcalculator
	

print("Working in the general directory:")
TaskfarmRunWriter(Species, T, dirout)                                                          #Writes the taskfarmer run file
TaskfarmMover(dirout)                                                              #writes the taskfarmer task list
IsothermExtractMover(Species, T, dirout)	                                                           #Writes isothermextractor.py into place
