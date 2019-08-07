######This python script takes in some simulation information, and spits out a music_post control file to get an isotherm
#It requires the following variables: 
#species - the name of your sorbent species. Since we're making a map, this is actually one of your sorbent elements
#framework - the name of your MOF (or other porous material
#T - your simulaiton temperature
#n - your number of pressure points
#################
from math import sqrt, log
import logging
from pathlib import Path
import datetime
now = datetime.datetime.now()
import random
#################
####This section sorts out your messages from this script to the console and a log file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler('{0}/sorb_sorb_maker.log'.format('.'))
handler.setLevel(logging.DEBUG)
handler2 = logging.StreamHandler()
handler2.setLevel(logging.INFO)
logger.addHandler(handler)
logger.addHandler(handler2)
logger.debug('-----------------------------------')
logger.debug('Generating new files on {0}-{1}-{2} at {3}:{4}'.format(now.year, now.month, now.day, now.hour, now.minute))
logger.debug('-----------------------------------')
##################
species = ['DMF_v4']
framework = 'IRMOF1'
#T = 298
n = 20
##################
#Writes a .ctr file for your postprocessing. There;s orobably a better way than using music_post, but I'm not there yet
def PostControlChanger(Species,  n, framework, dirout, directory):
    with open("{0}post.ctr" .format(dirout), 'w') as file:
        file.write("""####This section is apparently required for working with any post code. Who knows why? not me!
#
#
------------------------------------------------------------
   ### Required section ######
-- Post Processor Information ------------
GCMC                            # Type of simulation GCMC, NVTMC , MD ....
./{1}.{0}.con                    # basename for config files
1, {2}                          # first and last file numbers
post.ctr.out                       # name for new ctrlfile that will regenerated
{3:02d}.postfile          # Base name for output files
0, 0                         # Percentages of data to skipped at start and end 


# The sections below are necessary only if you want the corresponding 
# analysis performed
# ---------------- ALL OF THEM ARE OPTIONAL ------------------------


####    This section is reqd for energy averages in your post code output files
####    as of now only total enrgies vs sim. step
------ Post : Energy Average Info -----------------------------------
100       # Number of blocks into which data should be divided for stats

####    This section is reqd for Loading averages in your post code outputfiles
####    as of now only species loading vs sim. step (for all species)
------ Post : Loading Average Info -----------------------------------
100       # Number of blocks into which data should be divided for stats""".format(Species, framework, n, directory))
    logger.debug("Post control file written!")
#############
for directory in range(1, 17):
    PostControlChanger(species[-1],  n, framework, './{0:02d}/'.format(directory), directory)
