# Music_python
This git repository is full of python tools for using Music!
I'll try to improve the documentation soon, but for now, here are the headlines:

-------Repository Structure---------
-
-------The foundations---------
-
The repository is made up of a main python script directory (Hi!), a directory called musicpy where all the useful bits are hidden, and a directory full of old spaghetti code for reference/cannibalising.

In the 'musicpy' directory there are 4 .py files: Antoine, Forcefield, setup, and analysis.

Antoine contains one dictionary of Antoine equation parameters extracted from the NIST online database, for compounds I've used in the past. This provides an estimate of sorbent vapour pressure at a certain temperature, allowing you to calculate fugacities relative to this quantity for making isotherms.
Forcefield is similar, but contains two dictionsaries. The first - 'IntParams' - is a dictionary of forcefield parameters for different atom types in my models. The keys of this dictionary take the format A_B, where A is the atom type (e.g. H) and B is the forcefield used (e.g. TRAPPE). Also in this file is a dictionary called MOF_el_list, which contains lists of element types within pure MOF structures.
These two are expected to be living documents that get added to and pruned as and when is required for your specific simulations.

setup.py and analysis.py contain a list of functions for simulation setup and analysis, respectively. (As of 29-August-2019, analysis is only half-finished, sorry!) Using these, you can go from essentially nothing to a working simulaiton in very little time at all, and postprocess the ensuing results in a way to get immediately graphable data from the simulation outputs. 
These files are intended not to be used directly, but instead to be invoked by other scripts to do some (or all) of the simulation manuipulation.

-------The simulation setup scripts---------
-
In the main directory, there are a group of numbered python scripts from 00 to 11, which are intended to be (roughly) followed in order to make a quick simulation and grab an isotherm from it. These are described below:

00 - molfilemaker.py: turns a .pdb of your sorbent into a basic Music .mol file.

01 - atom_atom_maker.py: using the new .mol file, some of your sorbate information, and the Forcefield database, creates an atom_atom_interaction file. By default (29-aug-19), it does LB mixing for your lennard-jones interactions and Ewald summation for sorbent-sorbent coulombic interactions.

02 - sorb_sorb_writer.py: similar to the above, this file takes in the name of your sorbent and sorbate, and the atoms involved, and creates a molecule-molecule interactions file for the above. By default (29-aug-19), it does LB mixing for all lennard-jones interactions, ewald summation for all coulombic interactions, and expects to use pmap and emap files of precalculated interations for all sorbent-sorbate interactions (lennard-jones and coulombic).

03 - Intrawriter.py: this file writes an intramolecular interactions file for the above, including all sorbent atoms (useful for map generation!). By default (29-aug-19) it lists no intramolecular flags for anything.

04 - MapCtrlWriter.py: this file writes Music control scripts for making all the pmaps in your sorbent molecules, to be used with music_mapmaker. You could probably repurpose this to also make emap control files too, but there are some ongoing issues (29-aug-19) with making emaps on Balena.

05 - MapMakeWriter.py: this file makes a bash script for running music_mapmaker for all of the above map control files, one after the other. 

06 - IsothermGenerator.py: this file takes in Antoine data from Antoine.py, estimates the saturation pressure (pSat), and creates a log-linear sequence of fugacities for Music to work through. Because only a short section of the isotherm is usually interesting, it has maximum and minimum relative pressure variables (equal to pSat\*10^(-x)) to ensure you're simulating only the good bits.

07 - GCMCCtrlWriter.py: this file takes in some useful information about your simulation, sorbent, and isotherm, and creates a control file for running your music_gcmc simulation. Optionally, it can create a series of restartfiles too to give you a final snapshot from your simulation results.

08 - PostCtrlWriter.py: this file does much the same as 07, but creates control files for music_post analysis of your final .con files producecd using music_gcmc.

09 - GCMCRunWriter.py: this file writes a bash script for executing the above music_gccm and music_post instances. This can then be run  to produce a single isotherm directly.

The following two scripts are only useful if you want to run multiple parallel simulations using taskfarmer. The University of Bath's HPC, Balena, assigns a minimum of a complete node (16 cores) per job submitted, therefore it's very wasteful to run Music for a single job at once. By default (29-aug-19), the following scripts are useful for running repeats of the above to reduce uncertainty per job.

10 - IsothermExtractWriter.py: this python script writes a second python script, whose job is to read all of the individual isotherms from each and present the mean and standard deviation of all the above files as a .csv file for graphing.

11 - taskfarmRunWriter.py: this file creates a bash script for running multiple Music simulations in parallel, using taskfarmer.

-------A (very short) quickstart guide---------
-
In Music, the majority of the differnt information are held in separate files, and all brought together the for the full simulation.

Generally, the files required for Music are: 

1. Atom and molecule files for each of the species you're simulating (including dummy atoms)

2. Interactions databases describing all the possible interactions at an atom-atom scale, molecule-molecule scale, and intramolecular

3. Map files of pre-calculated interaction strengths (i.e. between a sorbent and fixed sorbate)

4. The simulation control files, including any pressure point files and postprocessing control files too

The python scripts in this repository run through your setup in roughly this order. First, manually make all the atom files you need, including for any dummy atoms you'll be using . Second, make sure you have a .mol file of your sorbate and .pdb file(s) for your sorbent(s). Run script 00 to make a basic .mol file of your sorbent(s). This should get you to the end of (1). 

Make sure that all of the key interactions from your forcefield and your framework atomtype list are in Forcefield.py, then manually make .mol files for each of the atoms you'll later want to make maps with (generally all the atoms in yoru sorbent). Once done, run through python scripts 01-03 to produce all the interactions files you'll be needing. Since the species variable is a list, you can (probably - untested as of 29-aug-19) make a single set of interactions files that'll both amke your maps and run your simulations. This gets you to the end of (2). 

Run through scripts 04-05 to make control and runfiles for your mapmaking. Then, arranging the above documents in the right directories (I'll put an expected file tree in the resources directory), run the runscript on your HPC. Expect each map to take ~10 minutes to run, so use an interactive node or submit it properly! Once done, you should be at the end of (3).

Finally, run through scripts 06-09 to make the appropriate runscripts and control files for making and analysing your isotherm. Once all files are in the correct place (see attached tree), you'll be ready to go! 

If you want to run through a taskfarmer to maximise your HPC resources, run through scripts 10-11 to make your taskfarmer run file and python isotherm extractor, then you'll be ready to run your Music simulations in parallel!

-------Some tips and tricks---------
-
Using He as your dummy atom element type makes it a lot easier to visualise later on using VMD etc.
