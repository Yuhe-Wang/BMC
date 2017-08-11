# Grammar instruction
# 1) All characters after # or ! or % are comments in each line
# 2) The configuration is made up of key - value pairs and sub blocks
# 3) Key - value pair format : key = value, if value consists of multiple data, they should have the same type,
#     and be segmented by space character.For example, position = 2.0 3.0 0.0
# 4) Sub block format :
#     blockname = # here mustn't follow any recordable string
#     {
#          block content, same grammar as parent content, i.e.potentially has sub blocks too.
#     }
# 5) It's allowed to have duplicated keys, which are recorded in their sequence appearing in file.
# 6) It's flexible to add empty line, space, tab, or other unprintable characters.
# 7) The file should be ended with an EMPTY line('\n').
# 8) Strongly recommend editing this file with Notepad++ if on windows platform. You can modify the extension
#    to ".py, .sh, .bsh .pl .rb" to utilize block comment based on #, 
#    to ".m" to utilize block comment based on %
#    to ".f90" to utilize block comment based on ! (not well supported)
#    and remember to associate this extension in Notepad++ with administrative right.

# note the directory is counted from where the executable file locates
# $cdir$ means current directory where this config file locates
# How many particles to be simulated
NSIMU = 1e6
# note this file name should not contain extension. It will automatically add ".dose"
output file name = $cdir$PENELOPE
log file directory = $cdir$
log file append = no # no to delete all previous log content if it detects the same file name
# a short description stated in the log header, recommended if you choose to append log content
log description = NA

# GPU thread division setting
GPU Query = no   # if yes, will print GPU card info before running
GPU =
{
	GPU Index = 0 1   #set multiple GPU index to use, starting from 0.
	GPU Block Num = 128
	GPU Block Dim = 256
	GPU Batch Num = 10
	GPU RNG Statistic = yes # find how often it naturally refill RNG buffer, output file name is rngstat.txt
	GPU Refill Period = 70 # a negative number will disable periodic refill of RNG
	Source Reuse Times = 10
}
# Note the GPU will copy (GPU Batch Num)* (Block Num)* (Block Dim) particles from Source Head
# in order to save copy time


# PENELOPE information
PENELOPE =
{
	material database = PENELOPE_DataBase.7z # Will automatically extract the material files from database	
        # if use last generated material data to boost up simulation.
	# Must turn it off during first run after modifying Emax or MAT#
	useLastMats = yes 
	# last dose file = $cdir$$fname$.dose # if it exists, we will try to reuse it
	#EMax = 1.35e6   #unit eV
	DSMax = 0.5  #unit cm
	
	particle stack depth = 40
	# possible event to switch off if you are sure their probability can be ignored
	# set 0 to switch off
	# Photon events: Rayleigh scattering, Compton scattering, Photoelectric effect, Pair production
	GRA = 1
	GCO = 1
	GPH = 1
	GPP = 1
	
	SEC = 1 # top switch: whether to simulate secondary particles
	# Electron events: Elastic scattering, Inelastic scattering, Bremsstrahlung radiation, Inner-shell impact ionisation
	EEL = 1
	EIN = 1
	EBR = 1
	ESI = 1
	
	POS = 1 # whether to simulate positron events, bigger switch
	# Positron events: Elastic scattering, Inelastic scattering, Inner-shell impact ionisation, Annihilation
	PEL = 1 
	PIN = 1 
	PSI = 1
	PAN = 1
	
	# The program will search material block name in sequence as MAT1, MAT2, MAT3....
	# and stop when failing to find MAT#.So you can put these blocks in any order, but
	# there mustn't be any discontinued integer #.

	# sub block describing the properties of material #1
	MAT =
	{
		ID number= 278    # PENELOPE material index, Water
		Eabs_e = 200.0e3  # unit eV
		Eabs_g = 50.0e3   # unit eV
		Eabs_p = 200.0e3  # unit eV
		C1 = 0.2
		C2 = 0.2
		Wcc = 200.0e3     # unit eV
		Wcr = 50.0e3	  # unit eV
	}
}


SOURCEHEAD =
{
	# phase space file name
	DataDir = ViewRayData/
	CoPhaseSpace = ViewRayCo.phsp
	Geometry file = vrSourceDC.py
	
	# thread num used to generate particle from source head
	NThread = 1
	Sample Stack Depth = 400
	
	# two kinds of beam files are supported: (1)mfQAData.beams (2)PlanOverview.txt
	Beams file = $cdir$mfQAData.beams
	
	# If Beams file doesn't exists, the following customized beams will be used
	#{{
	Leaf spacing adjustment = $LeafAdjustment$ # unit cm
	isocenter = 0 0 0 #unit cm, iso center's position relative to the patient phantom center
	# Define the beams.
	# This block can be repeated as many times as needed to define all beams.
	# In each beam block there can be an arbitrary number of segment blocks to define the segments of the beam.
	
	Beam =
	{
		gantry angle = 90  #in degree
		
		head index = 1
		# the segment search will be like segment1, segment 2, segment 3, ...
		# the key should contain continued integer.
		segment =
		{
			on time = 10  #uint second
			leaf position = 0 29  0  0
			leaf position = 13 16 -2.10 2.10 #the later one will overlap part of previous defined
		}
	}
	
	#}}
}


# phantom information, currently we only deal with cuboid
PHANTOM =
{	
	# if define phantomName, the next block is ignored. It can be .phtm or .red format
	phantomName = $cdir$mfQAData.red
	
	# if the voxel density is less than threshold, the dose will be trimmed to zero
	trim dose threshold = 0.08 # unit g/cm^3
	
	# If phantomName is not defined, it'll use following parameters to define the cuboid instead.
	#{{
	DX DY DZ = 0.3 0.3 0.3 #unit cm
	
	# the later will repalce the overlapped part in the previous one
	# the first one should be larger than any others
	# meaning: NX, NY, NZ, x,y,z, density
	# x,y,z is the center coordinates of cuboid in isocenter coordinate system
	matrix = 201 201 201 0 0 0 0.002
	matrix = 101 27 101 0 0 0 1.0
	
	#if the variable SSD exist, it will redefine the position in z direction
	# SSD = $SSD$ #unit cm
	#if the variable centerXZ exist, it will redefine the position in XZ plane
	centerXZ = 0 0 #unit cm, coordinates of the center of the cuboid relative to isocenter
	# MatID = 1 #index of material listed above
	#}}
	
	#this is necessary if using RED file as input
	#couch density factor = 1.345 # coefficient to modify the density in couchConfig file
	# only replace couch density when couchConfig exists
	#couchConfig = couch.ced
	# if couchXYZ doesn't exist, use position definition in ViewRay's plan 
	couchXYZ = 0 4.05 0 # the top surface center's coordinates in isocenter coordinate system
	
	magnetic field switch    = off #other from "on" will turn B off
	magnetic field strength  = 0.35 #unit Tesla
	# direction can be specified by theta, phi or unit vector such as 0 0 1
	# the parameter number determines the type
	magnetic field direction = 0 0 1
}

