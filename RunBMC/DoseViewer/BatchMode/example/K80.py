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
# 8) Strongly recommend editing this file with Notepad++ if on windows platform.

# Note the working directory is where the executable file locates(../../) rather than where this config file locates
# $cdir$ means the directory where this config file locates
# Use $Var$ to acess the marco Var defined in ../../jobList.py

# How many particles to be simulated
NSIMU = 1e7
target uncertainty = -0.01 # This option will override NSIMU. A negtive value will cancel this option.
threshold region = 0.10 # The percent of max dose as the threshold to calculate uncertainty
# note this file name should not contain extension. It will automatically add ".dose"
output file name = $cdir$PENELOPE
ViewRay format = no # whether to output the data in format of ViewRay
proceed last simulation = no # whether to continue last simulation with existing dose file
peek output num = 0 # how many intermediate dose files to output. Useful if the simulation takes too long.
execution type = 0 # 0 means kernel-by-kernel method(default), and 1 means one-thread-one-history method
log file directory = $cdir$
log file append = no # no to delete all previous log content if it detects the same file name
log short mode = no # short mode means to show the simulation process in one line
# a short description stated in the log header, recommended if you choose to append log content
log description = NA


# GPU setting
GPU =
{
	GPU Query = no   # if yes, will print GPU card info before running
	#GPU Index = 0 1   #set multiple GPU index to use, starting from 0.
	GPU Num = 2 # if "GPU Index" is not listed, it will try to select "GPU Num" GPU(s) to run
	GPU Block Num = 128
	GPU Block Dim = 256
	GPU Batch Num = 10
	Source Reuse Times = -1 # set -1 to determine it automatically
}

# PENELOPE parameters. This config file can be shared with gPENELOPE
PENELOPE =
{
	material database = C:\PENELOPE_DataBase.7z # Will automatically extract the material files from database
	# Whether to use last generated material data to boost up simulation.
	useLastMats = yes # Suggest you to keep yes. Delete common.dat to enforce no once if you modified MAT block(s) 
	#EMax = 1.35e6   #unit eV, if it's specified, it will override the value returned from source head
	DSMax = 0.5  #unit cm
	
	particle stack depth = 40
	# Frequency to simulate the following events. Set 0 to turn off
	# Photon events: Rayleigh scattering, Compton scattering, Photoelectric effect, Pair production
	GRA = 1
	GCO = 1
	GPH = 1
	GPP = 1
	
	SEC = 1 # Top switch: whether to simulate secondary particles
	# Electron events: Elastic scattering, Inelastic scattering, Bremsstrahlung radiation, Inner-shell impact ionisation
	EEL = 1
	EIN = 1
	EBR = 1
	ESI = 1
	
	POS = 1 # Bigger switch: whether to simulate positron events
	# Positron events: Elastic scattering, Inelastic scattering, Inner-shell impact ionisation, Annihilation
	PEL = 1 
	PIN = 1 
	PSI = 1
	PAN = 1
	
	# Sub block describing the simulation parameters of each material
	# You can repeat this block as many as necessary, and the first material will be the default material
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
	MAT =
	{
		ID number= 74    # PENELOPE material index, Tungsten
		Eabs_e = 200.0e3  #unit eV
		Eabs_g = 50.0e3  #unit eV
		Eabs_p = 200.0e3  #unit eV
		C1 = 0.2
		C2 = 0.2
		Wcc = 200.0e3     #unit eV
		Wcr = 50.0e3	    #unit eV
	}
}

#Zeus parameters. This config file can be shared with gZeus
ZEUS =
{
	Simulate electron = yes # yes or no
	Forward detection = yes # whether to freely advance the photon to high density area (skip the air)
	Particle stack depth = 40 # Stack depth for each GPU thread. It may not be used in some versions
	NMaxSplit = 1 # default 50. Set it to 1 to cancel splitting
	Fixed Split = yes # yes or no. If no, the splitting number will vary with particles' energy
	EAbsPhoton = 50e3 #unit eV. Default 50e3
	EAbsElectron = 50e3 #unit eV. Default 50e3
	EMaxCSDA = 200e3 # unit eV. Default 200e3
}

SOURCEHEAD =
{
	# two kinds of beam files are supported: (1)mfQAData.beams (2)PlanOverview.txt
	Beams file = $cdir$mfQAData.beams
	
	# Note the energy range of Co60 is < 1.35 MeV (1.17 MeV and 1.33 MeV primary decay)
	Energy scale factor = 1.0 # Simulate higher or lower energy by simply scaling the output energy
	
	DataDir = ViewRayData/
	CoPhaseSpace = ViewRayCo.phsp # phase space file name. It is given in a seperate path due to its big size
	
	# thread num used to generate particle from source head
	NThread = 16
	Sample Stack Depth = 400
	
	overwrite isocenter = no # Change the isocenter position using following definition in above beam file

	Leaf spacing adjustment = $LeafAdjustment$ # unit cm
	
	# If Beams file doesn't exists, the following customized beams will be used
	#{{
	
	isocenter = 0 0 0 #unit cm, iso center's position relative to the patient phantom center
	# Define the beams.
	# This block can be repeated as many times as needed to define all beams.
	# In each beam block there can be an arbitrary number of segment blocks to define the segments of the beam.
	
	Beam =
	{
		gantry angle = 0  #in degree
		head index = 1
		# the segment search will be like segment1, segment 2, segment 3, ...
		# the key should contain continued integer.
		segment =
		{
			on time = 100  #uint second
			leaf position = 0 29  0  0
			leaf position = 13 16 -2.10 2.10 #the later one will overlap part of previous defined
		}
	}
	
	#}}
}

# phantom information, currently we only deal with cuboid voxel matrix
PHANTOM =
{	
	# if define phantomName, the next block is ignored. It can be .phtm or .red format
	phantomName = $cdir$mfQAData.red
	
	# if the voxel density is less than threshold, the dose will be trimmed to zero
	trim dose threshold = 0.004 # unit g/cm^3. Note air density = 0.00129
	
	# If phantomName is not defined, it'll use following parameters to define the cuboid instead.
	#{{
	DX DY DZ = 0.3 0.3 0.3 #unit cm
	
	# the later will repalce the overlapped part in the previous one
	# the first one should be larger than any others
	# meaning: NX, NY, NZ, x,y,z, relative density(optional, default 1), PENELOPE mat index(optional, default 278 water)
	# x,y,z is the center coordinates of cuboid in isocenter coordinate system
	matrix = 101 101 101 0 0 0 1.0 # the water
	matrix = 101 33 101 0 0 0 0.3 # the lung
	
	#if the variable SSD exist, it will redefine the position in z direction
	SSD = $SSD$ #unit cm
	#if the variable centerXZ exist, it will redefine the position in XZ plane
	centerXZ = 0 0 #unit cm, coordinates of the center of the cuboid relative to isocenter
	#}}
	
	# This is a necessary correction if the RED file hasn't been scaled
	#couch density factor = 1.345 # coefficient to modify the density in couchConfig file
	# It will replace couch density only when couchConfig exists
	#couchConfig = couch.ced
	# if couchXYZ doesn't exist, use position definition in ViewRay's plan 
	couchXYZ = 0 4.05 0 # the top surface center's coordinates in isocenter coordinate system
	
	magnetic field strength  = 0 #unit Tesla
	# direction can be specified by theta, phi (e.g. 90 180) or unit vector (e.g. 0 0 -1)
	magnetic field direction = 0 0 1
}

OCTREE =
{
	go through = no # whether to go through the octree
	3D mesh file = $cdir$block.obj
	max level = 9 # Max level of the octree. It should be no more than the internal MAX_LEVEL
	max triangles = 6 # Max number of triangles in each spacial box
}
