# Grammar instruction:
# 1) All characters after # or ! or % are comments in each line
# 2) The configuration is made up of key - value pairs and sub blocks
# 3) Key - value pair format : key = value, if value consists of multiple data, they should have the same type (or can
#    be casted to the same type), and be segmented by space character. For example, position = 2.0 3.0 0.0
# 4) Sub block format :
#     blockname = # here must NOT be followed by any recordable string
#     {
#          block content, same grammar as parent content, i.e.potentially has sub blocks too.
#     }
# 5) It's allowed to have duplicated keys, which are recorded in their sequence appearing in file.
# 6) It's flexible to add empty line, space, tab, or other unprintable characters.
# 7) The file should be ended with an EMPTY line('\n').
# 8) Strongly recommend editing this file with Notepad++ if working on windows platform. You can modify the extension
#    to ".py, .sh, .bsh .pl .rb" to utilize block comment based on #, 
#    to ".m" to utilize block comment based on %
#    to ".f90" to utilize block comment based on ! (not well supported)
#    and remember to associate this extension in Notepad++ with administrative right.

# Note the relative directory is calculated from where the executable file locates.
# $cdir$ means the current directory where this config file locates.
# You can also access the macro definitions in jobList.py by using $variable$

# How many particles to be simulated
NSIMU = 1e7
target uncertainty = -1 #0.01 # This option will override NSIMU. A negtive value will cancel this option.
threshold region = 0.50 # The percent of max dose as the threshold to calculate uncertainty
# Note the output file name should not contain extension. It will automatically add ".dose"
output file name = $cdir$gZeus
ViewRay format = no
proceed last simulation = no # whether to continue last simulation with existing dose file
log file directory = $cdir$
log file append = no # yes or no. Set no to delete all previous log content if it detects the same file name
log short mode = no # short mode means to show the simulation process in one line
# A short description stated in the log header, recommended if you choose to append log content
log description = 4.2cm square field

# GPU thread division setting
GPU =
{
	GPU Query = no # yes or no. If yes, it will print the GPU information and pause.
	#GPU Index = 0 1   # Set multiple GPU indices to use, starting from 0.
	GPU Num = 2 # if "GPU Index" is not listed, it will try to select "GPU Num" GPU(s) to run
	GPU Block Num = 64
	GPU Block Dim = 64
	GPU Batch Num = 16 # Number of incident photon each thread processes
	Source Reuse Times = 1 # Default to 1. If the source generating speed lags the GPU simulation, increase it 
}
# Note the GPU will copy (GPU Batch Num)* (Block Num)* (Block Dim) particles from Source Head in order to save copy time

#Zeus parameters
ZEUS =
{
	Simulate electron = yes #yes # yes or no
	Forward detection = yes # whether to freely advance the photon to high density area (skip the air)
	Particle stack depth = 40 # Stack depth for each GPU thread. It may not be used in some versions
	NMaxSplit = 50 # default 50. Set it to 1 to cancel splitting
	Fixed Split = yes # yes or no. If no, the splitting number will vary with particles' energy
	EAbsPhoton = 50e3 #unit eV. Default 50e3
	EAbsElectron = 50e3 #unit eV. Default 50e3
	EMaxCSDA = 200e3 # unit eV. Default 200e3
}

SOURCEHEAD =
{
	# Note the energy range of Co60 is < 1.35 MeV (1.17 MeV and 1.33 MeV primary decay)
	Energy scale factor = 1.0 # Simulate higher or lower energy by simply scaling the output energy
	# Phase space file name
	DataDir = ViewRayData/
	CoPhaseSpace = ViewRayCo.phsp
	Geometry file = vrSourceDC.py

	# Number of threads used to generate particles from source head.
	NThread = 8 # It should vary with your machine's capability. Set it as 0 to determine automatically.
	Sample Stack Depth = 400 # Increase this number if the sourch head report stack overflow. 200 will be enough.
	
	# Slightly adjust the MLC open width to match the experiment
	Leaf spacing adjustment = -0.00 # unit cm. Currently default to -0.07
	
	# Two kinds of beam files are accepted: (1)mfQAData.beams (2)PlanOverview.txt
	Beams file = $cdir$mfQAData.beams
	
	overwrite isocenter = no # Change the isocenter position using following definition in above beam file
	
	# If Beams file doesn't exists, the following customized beams will be used
	#{{
	isocenter = 0 0 0 #unit cm. Isocenter's position relative to the patient phantom center
	# Define the beams.
	# The "Beam" block can be repeated as many times as needed to define all beams.
	# In each beam block there can be an arbitrary number of segment blocks to define the segments of the beam.
	Beam =
	{
		gantry angle = 0  #in degree
		head index = 1 # The head index may be used in future version
		# The "segment" block can be repeated as many times as needed
		segment =
		{
			on time = 10  #uint second
			leaf position = 0 29  0  0
			leaf position = 13 16 -2.10 2.10 #the later one will overlap part of previous defined
		}
	}
	#}}
}


# Phantom information
PHANTOM =
{	
	# if define phantomName, the next block is ignored. It can be .phtm or .red format
	phantomName = $cdir$mfQAData.red
	
	importantRegion = $cdir$ImportantRegion.txt
	
	# if the voxel density is less than threshold, the dose will be trimmed to zero
	trim dose threshold = 0.004 # unit g/cm^3. Note air density = 0.00129
	
	# If phantomName is not defined, it'll use following parameters to define the cuboid instead.
	#{{
	DX DY DZ = 0.3 0.3 0.3 #unit cm
	
	# The later will repalce the overlapped part in the previous one.
	# The first one should be larger than any others
	# Meaning: NX, NY, NZ, x,y,z, density
	# x,y,z is the center coordinates of cuboid in isocenter coordinate system
	matrix = 101 101 101 0 0 0 1.0
	
	# If the variable SSD exist, it will redefine the position in z direction
	#SSD = 104.5 #unit cm
	
	# If the variable centerXZ exist, it will redefine the position in XZ plane
	centerXZ = 0 0 #unit cm, coordinates of the center of the cuboid relative to isocenter
	#}}
	
	
	#this is necessary if using RED file as input
	#couch density factor = 1.345 # coefficient to modify the density in couchConfig file
	# only replace couch density when couchConfig exists
	#couchConfig = couch.ced
	# if couchXYZ doesn't exist, use position definition in ViewRay's plan 
	couchXYZ = 0 4.05 0 # the top surface center's coordinates in isocenter coordinate system
	
	magnetic field strength  = 0.35 #unit Tesla
	# direction can be specified by theta, phi (e.g. 90 180) or unit vector (e.g. 0 0 -1)
	magnetic field direction = 0 0 -1
}
