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
NSIMU = 1e6
# Note the output file name should not contain extension. It will automatically add ".dose"
output file name = $cdir$gZeus
ViewRay's format = yes
log file directory = $cdir$
log file append = no # yes or no. Set no to delete all previous log content if it detects the same file name
# A short description stated in the log header, recommended if you choose to append log content
log description = 4.2cm square field

Simulate electron = yes # yes or no
Particle stack depth = 100 # Stack depth for each GPU thread. It may not be used in some versions
NMaxSplit = 50 # default 50. Set it to 1 to cancel splitting
Fixed Split = no # yes or no. If no, the splitting number will vary with particles' energy
EAbsPhoton = 50e3 #unit eV. Default 50e3
EAbsElectron = 50e3 #unit eV. Default 50e3
EMaxCSDA = 200e3 # unit eV. Default 200e3

# GPU thread division setting
GPU =
{
	GPU Query = no # yes or no. If yes, it will print the GPU information and pause.
	#GPU Index = 0 1   # Set multiple GPU indices to use, starting from 0.
	GPU Num = 2 # if "GPU Index" is not listed, it will try to select "GPU Num" GPU(s) to run
	GPU Block Num = 32 #128
	GPU Block Dim = 32 #256
	GPU Batch Num = 20 # Number of incident photon each thread processes
	Source Reuse Times = 1 # Default to 1. If the source generating speed lags the GPU simulation, increase it 
}
# Note the GPU will copy (GPU Batch Num)* (Block Num)* (Block Dim) particles from Source Head in order to save copy time


SOURCEHEAD =
{
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
	#Beams file = $cdir$mfQAData.beams
	
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
	#phantomName = $cdir$mfQAData.red
	
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
	# MatID = 1 #index of material listed above
	#}}
	
	
	#couchConfig = couch.ced  # It will only replace couch density when couchConfig exists
	# if couchXYZ doesn't exist, it will use position definition in ViewRay's plan 
	couchXYZ = 0 4.05 0 # The coordinates of top surface center in isocenter coordinate system
	#couch density factor = 1.345 # Density coefficient of the couch. It won't affect if couch doesn't exist in phantom
	
	magnetic field strength  = 0 #0.35 #unit Tesla
	# Direction can be specified by theta, phi or unit vector such as 0 0 1
	# and the parameter number determines the type
	magnetic field direction = 0 0 -1
}
