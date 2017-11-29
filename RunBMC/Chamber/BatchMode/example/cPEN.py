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
NBatch  = 50 # Default 50
NThread = 1 # If not defined, the maximum cores will be used
# note this file name should not contain extension. It will automatically add ".dose"
output file name = $cdir$PENELOPE
log file directory = $cdir$
log file append = no # no to delete all previous log content if it detects the same file name
# a short description stated in the log header, recommended if you choose to append log content
log description = NA

# PENELOPE information
PENELOPE =
{
	material database = C:\PENELOPE_DataBase.7z
	# if use last generated material data to boost up simulation.
	# Must turn it off during first run after modifying Emax or MAT#
	useLastMats = yes 
	# last dose file = $cdir$$fname$.dose # if it exists, we will try to reuse it
	#EMax = 1.35e6   #unit eV
	DSMax = 0.1  #unit cm

	# Sub block describing the simulation parameters of each material
	# You can repeat this block as many as necessary
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
		ID number= 29    # PENELOPE material index
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
		ID number= 228    # PENELOPE material index
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
		ID number= 126    # PENELOPE material index, C552 air equivalent
		XCSE = 512
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
		ID number= 225    # PENELOPE material index
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
		ID number= 104    # PENELOPE material index, air
		XCSE = 512
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
	CoPhaseSpace = C:\ViewRayCo.phsp
	Geometry file = vrSourceDC.py
	
	# thread num used to generate particle from source head
	NThread = 1
	Sample Stack Depth = 400
	
	# two kinds of beam files are supported: (1)mfQAData.beams (2)PlanOverview.txt
	#Beams file = $cdir$mfQAData.beams
	
	# If Beams file doesn't exists, the following customized beams will be used
	#{{
	Leaf spacing adjustment = $LeafAdjustment$ # unit cm
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
	geometry file = $cdir$A18.geom # the geometry config file of EGS++
	magnetic field = -0.35 #unit Tesla. Possitive +z direction, negative -z direction
}


