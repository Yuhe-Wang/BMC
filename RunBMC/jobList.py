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

# This jobList.py is an overall configuration of this program. Please configure each job separately.

# This block works like macro definitions in C language. Use $variable$ to access them in job config files
# If you define SSD = 100 here, for example, the config item SSD = $SSD$ in your job config file is now equivalent to SSD = 100 
ConfigVariables= 
{
	dir = patients/
	#fname = gZeus
	SSD = 100 #unit cm
	LeafAdjustment = -0.07 # unit cm
}

# automatically search and execute jobs in the directory named BatchMode
# Each directory under BatchMode is treated as a standalone job and must contain only one config file
# If a file with extension of ".skip" exists in the job directory, this job will be skipped
# If BatchMode is turned on, the following job list will be ignored
BatchMode = yes


# Manually specify the job configurations. It can be repeated as many times as needed.
# Jobs will be executed by their order in this list.
# Each job has several "pre command", one "config file", several "post command" and one "rest time"
# "config file" is required, while others are optional.
# "pre command" and "post command" depends on the system: windows or linx, which will execute in sequence.
# Note the relative directory is calculated from where the executable file locates.
job=
{
	repeat times = 1 # This is meaningful only when you set to reuse last result. 0 to skip this job
	# pre command = 
	#config file = $dir$T750.py
	config file = $dir$K80.py
	# post command = move /y *.dose $dir$
	# rest time = 0.0 # Time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = 4.2cm/K80.py
	# post command = move /y *.dose 14.7cm/
	# rest time = 0.0 # Time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = $dir$config-10.5.py
	# post command = move /y *.dose 14.7cm/
	# rest time = 0.0 # Time to rest the CPU, unit minutes
}


# ------------------------------------------------------------------------------------------------- #
# Only useful when you run it on a cluster with many computation nodes.
confirm to exit = no # If yes, it will print "Press enter to exit..." and wait for user's confirmation before exiting
test node speed = no 

email= # Currently we only support gmail notification
{
	email notification = no # yes or no
	email encrypted file = myEmail.encrypt # Please use the program ConfigEmail.exe to create your own encryped file
	#If a valid "email encrypted file" exists, the following address/password will be ignored
	email address = your email
	email password = secret
}



