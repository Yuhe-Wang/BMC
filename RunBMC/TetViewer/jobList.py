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

# jobs will be executed by their order in this list
# they should have the same blocks title "Job"
# each job has several "pre command", one "config file", several "post command" and one "rest time"
# "config file" is required, while others are optional.
# "pre command" and "post command" depends on the system: windows or linx, which will execute in sequence.
# note current working directory is where the executable file locates

ConfigVariables= # it works like macro definitions in C language. Use $var$ to access them in config files
{
	dir = patients/
	#fname = perfect
	SSD = 100 #unit cm
	LeafAdjustment = -0.07 # unit cm
}

BatchMode = yes #on to execute jobs in folder BatchMode

job=
{
	repeat times = 1
	# pre command = 
	#config file = $dir$T750.py
	config file = $dir$K80.py
	# post command = move /y *.dose $dir$
	# rest time = 0.0 # time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = 4.2cm/K80.py
	# post command = move /y *.dose 14.7cm/
	# rest time = 0.0 # time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = 4.2cm/90.py
	# post command = move /y *.dose 14.7cm/
	# rest time = 0.0 # time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = 4.2cm/180.py
	# post command = move /y *.dose 14.7cm/
	# rest time = 0.0 # time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = 4.2cm/270.py
	# post command = move /y *.dose 14.7cm/
	# rest time = 0.0 # time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = 27.3cm/K80.py
	# post command = move /y *.dose $dir$
	# rest time = 0.0 # time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = 10.5cm/K80.py
	# post command = move /y *.dose 14.7cm/
	# rest time = 0.0 # time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = $dir$config-6.py
	# post command = move /y *.dose $dir$
	# rest time = 0.0 # time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = $dir$config-8.py
	# post command = move /y *.dose 14.7cm/
	# rest time = 0.0 # time to rest the CPU, unit minutes
}

job=
{
	repeat times = 0
	# pre command = 
	config file = $dir$config-10.5.py
	# post command = move /y *.dose 14.7cm/
	# rest time = 0.0 # time to rest the CPU, unit minutes
}

confirm to exit = no #if yes, will print "Press enter to exit..." and wait 
test node speed = no #if your node has different hardware configuration, 

email= #currently we only support Gmail notification
{
	email notification = no # yes or no
	email encrypted file = myEmail.encrypt
	#If a valid "email encrypted file" exists, the following address/password will be ignored
	email address = yuhewang.ustc@gmail.com
	email password = secret
}



