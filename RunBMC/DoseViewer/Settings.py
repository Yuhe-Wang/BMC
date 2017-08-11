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

# DCM_RescaleSlope = 
# DCM_RescaleIntercept = 
# DCM_BitsStored = 
# DCM_Rows = 
# DCM_Columns = 
# DCM_SliceThickness = 
DCM_KVP = 120
# DCM_PixelSpacingX = 
# DCM_PixelSpacingY =

MC engines = gZeus gMeshMM # dose engine list. They are all in dll/so form





