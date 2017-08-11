if(NOT SET_UP_CONFIGURATIONS_DONE)
    set(SET_UP_CONFIGURATIONS_DONE 1)

    # No reason to set CMAKE_CONFIGURATION_TYPES if it's not a multiconfig generator
    # Also no reason mess with CMAKE_BUILD_TYPE if it's a multiconfig generator.
    if(CMAKE_CONFIGURATION_TYPES) # multiconfig generator?
        set(CMAKE_CONFIGURATION_TYPES "MyDebug;Release;" CACHE STRING "" FORCE) 
    else()
        if(NOT CMAKE_BUILD_TYPE)
            message("Defaulting to release build.")
            set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
        endif()
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY HELPSTRING "Choose the type of build")
        # set the valid options for cmake-gui drop-down list
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "MyDebug;Release")
    endif()
    # now set up the Profile configuration
	if(WIN32) # compiler flag specific to cl.exe
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP /MD")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /MD")
		set(CMAKE_C_FLAGS_MYDEBUG "/Zi /Ob0 /Od /RTC1")
		set(CMAKE_CXX_FLAGS_MYDEBUG "/Zi /Ob0 /Od /RTC1")
		set(CMAKE_EXE_LINKER_FLAGS_MYDEBUG "/debug /INCREMENTAL")
		set(CMAKE_SHARED_LINKER_FLAGS_MYDEBUG "/debug /INCREMENTAL")
	endif()
endif()