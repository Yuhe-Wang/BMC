if(MSVC)
cmake_minimum_required(VERSION 3.8.2)
endif()

project(BMC) 

#set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/RunBMC")

include(cotire)

include(setupConfigurations)

include_directories("${PROJECT_SOURCE_DIR}/../Build/include")

include_directories("${PROJECT_SOURCE_DIR}/../Build/include/matio")

include_directories("${PROJECT_SOURCE_DIR}/../Build/include/wx-3.1")


if(WIN32)
	#set(CMAKE_C_FLAG "${CMAKE_C_FLAGS} /MP")
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
    ENABLE_LANGUAGE(RC) # useful in WIN32 project
    include_directories(${PROJECT_SOURCE_DIR}/../Build/lib/mswu)
	# set nvcc complile flag
    set(
        CUDA_NVCC_FLAGS
        "-Xcudafe --diag_suppress=field_without_dll_interface -Xcompiler /openmp -Xcompiler /wd4819"
    )
	set(exeExt ".exe")
elseif(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -pthread -fopenmp -Wno-unused-result -fvisibility=hidden")
    include_directories(${PROJECT_SOURCE_DIR}/../Build/lib/wx/include/gtk2-unicode-3.1)
    set(CUDA_SDK_ROOT_DIR /usr/local/cuda)
    set(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda)
    set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath=." )
    set(CMAKE_SHARED_LINKER_FLAGS  "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath=." )
	# set nvcc complile flag
    set(
        CUDA_NVCC_FLAGS
        -shared -std=c++11 -Xcompiler -std=c++11  -Xcudafe --diag_suppress=field_without_dll_interface # it's a little weird, -Xcompiler -std=c++11 tries to replace gnu++11
    )
	set(exeExt ".app")
endif()

message(STATUS CMAKE_C_FLAGS=${CMAKE_C_FLAGS})
message(STATUS CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS})
message(STATUS CMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS})
message(STATUS CMAKE_SHARED_LINKER_FLAGS=${CMAKE_SHARED_LINKER_FLAGS})
message(STATUS CMAKE_STATIC_LINKER_FLAGS=${CMAKE_STATIC_LINKER_FLAGS})

message(STATUS CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG})
message(STATUS CMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG})
message(STATUS CMAKE_EXE_LINKER_FLAGS_DEBUG=${CMAKE_EXE_LINKER_FLAGS_DEBUG})
message(STATUS CMAKE_SHARED_LINKER_FLAGS_DEBUG=${CMAKE_SHARED_LINKER_FLAGS_DEBUG})
message(STATUS CMAKE_STATIC_LINKER_FLAGS_DEBUG=${CMAKE_STATIC_LINKER_FLAGS_DEBUG})

message(STATUS CMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE})
message(STATUS CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS CMAKE_EXE_LINKER_FLAGS_RELEASE=${CMAKE_EXE_LINKER_FLAGS_RELEASE})
message(STATUS CMAKE_SHARED_LINKER_FLAGS_RELEASE=${CMAKE_SHARED_LINKER_FLAGS_RELEASE})
message(STATUS CMAKE_STATIC_LINKER_FLAGS_RELEASE=${CMAKE_STATIC_LINKER_FLAGS_RELEASE})

link_directories(${PROJECT_SOURCE_DIR}/../Build/lib)

file(GLOB_RECURSE rw_module_cmake_list . "CMakeLists.txt")

foreach(rw_module_cmake ${rw_module_cmake_list})
    string(LENGTH ${PROJECT_SOURCE_DIR} module_prefix_length)
    string(SUBSTRING ${rw_module_cmake} ${module_prefix_length} -1 rw_module_name)
    string(SUBSTRING ${rw_module_name} 1 -1 rw_module_name)
    string(FIND ${rw_module_name} "/" rw_module_name_length)

    if(NOT (${rw_module_name_length} STREQUAL -1))
        string(SUBSTRING ${rw_module_name} 0 ${rw_module_name_length} rw_module_name)
        set(rw_current_module_name ${rw_module_name})
        add_subdirectory(${rw_module_name})
        message(STATUS "The module ${rw_module_name} has been added")
    endif()
endforeach()




