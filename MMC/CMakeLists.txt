
aux_source_directory(. DIR_SRCS)

file(GLOB headers *.h)
set(DIR_SRCS ${DIR_SRCS} ${headers})

get_filename_component(module ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_library(${module} SHARED ${DIR_SRCS})

#set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../Run/${module})

set_target_properties(${module} PROPERTIES PREFIX "")

add_dependencies(${module} SourceHead)

target_link_libraries(${module} SourceHead)

if(WIN32)
    
elseif(UNIX)
    target_link_libraries(${module} gomp)
endif()

add_dependencies(${module} MCShell)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/../MCShell/MCShell${exeExt}" DESTINATION . RENAME ${module}${exeExt})

install(TARGETS ${module} 
	RUNTIME DESTINATION .
	LIBRARY DESTINATION .
)

add_custom_command(TARGET ${module}
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                       $<TARGET_FILE:${module}>
                       ${CMAKE_INSTALL_PREFIX}
				    COMMAND ${CMAKE_COMMAND} -E copy
				   "${CMAKE_CURRENT_BINARY_DIR}/../MCShell/MCShell${exeExt}"
				   ${CMAKE_INSTALL_PREFIX}/${module}${exeExt}
                   COMMENT "Copying the target ${module} to the running directory ${CMAKE_INSTALL_PREFIX}")

set_target_properties(${module} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX})
