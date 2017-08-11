
aux_source_directory(. DIR_SRCS)

file(GLOB headers *.h)
set(DIR_SRCS ${DIR_SRCS} ${headers})

get_filename_component(module ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_library(${module} SHARED ${DIR_SRCS})

#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../Run)

add_dependencies(${module} Tools)

target_link_libraries(${module} Tools)

install(TARGETS ${module} 
	RUNTIME DESTINATION .
	LIBRARY DESTINATION .
)

add_custom_command(TARGET ${module}
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                       $<TARGET_FILE:${module}>
                       ${CMAKE_INSTALL_PREFIX}
                   COMMENT "Copying the target ${module} to the running directory ${CMAKE_INSTALL_PREFIX}")

set_target_properties(${module} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX})
