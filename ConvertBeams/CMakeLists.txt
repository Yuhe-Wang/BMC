
aux_source_directory(. DIR_SRCS)

get_filename_component(module ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_executable(${module} ${DIR_SRCS})

#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../Run)

add_dependencies(${module} Tools)

target_link_libraries(${module} Tools)

add_dependencies(${module} SourceHead)

target_link_libraries(${module} SourceHead)

set_target_properties(${module} PROPERTIES SUFFIX ${exeExt})

if(WIN32)

elseif(UNIX)
    target_link_libraries(${module} -lstdc++fs)
endif()

install(TARGETS ${module} 
	RUNTIME DESTINATION .
)

add_custom_command(TARGET ${module}
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                       $<TARGET_FILE:${module}>
                       ${CMAKE_INSTALL_PREFIX}
                   COMMENT "Copying the target ${module} to the running directory ${CMAKE_INSTALL_PREFIX}")

set_target_properties(${module} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX})
