
aux_source_directory(. DIR_SRCS)

file(GLOB headers *.h)
set(DIR_SRCS ${DIR_SRCS} ${headers})

get_filename_component(module ${CMAKE_CURRENT_SOURCE_DIR} NAME)

#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../Run/${module})

if(WIN32)
    add_executable(${module} WIN32 ${DIR_SRCS} Resource.rc)
    add_definitions(-D_FILE_OFFSET_BITS=64 -DWXUSINGDLL)
	target_link_libraries(${module} opengl32.lib glu32.lib wxbase31u.lib wxbase31u_xml.lib wxmsw31u_adv.lib wxmsw31u_core.lib wxmsw31u_gl.lib wxmsw31u_html.lib wxmsw31u_xrc.lib)
elseif(UNIX)
    add_executable(${module} ${DIR_SRCS})
    add_definitions(-D_FILE_OFFSET_BITS=64 -DWXUSINGDLL -D__WXGTK__)
    include_directories("${PROJECT_SOURCE_DIR}/../Build/lib/wx/include/gtk3-unicode-3.1")
    target_link_libraries(${module} X11 GL GLU stdc++fs wx_gtk3u_xrc-3.1 wx_gtk3u_html-3.1 wx_gtk3u_gl-3.1 wx_gtk3u_adv-3.1 wx_gtk3u_core-3.1 wx_baseu_xml-3.1 wx_baseu-3.1)
endif()

add_dependencies(${module} Tools)

target_link_libraries(${module} Tools)

add_dependencies(${module} DcmTK)

target_link_libraries(${module} DcmTK)

add_dependencies(${module} GPUGamma)

target_link_libraries(${module} GPUGamma)

target_link_libraries(${module} mgl)

target_link_libraries(${module} dcmdata dcmimgle dcmrt ofstd oflog)

set_target_properties(${module} PROPERTIES SUFFIX ${exeExt})

set_target_properties(${module} PROPERTIES COTIRE_CXX_PREFIX_HEADER_INIT "PreCompile.h")
cotire(${module})

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
