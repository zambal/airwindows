function(add_airwindows_plugin name)
    file(GLOB plug_src
            src/${name}/*.h
            src/${name}/*.cpp)
    add_library(${name} MODULE ${plug_src} ${VSTSDK_SOURCES} ${VCL_SOURCES})
    target_include_directories(${name} PRIVATE ${VSTSDK_ROOT} ${VCL_ROOT})
    set_target_properties(${name} PROPERTIES PREFIX "")
    target_link_libraries(${name} -lm)
endfunction(add_airwindows_plugin)
