cmake_minimum_required(VERSION 3.16)
project(ecalSim)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(ecalSim ecalSim.cc ${sources} ${headers})
target_include_directories(ecalSim PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ecalSim PRIVATE ${Geant4_LIBRARIES})

foreach(_script
    macros/run.mac
    macros/init_vis.mac
    macros/vis.mac
)
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
  )
endforeach()
