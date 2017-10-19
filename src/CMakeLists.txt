cmake_minimum_required (VERSION 2.6)
project(SegCombineReconstruct)
find_package(ITK REQUIRED)
find_package(VTK REQUIRED)
include(${ITK_USE_FILE})
include(${VTK_USE_FILE})

add_executable(SegCombineReconstruct MACOSX_BUNDLE SegCombineReconstruct.cxx four-neighbour.h)
target_link_libraries(SegCombineReconstruct ${ITK_LIBRARIES})
target_link_libraries(SegCombineReconstruct ${VTK_LIBRARIES})


if(BUILD_TESTING)
  add_subdirectory(test)
endif()