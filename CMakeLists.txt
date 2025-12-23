cmake_minimum_required(VERSION 3.20)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(PracticalWork_Semester3 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED COMPONENTS CXX)

set(VTK_GROUP_ENABLE_VR OFF CACHE BOOL "" FORCE)
set(VTK_MODULE_ENABLE_VTK_RenderingOpenVR NO CACHE STRING "" FORCE)
set(VTK_MODULE_ENABLE_VTK_RenderingVR NO CACHE STRING "" FORCE)

find_package(VTK REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  FiltersCore
  FiltersGeometry
  RenderingCore
  RenderingOpenGL2
  InteractionStyle
  IOXML
)

add_executable(my_app
 src/main.cpp
 src/grid.cpp
 src/math_utils.cpp
 src/visualize_split.cpp
)

target_link_libraries(my_app PRIVATE
  VTK::CommonCore
  VTK::CommonDataModel
  VTK::FiltersCore
  VTK::FiltersGeometry
  VTK::RenderingCore
  VTK::IOXML
  VTK::RenderingOpenGL2
  VTK::InteractionStyle
)

vtk_module_autoinit(
  TARGETS my_app
  MODULES ${VTK_LIBRARIES}
)
