cmake_minimum_required (VERSION 2.8)
#set project environment parameters
#set project name
set (project_name SignalModel)
#set source files list
set (SRC_LIST SignalModel.cxx)
#set include file path
set (INCLUDE_DIR  ${PROJECT_SOURCE_DIR})

project(${project_name})

find_package(VTK)
include(${VTK_USE_FILE})
message("VTK_USE_FILE is ${VTK_USE_FILE}")

find_package(ITK)
include(${ITK_USE_FILE})
message("ITK_USE_FILE is ${ITK_USE_FILE}")

add_executable(${project_name} ${SRC_LIST})

if(VTK_LIBRARIES)
	target_link_libraries(${project_name} ${VTK_LIBRARIES})
	MESSAGE("link VTK libraries done")
else()
	target_link_libraries(${project_name} $vtkHybrid vtkWidgets)
endif()

if(ITK_LIBRARIES)
	target_link_libraries(${project_name} ${ITK_LIBRARIES})
	MESSAGE("link ITK libraries done")
else()
	target_link_libraries(${project_name} $vtkHybrid vtkWidgets)
endif()
