cmake_minimum_required(VERSION 3.0)


project(mainfunction)


include_directories("${PROJECT_SOURCE_DIR}")
#add_subdirectory(cxx)
#add the executable files
add_executable(mainfunction main.cxx)



SET(LOCAL_MODULE_SRCS
vtkMultipleInputsImageFilter.cxx
vtkActivationEstimator.cxx
vtkGLMEstimator.cxx
vtkActivationDetector.cxx
vtkGLMDetector.cxx
vtkActivationVolumeGenerator.cxx
vtkGLMVolumeGenerator.cxx
vtkCDF.cxx
vtkActivationVolumeCaster.cxx
vtkActivationRegionStats.cxx
vtkLabelMapWhitening.cxx
vtkActivationFalseDiscoveryRate.cxx
vtkIsingActivationThreshold.cxx
vtkIsingActivationTissue.cxx
vtkParzenDensityEstimation.cxx
vtkIsingConditionalDistribution.cxx
vtkIsingMeanfieldApproximation.cxx
GeneralLinearModel.cxx
SignalModel.cxx
)

configure_file(
	"${PROJECT_SOURCE_DIR}/vtkFMRIEngineConfigure.h.cin"
	"${PROJECT_SOURCE_DIR}/vtkFMRIEngineConfigure.h"
)


message("${PROJECT_SOURCE_DIR}")

add_library(fMRIEngine ${LOCAL_MODULE_SRCS})
target_link_libraries(mainfunction fMRIEngine)

find_package(ITK)
include (${ITK_USE_FILE})
message("ITK_USE_FILE IS ${ITK_USE_FILE}")

find_package(VTK)
include (${VTK_USE_FILE})
message("ITK_USE_FILE IS ${VTK_USE_FILE}")

if(VTK_LIBRARIES)
	target_link_libraries(mainfunction ${VTK_LIBRARIES})
	message("use vtk_library")
ELSE()
	target_link_libraries(mainfunction vtkHybrid vtkWidgets)
	message("vtklib not found")
endif()

if(ITK_LIBRARIES)
	target_link_libraries(mainfunction ${ITK_LIBRARIES})
ELSE()
	target_link_libraries(mainfunction vtkHybrid vtkWidgets)
endif()


