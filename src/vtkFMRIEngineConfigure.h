/* 
 * Here is where system computed values get stored.
 * These values should only change when the target compile platform changes.
 */

/* #undef VTKFMRIENGINE_BUILD_SHARED_LIBS */
#ifndef VTKFMRIENGINE_BUILD_SHARED_LIBS
#define VTKFMRIENGINE_STATIC
#endif

#if defined(WIN32) && !defined(VTKFMRIENGINE_STATIC)
#pragma warning ( disable : 4275 )

#if defined(vtkFMRIEngine_EXPORTS)
#define VTK_FMRIENGINE_EXPORT __declspec( dllexport ) 
#else
#define VTK_FMRIENGINE_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_FMRIENGINE_EXPORT
#endif
