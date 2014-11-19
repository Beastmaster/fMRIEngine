/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkActivationVolumeCaster.h,v $
  Date:      $Date: 2006/01/06 17:57:36 $
  Version:   $Revision: 1.4 $

=========================================================================auto=*/

#ifndef __vtkActivationVolumeCaster_h
#define __vtkActivationVolumeCaster_h


#include <vtkFMRIEngineConfigure.h>
#include "vtkSimpleImageToImageFilter.h" 

class  VTK_FMRIENGINE_EXPORT vtkActivationVolumeCaster : public vtkSimpleImageToImageFilter
{
public:
    static vtkActivationVolumeCaster *New();
    vtkTypeMacro(vtkActivationVolumeCaster, vtkSimpleImageToImageFilter);

    // Description:
    // Gets the low/high range 
    vtkGetMacro(LowRange, short);// short GetLowRange() {return LowRange;}
    vtkGetMacro(HighRange, short);//void  GetHighRange(){return HighRange;}

    // Description:
    // Sets the lower/upper threshold  
    vtkSetMacro(LowerThreshold, float);//void SetLowerThreshold(float f) {LowerThreshold=f;}
    vtkSetMacro(UpperThreshold, float);//void SetUpperThreshold(float f) {UpperThreshold=f;}

protected:
    vtkActivationVolumeCaster();

    float LowerThreshold;
    float UpperThreshold;
    short LowRange;
    short HighRange;

    void SimpleExecute(vtkImageData *input,vtkImageData *output);
};


#endif
