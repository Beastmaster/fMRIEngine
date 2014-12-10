//head files for fMRIEngine
#include <vtkFMRIEngineConfigure.h>
#include "vtkMultipleInputsImageFilter.h"
#include "vtkImageData.h"
#include "FMRIEngineConstants.h"              
#include "vtkGLMDetector.h"
#include "GeneralLinearModel.h"
#include "vtkGLMEstimator.h"
#include "vtkActivationDetector.h"
#include "vtkGLMVolumeGenerator.h"
#include "vtkActivationEstimator.h"
#include "vtkIsingActivationThreshold.h"
#include "vtkActivationFalseDiscoveryRate.h"
#include "vtkIsingActivationTissue.h"
#include "vtkActivationRegionStats.h"
#include "vtkIsingConditionalDistribution.h"
#include "vtkActivationVolumeCaster.h"//zero out values out of range of low/high threshold
#include "vtkIsingMeanfieldApproximation.h"
#include "vtkActivationVolumeGenerator.h"
#include "vtkLabelMapWhitening.h"
#include "vtkCDF.h"            
#include "vtkMultipleInputsImageFilter.h"
#include "vtkFMRIEngineConfigure.h"
#include "vtkParzenDensityEstimation.h"
#include "SignalModel.h"
#include <iostream>

//Descirption:headfile used by ReadFile()
//Info : by qinshuo
#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkDICOMImageReader.h>
#include <vtkBMPWriter.h>
//Descirption:headfile used by ReadFile()
//Info : by qinshuo
#include <vtkImageCast.h>

//headfiles for write designmatrix as test
#include <fstream>

//headfiles for filename generator
#include "itkNumericSeriesFileNames.h"

//thresholding
#include "vtkImageThreshold.h"
//gaussian smoothing
#include "vtkImageGaussianSmooth.h"

//Description:thresholding
//Info:by qinshuo --- 2014.12.9
vtkImageData* Pre_Threshold(vtkImageData*,float lower,float upper);

//Description:gaussian smoothing
//Info:by qinshuo --- 2014.12.9
vtkImageData* Smooth(vtkImageData*);


//Description:generate filenames
//Info:by qinshuo
char** FilenameGenerator(int first,int last,char* prefix);


//Description:read/write file
//Info:by qinshuo 
//void fileout(vtkFloatArray* input,char* name);
vtkImageData * ReadmhaFile(char* filename);
vtkImageData * ReadvtiFile(char* filename);
vtkImageData * ReaddcmFile(char* filename);
void WritemhaFile(vtkImageData* ImageToWrite,char* filename);
void WritevtiFile(vtkImageData* ImageToWrite,char* filename);
void WritebmpFile(vtkImageData* ImageToWrite,char* filename);
//Description: Set initial designMat
//Info: by qinshuo 9.5
vtkFloatArray * SetdesignMat();
vtkImageData* SetGLMEstimator(vtkFloatArray* designMat,char * prefix);
vtkImageData* SetGLMVolumeGenerator(vtkFloatArray* designMat,vtkMetaImageReader* reader);
vtkImageData* ActivationThreshold(vtkImageData * image);

//Description: Set initial designMat
//Info: by qinshuo 9.5
void write_txt(vtkImageData*,char*);

int main (int argc,char ** argv)
{
	/****filenames input***/
	char **filename;//filename for reader
	//char* prefix="dicom1/IMG-0002-";//filename prefix
	char* prefix="lft-mha/20130610_144057LTHANDMOTORs701a1007_";
	int start = 1;
	//char* thresholded_name = "threshold_byme.bmp";
	char* thresholded_name = "threshold_byme.mha";
/***************************************
Description:
	Model a signal and generate DesignMatrix
Input:
	vtkFloatArray*: onset, duration, contrast_vector
	int    numVolume //number of input images
	float  TR        //time of repetition
Output:
	vtkFloatArray* designMat
	vtkFloatArray* contrast(=contrast_vector)
Function:
	SignalModeling::SetLen(numVolume)
	SignalModeling::SetTR(TR)
	SignalModeling::SetOnset(onset)
	SignalModeling::SetDuration(duration)
	GenerateDesignMatrix(SignalModeling) ////return designMat
Tips:
	number of tuple(designMat)=number of component(contrast)
	contrast used in next section
*****************************************/
	//paradigm parameters
	int numVolume=60;
	float TR=3;
	vtkFloatArray* onset=vtkFloatArray::New();
	vtkFloatArray* duration=vtkFloatArray::New();
	vtkIntArray* contrast=vtkIntArray::New();

	onset->SetNumberOfComponents(1);
	duration->SetNumberOfComponents(1);
	onset->InsertComponent(0,0,10);
	//onset->InsertComponent(1,0,40);
	//onset->InsertComponent(2,0,70);
	onset->InsertComponent(1,0,30);
	onset->InsertComponent(2,0,50);
	//onset->InsertComponent(2,1,80);
	duration->InsertComponent(0,0,10);
	//duration->InsertComponent(1,0,10);
	//duration->InsertComponent(2,0,10);
	duration->InsertComponent(1,0,10);
	duration->InsertComponent(2,0,10);
	//duration->InsertComponent(1,1,10);
	//duration->InsertComponent(2,1,10);

	int xx_t = onset->GetNumberOfTuples();
	int xx_c = onset->GetNumberOfComponents();
	//parameter set up
	SignalModeling* siglm=new SignalModeling;
	siglm->SetLen(numVolume);
	siglm->SetTR(TR);
	siglm->SetOnset(onset);
	siglm->SetDuration(duration);

	//generate design matrix
	vtkFloatArray* designMat=GenerateDesignMatrix(siglm);

	std::ofstream file ("designmatrix.txt");
	
	int tuples=designMat->GetNumberOfTuples();
	int comp  =designMat->GetNumberOfComponents();
	
	if(file.is_open())
	{
		file<<"point data"<<std::endl;
		file<<"number of tuples:"<<tuples<<endl;
		file<<"number of component:"<<comp<<endl;

		for(int j=0;j<tuples;j++)
		{
			//file<<"\ncomponent:"<<j<<"\n";
			for(int i=0;i<comp;i++)
		   {
			 file<<designMat->GetComponent(j,i)<<" ";
			}
			file<<";";
		}
	}
	file.close();


	contrast->SetNumberOfComponents(1);
	for(int i=0;i<designMat->GetNumberOfComponents();i++)
	{
		contrast->InsertComponent(i,0,0);
	}
	contrast->InsertComponent(0,0,1);
	contrast->InsertComponent(1,0,0);


/***************************************
Description:
	estimate the model and compute the best fit of the design matrix
Input:
	vtkFloatArray* designMat
	vtkImageData* 
Output:
	vtkImageData* estimatedImage
Function:
	vtkGLMDetector::SetDesignMatrix(vtkFloatArray* designMat)
	vtkGLMEstimator::SetDetector(vtkGLMDetector* )
	vtkGLMEstimator::AddInput(vtkImageData* )//add all fMRImages
	vtkGLMEstimator::update()
	vtkImageData* vtkGLMEstimator::GetOutput();
Tips:
	number of fMRImages= numVolume
*****************************************/


	//class:vtkGLMDector
	//superclass: vtkActivationDector
	vtkGLMDetector *GLMDetector=vtkGLMDetector::New();
	GLMDetector->SetDesignMatrix(designMat);
	// Description:
    // Gets/Sets the activation detection method (GLM = 1; MI = 2).
	GLMDetector->SetDetectionMethod(ACTIVATION_DETECTION_METHOD_GLM);//default parameter

	/****setup parameters for estimator*****/
	//class:vtkGLMEstimator
	///superclass:vtkActivationEstimator
	////superclass:vtkMultipleInputImageFilter
	//Description: the estimation must occur before any activation volumes can be generated
	vtkGLMEstimator* GLMEstimator=vtkGLMEstimator::New();
	//set detector first
	GLMEstimator->SetDetector(GLMDetector);
	
		//write_txt(ReadvtiFile("xxx.vti"),"xxx.txt");
		//write_txt(ReadvtiFile("estimated.vti"),"estimated.txt");
	/*get input images for test*/
	//char **filename;//filename for reader
	//char* prefix="dicom1/IMG-0002-000";//filename prefix
	//int start = 1;
	//generate filename series
	filename=FilenameGenerator(start,start+numVolume-1,prefix );	
	//add input to GLMEstimator
	/*ReadmhaFile(filename[i])for test*/
	int i;
	for(i=0;i<numVolume;i++)
	{
		GLMEstimator->AddInput(Pre_Threshold(Smooth(ReadmhaFile(filename[i])),150,1000));
		//GLMEstimator->AddInput(Pre_Threshold(ReadmhaFile(filename[i]),150,1000));
	}

	//Description:
	//Whether prewhiten or not(1/0)
	GLMEstimator->SetPreWhitening(0);//default
    // Description:
    // Sets the lower threshold.
	//GLMEstimator->SetLowerThreshold(0);
    // Description:
    // Sets the cutoff frequency.
    //GLMEstimator->SetCutoff(0);
    // Description:
    // Enables or disables high-pass filtering.0:disable
    GLMEstimator->EnableHighPassFiltering(0);//default
    // Description:
    // Sets/Gets global effect.1:grand mean;2:global mean,3:pre-whiten data
	GLMEstimator->SetGlobalEffect(1);//default parameter
	
	GLMEstimator->Update();
	//return pointer
	vtkImageData *estimatedImage=vtkImageData::New();
	estimatedImage=GLMEstimator->GetOutput();
	WritevtiFile(estimatedImage,"estimated.vti");
	write_txt(estimatedImage,"estimatedx.txt");
	int imageDim[3];
	estimatedImage->GetDimensions(imageDim);
	std::cout<<"estimatedImage"<<std::endl;
	std::cout<<imageDim[0]<<imageDim[1]<<imageDim[2]<<std::endl;

/***************************************
Description:
	compute an activation volume(t-map)
Input:
	vtkImageData* estimatedImage //
Output:
	vtkImageData* GeneratedImage
Function:
	vtkGLMVolumeGenerator::SetContrastVector(vtkIntArray* contrast)
	vtkGLMVolumeGenerator::SetDesignMatrix(vtkFloatArray* designMat)
	vtkGLMVolumeGenerator::GetOutput();
Tips:
	contrast in the first step added here
*****************************************/
	vtkGLMVolumeGenerator *GLMVolmeGenerator=vtkGLMVolumeGenerator::New();
	//Descirption:
	//Determine whether to prewhiten 0/1
	GLMVolmeGenerator->SetPreWhitening(0);//default
    // Description:
    // Sets the contrast vector. 
    GLMVolmeGenerator->SetContrastVector(contrast);
    // Description:
    // Sets the design matrix 
    GLMVolmeGenerator->SetDesignMatrix(designMat);
	//input estimatedimage
	vtkImageData* estimatedImage1 = vtkImageData::New();
	estimatedImage1 = ReadvtiFile("xxx.vti");
	/////////////
	GLMVolmeGenerator->SetInput(estimatedImage);
	GLMVolmeGenerator->Update();
	vtkImageData *GeneratedImage=GLMVolmeGenerator->GetOutput();
/***************************************
Description:
	thresholding. final output of the function
Input:
	vtkImageData *GeneratedImage
Output:
	vtkImageData* thresholded
Function:
	vtkImageData* IsingActivationThreshold::GetOutput()
Tips:
	posactive voxel : 300
	noneactive voxel: 0
	negactive voxel : 150
	change in vtkIsingActivationThreshold.cxx #--38
*****************************************/
	//threshold is input as t
	//input: p,dof;
	//output:t
	double t;
	double p=0.001;
	long dof=59;//default=numVolume-1

	vtkCDF * CDF=vtkCDF::New();  
	t=CDF->p2t(p,dof);
	std::cout<<"value of t: "<<t<<std::endl;

	vtkIsingActivationThreshold *IsingActivationThreshold=vtkIsingActivationThreshold::New();
	//t-map added here
	IsingActivationThreshold->AddInput(GeneratedImage);
	//set threshold
	IsingActivationThreshold->Setthreshold(t);
	IsingActivationThreshold->SetthresholdID(1);//default
	IsingActivationThreshold->Update();
	//output image
	vtkImageData* thresholded= IsingActivationThreshold->GetOutput();
	//write an output for test
	vtkImageCast* castFilter = vtkImageCast::New();
	castFilter->SetInput(thresholded);
	castFilter->SetOutputScalarTypeToUnsignedChar();
	castFilter->Update();
	//WritebmpFile(castFilter->GetOutput(),thresholded_name);
	WritemhaFile(thresholded,thresholded_name);


	system("pause");
	return EXIT_SUCCESS;
}


//Description: visualize image
//Info: by qinshuo

vtkImageData * ReadmhaFile(char* filename)
{
	vtkMetaImageReader *reader=vtkMetaImageReader::New();
	reader->SetFileName(filename);
	reader->Update();
	std::cout<<"reader "<<filename<<" successfully"<<std::endl;
	return reader->GetOutput();
}
vtkImageData * ReadvtiFile(char* filename)
{
	vtkXMLImageDataReader *reader=vtkXMLImageDataReader::New();
	reader->SetFileName(filename);
	reader->Update();
	std::cout<<"reader "<<filename<<" successfully"<<std::endl;
	return reader->GetOutput();
}
vtkImageData * ReaddcmFile(char* filename)
{
	vtkDICOMImageReader *reader=vtkDICOMImageReader::New();
	reader->SetFileName(filename);
	reader->Update();
	std::cout<<"reader "<<filename<<" successfully"<<std::endl;
	return reader->GetOutput();
}

void WritemhaFile(vtkImageData* ImageToWrite,char* filename)
{
	vtkMetaImageWriter *writer2 =vtkMetaImageWriter::New();
	  writer2->SetFileName(filename);
	 writer2->SetInput(ImageToWrite);
	  writer2->Write();
	std::cout<<"write "<<filename<<" successfully"<<std::endl;
}
void WritevtiFile(vtkImageData* ImageToWrite,char* filename)
{
	vtkXMLImageDataWriter *writer2 =vtkXMLImageDataWriter::New();
	  writer2->SetFileName(filename);
	 writer2->SetInput(ImageToWrite);
	  writer2->Write();
	std::cout<<"write "<<filename<<" successfully"<<std::endl;
}
void WritebmpFile(vtkImageData* ImageToWrite,char* filename)
{
	vtkBMPWriter * writer = vtkBMPWriter::New();
	writer->SetFileName(filename);
	writer->SetInput(ImageToWrite);
	writer->Write();
	std::cout<<"write "<<filename<<" successfully"<<std::endl;
}


vtkFloatArray* SetdesignMat()
{
	int numVolume=90;
	vtkFloatArray* onset=vtkFloatArray::New();
	onset->SetNumberOfComponents(2);
	//condition 1
	onset->InsertComponent(0,0,10);
	onset->InsertComponent(1,0,40);
	onset->InsertComponent(2,0,70);
	//condition 2
	onset->InsertComponent(0,1,20);
	onset->InsertComponent(1,1,50);
	onset->InsertComponent(2,1,80);
	std::cout<<onset->GetNumberOfTuples()<<std::endl;

	vtkFloatArray* duration=vtkFloatArray::New();
	duration->SetNumberOfComponents(2);
	//condition 1
	duration->InsertComponent(0,0,10);
	duration->InsertComponent(1,0,10);
	duration->InsertComponent(2,0,10);
	//condition 2
	duration->InsertComponent(0,1,10);
	duration->InsertComponent(1,1,10);
	duration->InsertComponent(2,1,10);
	std::cout<<duration->GetNumberOfTuples()<<std::endl;
	
	SignalModeling* siglm=new SignalModeling;
	//parameter set up
	siglm->SetLen(numVolume);
	siglm->SetTR(2);
	siglm->SetOnset(onset);
	siglm->SetDuration(duration);
	
	vtkFloatArray* init=vtkFloatArray::New();
	//compute work flow

	init=GenerateDesignMatrix(siglm);

	system("pause");
	return init;
}


char** FilenameGenerator(int first,int last,char* prefix )
{
	char** returnName;
	//char* suffix="%05d.dcm";
	char* suffix="%02d.mha";
	char* name2=new char[strlen(suffix)+strlen(prefix)];
	strcpy(name2,prefix);
	strcat(name2,suffix);
  typedef itk::NumericSeriesFileNames    NameGeneratorType;
  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
  nameGenerator->SetSeriesFormat(name2);//( "BinaryImage%03d.png" );
  nameGenerator->SetStartIndex( first );
  nameGenerator->SetEndIndex( last );
  nameGenerator->SetIncrementIndex( 1 );

     std::vector< std::string > fileNames = nameGenerator->GetFileNames();
	
	 returnName=new char* [fileNames.size()];
  for(unsigned int i = 0; i < fileNames.size(); ++i)
    {
    returnName[i]=new char[fileNames[i].size()];
	strcpy(returnName[i],fileNames[i].c_str());
	std::cout<<returnName[i]<<std::endl;
    }
  returnName[fileNames.size()]=NULL;

  return returnName;
}

void write_txt(vtkImageData* input,char* filename)
{
    std::cout<<"writing "<<filename<<std::endl;
	vtkDataArray *scalarsInOutput = input->GetPointData()->GetScalars(); 
	ofstream file (filename);
	if(file.is_open())
	{	
		int vox=0;
		int imgDim[3];
		int noOfRegressors=4;
		input->GetDimensions(imgDim);

		for (int kk = 0; kk < imgDim[2]; kk++)
		{
			for (int jj = 0; jj < imgDim[1]; jj++)
			{
				for (int ii = 0; ii < imgDim[0]; ii++)
				{
					file<<"\n voxel: ii="<<ii<<",jj="<<jj<<",kk="<<kk<<endl;
					// put values into output volume
					 int yy = 0;
					 // betas
					 file<<"beta:\t";
					 for (int dd = 0; dd < noOfRegressors; dd++)
					 {
						 file<<scalarsInOutput->GetComponent(vox, yy++);
						 file<<"\t";
					 }
					  // chisq and p
					 file<<"chisq:"<<scalarsInOutput->GetComponent(vox, yy++)<<"\t";
					 file<<"p value:"<<scalarsInOutput->GetComponent(vox, yy++)<<"\t";
					 // % signal changes
					  file<<"pSigChanges:";
					 for (int dd = 0; dd < noOfRegressors; dd++)
					 {
					  file<<scalarsInOutput->GetComponent(vox, yy++)<<"\t";
					 }
					 vox++;
				 }
			 }
		  }
	 }
	file.close();
	std::cout<<"writing done: "<<filename<<std::endl;
}


vtkImageData* Pre_Threshold(vtkImageData* image_in,float lower,float upper)
{
	double out_value = 0;
	
	vtkImageThreshold* threshold_hd = vtkImageThreshold::New();
	
	//set input
	threshold_hd->SetInput(image_in);

	//set threshold range
	threshold_hd->ThresholdBetween(lower,upper);
	//replace the pixel out of range with out value
	threshold_hd->ReplaceOutOn();
	threshold_hd->SetOutValue(out_value);
	threshold_hd->Update();
	//get out thresholded image pointer
	return threshold_hd->GetOutput();
}

vtkImageData* Smooth(vtkImageData* image_in)
{
	double deviations[3];
	double radiusfactors[3];
	vtkImageGaussianSmooth* gaussianSmooth = 
		vtkImageGaussianSmooth::New();
	gaussianSmooth->SetInput(image_in);
	//set gassian parameters
	//gaussianSmooth->SetStandardDeviations(deviations);
	//gaussianSmooth->SetRadiusFactors(radiusfactors)
	gaussianSmooth->Update();
	return gaussianSmooth->GetOutput();
}