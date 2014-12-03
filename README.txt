/****************************************/
//class hierarchy
vtkSimpleImageToImageFilter
|___vtkActivationVolumeCaster
|___vtkLabelMapWhitening
|___vtkActivationVolumeGenerator
|     |___vtkGLMVolumeGenerator
|___vtkMultipleInputImageFilter
     |___vtkActivationRegionStats
     |___vtkIsingConditionalDistribution
     |___vtkIsingActivationThreshold
     |___vtkIsingActivationTissue
     |___vtkIsingMeanfildApproximation
     |___vtkParzenDensityEstimatoion
     |___vtkActivationEstimator
            |___vtkActivationEstimator
                |___vtkGLMEstimator
vtkObject
|___vtkActivationDector
      |___vtkGLMDector
	  
SignalModeling

/****************************************/
//workflow
//setp 1:
/***************************************
Description:
	Model a signal and generate DesignMatrix
Class:
	SignalModeling
Default Properties:
	Waveform: BoxCar
	Convolution: HRF
	Derivatives: None
	Trend Model: Discrete Cosine
Input:
	vtkFloatArray*: onset, duration, contrast_vector
	// number of components(onset or duration) = number of conditions
	// onset tuples  : onset time points
	//duration tuples: durations after the onsets
	int    numVolume //number of total input images
	float  TR        //time of repetition, 2s by default
Output:
	vtkFloatArray* designMat
	vtkFloatArray* contrast(=contrast_vector)
Function:
	void SignalModeling::SetLen(numVolume)
	void SignalModeling::SetTR(TR)
	void SignalModeling::SetOnset(onset)
	void SignalModeling::SetDuration(duration)
	void GenerateDesignMatrix(SignalModeling) ////return designMat (vtkFloatArray*)
Tips:
	number of components (designMat)=number of components(contrast)
	contrast used in next section
*****************************************/

//step 2:
/***************************************
Description:
	estimate the model and compute the best fit of the design matrix
Class:
	vtkGLMDetector
	vtkGLMEstimator
Input:
	vtkFloatArray* designMat
	vtkImageData* input  //use AddInput()
Output:
	vtkImageData* estimatedImage //beta map
Function:
	vtkGLMDetector::SetDesignMatrix(vtkFloatArray* designMat)
	vtkGLMEstimator::SetDetector(vtkGLMDetector* )
	vtkGLMEstimator::AddInput(vtkImageData* )//add all fMRImages
	vtkGLMEstimator::update()
	vtkImageData* vtkGLMEstimator::GetOutput() //get output estimated statistic map
Tips:
	number of fMRImages= numVolume
*****************************************/

//step 3;
/***************************************
Description:
	compute an activation volume(t-map)
Class:
	vtkGLMVolumeGenerator
Default Properties:
	No preWhitening
Input:
	vtkIntArray* contrast		//get from step 1
	vtkFloatArray* designMat	//get from step 1
	vtkImageData* estimatedImage // get from setp 2
Output:
	vtkImageData* GeneratedImage
Function:
	void vtkGLMVolumeGenerator::SetContrastVector(vtkIntArray* contrast)
	void vtkGLMVolumeGenerator::SetDesignMatrix(vtkFloatArray* designMat)
	void vtkGLMVolumeGenerator::Update()  //use this function for calculation
	vtkImageData* vtkGLMVolumeGenerator::GetOutput();// get output image
Tips:
	contrast in the first step added here
*****************************************/

//setp 4
/***************************************
Description:
	thresholding. final output of the function
	Determine activation status of voxels
Class:
	IsingActivationThreshold
Input:
	float p     //0.001 for high dependency
	int   dof   //numVolume-1 by default
	vtkImageData *GeneratedImage  //t-map in step 3
Output:
	vtkImageData* thresholded
Function:
	float  CDF->p2t(float p,int dof);      //generate t threshold
	void IsingActivationThreshold::Setthreshold(t); //set t threshold
	void IsingActivationThreshold::Update();  //use this function for calculate
	vtkImageData* IsingActivationThreshold::GetOutput(); //use this function to get output image
Tips:
	output properties:
		posactive voxel : 300
		noneactive voxel: 0
		negactive voxel : 150
		(change in vtkIsingActivationThreshold.cxx #--38)
*****************************************/

//timecourse
/**********************************
Description:
	get voxel/region timecourse
Class:
	vtkGLMEstimator	
Input: 
	float coordinate[3];   //coordinate of the certain voxel
	vtkFloatArray* region; //an array that holds all coordinates of the region
	vtkImageData*  input;
Output:
	vtkFloatArray*  //return value of the function
Function:
	void vtkGLMEstimator::AddInput(vtkImageData*  input);  //add all volumes 
	vtkFloatArray* vtkGLMEstimator::GetTimecourse(float* coordinate)
	vtkFloatArray* vtkGLMEstimator::GetRegionTimecourse(vtkFloatArray* region)
Tips:
***********************************/



/****************************************/
//main class details


//Class SignalModeling
/***************************************
Class SignalModeling
public:
	SignalModeling();
	~SignalModeling();
	//Description: 
	//model convolove with the FRF func
	vtkFloatArray* ConvolveWithHRF();
	//Description:
	// takes a list in, subsamples it to a new length and returns the new list.
	vtkFloatArray* GaussianDownsampleList(vtkFloatArray* inputList);
	//Description:
	//Generate DCBasis
	vtkFloatArray* BuildDiscreteCosineBasis();
	//Description:
	//compute baseline
	vtkFloatArray* Baseline();
	//Description:
	//Parameters input
	void SetLen(int len) {this->numTimePoints=len;}
	void SetTR(float TR) {this->TR=TR;}
	void SetOnset(vtkFloatArray* onset) {this->onset=onset;}
	void SetDuration(vtkFloatArray* duration) {this->duration=duration;}
private:
//parameters:
	float TR;      //2s by default
	int numTimePoints;   //number of input volumes
	float TimeIncrement; //0.1 by default
	vtkFloatArray* flipHRF;
	vtkFloatArray* onset;
	vtkFloatArray* duration;
	int NumberOfCosineBasis;  //get from FindNumberOfCosineBasis
// functions:
	//Description:
	//stimuli signal generator. wave-form: BoxCar
	//sampled by 0.1s, use this array to convolve with the HRF
	vtkFloatArray* ComputeBoxCar(vtkFloatArray* onset,vtkFloatArray* duration);
	//Description:
	//HRF point data. Generate by 2 gamma functions and sampled by 0.1s and flip the list to convolove
	vtkFloatArray* ComputeHRF();
	//Description：
	//To downsample a list, setup a gaussian kernal for convolution. 
	//invoked by GaussianDownsampleList()
	vtkFloatArray* ComputeGaussianFilter();
	//Description:
	//func as name invoked by FindNumberOfCosineBasis
	float ComputeDefaultHighpassTemporalCutoff();
	int FindNumberOfCosineBasis();
	//Description:
	//rescale the array to range [-1.0,1.0]
	vtkFloatArray* Rescale(vtkFloatArray* Inputlist);
//external function
//Description:
//put the parameters together and generate a designmatrix
//Data structure:
//component[1..n]: condition[1...n];component[n+1]:baseline;component[n+2...n+2+NumberOfCosineBasis]:DCBasis[1...NumberOfCosineBasis]
vtkFloatArray* GenerateDesignMatrix(SignalModeling* siglm);
***************************************/

//Class: vtkGLMDetector
/***************************************
Class: vtkGLMDetector
//is used to compute voxel activation based on
//paradigm and detection method (GLM or MI).
//Default is GLM in this project
Main functions:
    // Description:
    // Gets the design matrix 
    vtkFloatArray *GetDesignMatrix();
    // Description:
    // Sets the design matrix 
    void SetDesignMatrix(vtkFloatArray *designMat);
    void SetAR1DesignMatrix ( );
    vtkFloatArray *GetAR1DesignMatrix ( );
    vtkFloatArray *GetResiduals ( );
    // Description:
    // Fits linear model (voxel by voxel) 
    void FitModel(vtkFloatArray *timeCourse, float *beta, float *chisq ); 
    // This uses the first estimated beta to subtract
    // the model from the data to compute errors like so:
    // Y = XB + e --> e = Y-XB_hat
    void ComputeResiduals ( vtkFloatArray *timeCourse, float *beta );
    // Description::
    // This uses the residuals to compute the correlation coefficient
    // at lag 1 used in pre-whitening for data and residuals.
    float ComputeCorrelationCoefficient ( );
    // Description:
    // This whitens the DesignMatrix and timeCourse.
    // Saves the whitened design matrix in AR1DesignMatrix,
    // and replaces the extracted timecourse by new values.
    void PreWhitenDataAndResiduals (vtkFloatArray *timeCourse, float corrCoeff);
    // Description:
    // Sets the AR1DesignMatrix
    // and turns on the whitening flag.
    void EnableAR1Modeling ( );
    // Description:
    // This sets the AR1DesignMatrix to NULL,
    // sets the WhiteningMatrix to NULL,
    // and turns off the whitening flag.
    void DisableAR1Modeling ( );
***************************************/

//Superclass vtkMultipleInputImageFilter
/***************************************
Superclass vtkMultipleInputImageFilter
Description:
//Multiple inputs and single output filter
Function:
	void AddInput(vtkImageData*)
	vtkImageData* GetInput(int)
	//Description:
	//this is a virtual function and subclass need to implement this function
	//while in the vtk5.10 this function does not func well in multiple inputs but cannot delete
	//so another function replace it and SimpleExecute() does nothing.
	virtual SimpleExecute(vtkImageData *inputs, vtkImageData* output)
***************************************/

//Class vtkGLMEstimator
/***************************************
Class vtkGLMEstimator
Superclass vtkMultipleInputImageFilter
Description:
	take the following functions instead of SimpleExecute() to enable multiple inputs
	int RequestData(vtkInformation *vtkNotUsed(request),vtkInformationVector **inputVector,vtkInformationVector *outputVector)
	int FillInputPortInformation( int port, vtkInformation* info )
and 
	Getinputx(int i,vtkInformationVector **inputVector) instead of GetInput(int)
***************************************/
	
	
	
	
	
	