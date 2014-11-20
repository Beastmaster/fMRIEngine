/************************************
Description : generate designmatrix
Info:		: by qinshuo, 2014,10,15
**************************************/
//
#include "SignalModel.h"


SignalModeling::SignalModeling()
{
	this->TR=2;
	this->numTimePoints=30;
	this->TimeIncrement=0.1;
	this->flipHRF=this->ComputeHRF();
	this->onset=NULL;
	this->duration=NULL;
}

//convolve stimuli with the HRF 
vtkFloatArray* SignalModeling::ConvolveWithHRF()
{
	//#--- Convolve image with canonical hemodynamic response
    //#--- HRF swings between -1.0 and 1.0
    //#--- want to convolve this with the signal
    //#--- which also swings from -1.0 to 1.0.
    //#--- Do this in the following way:
    //#--- 1. get or compute the flipped HRF
    //#--- 2. use the signal list, which contains a
    //#---     sample for every time increment,
    //#---     as does the HRF.
    //#--- 3. zeropad signal on both ends for convolution
    //#--- 4. convolve HRF with data
    //#--- 5. save signal   
    
	
	//#--- 1: compute or get flipped HRF
	vtkFloatArray* Signal=this->ComputeBoxCar(this->onset,this->duration);
	int HRFsamps=this->flipHRF->GetNumberOfTuples();
	
	int posvals=HRFsamps;//length of postive value of HRF function
	int negvals=posvals;
	int len = Signal->GetNumberOfTuples();//length of signal
	int condition=Signal->GetNumberOfComponents();

	//construct outputarray
	//float* convResult=new float [len+posvals];
	vtkFloatArray* convResult=vtkFloatArray::New();
	//convResult->SetNumberOfTuples(len+posvals);
	convResult->SetNumberOfComponents(condition);

	for(int con=0;con<condition;con++)
	{
		//#--- 2&3: make zeropadded signal list
		float* data=new float [negvals+posvals+len];
		int data_count=0;
		for(int t=0;t<posvals;t++)
		{
			//zero pad
			data[data_count++]=0.0;
		}

		for(int t=0;t<len;t++)
		{
			//#--- data from signal
			//here we define Boxcar
			data[data_count++]=Signal->GetComponent(t,con);//[t];
		}

		for(int t=0;t<posvals;t++)
		{
			//zero pad
			data[data_count++]=0.0;
		}

		//#--- 4: convolve
		//#--- start at first samp of $data; shift $HRF;
		//#--- compute imghit + $posvals shifts and function mults;
		//#--- append each function multiply to $convResult.
		int shift = 0;
		float sval,hval,sum=0;
	
		int convResult_count=0;
		for(int y =0;y<len+posvals;y++)
		{
			sum=0;
			for(int s = 0;s<HRFsamps;s++)
			{
				sval = data[s+shift];
				hval = this->flipHRF->GetComponent(s,0);//[s];
				sum  = sum+(sval*hval);
			}
			if(y<len)
			convResult->InsertComponent(convResult_count++,con,sum);
			shift++;
		}
		delete[] data;
	}
	fileout(convResult,"conv.txt");
	fileout(this->Rescale(this->GaussianDownsampleList(convResult)),"convRes.txt");
	//#--- 5: save signal
	//output: convResult!!!!!!!!!!
	return this->Rescale(this->GaussianDownsampleList(convResult));
}
//compute HRF func data at single time point
vtkFloatArray* SignalModeling::ComputeHRF()
{
	float time_increment =   0.1;				//set tinc $::fMRIModelView(Design,Run$r,TimeIncrement) 
    //bold signal last for 30s
	float seemsEnough    =   30.0/time_increment;    //set seemsEnough [ expr 30.0 / $tinc ]
    int   HRF_Samples    =   int(seemsEnough);   //set HRFsamps [ expr round ($seemsEnough) ]
    
	//parameters for HRF 
	float a1=6;	//set a1 6
    float a2=12;    //set a2 12
    float b1=0.9;    //set b1 0.9
    float b2=0.9;    //set b2 0.9
    float c =0.35;    //set c 0.35
    float d1=a1*b1;    //set d1 [ expr $a1 * $b1 ]
    float d2=a2*b2;    //set d2 [ expr $a2 * $b2 ]

	//#--- compute first gamma function g1
    float*  g1 =new float[HRF_Samples];		//array for gamma1 function
	float   v = 0;		//static variable
	float   t = 0;	//
	
	for(int x=0;x<HRF_Samples;x++)
	{
		v     = pow(t/d1,d1)*exp((d1-t)/b1);
		g1[x] = v;
		t     = t+time_increment;
	}
	

	 //#--- compute second gamma function g2
     float* g2=new float[HRF_Samples];		//array for gamma2 function
	 v=0;
	 t=0;
	
	 for(int x=0;x<HRF_Samples;x++)
	 {
		 v     = pow(t/d2,d2)*exp((d2-t)/b2);
		 g2[x] = v;
		 t     = t+time_increment;
	 }

	 //#--- set h as difference of g1 and g2
     v=0;
	 t=0;
	 float max = -100000.0;
	 float min =  100000.0;
	 float v1,v2=0;
	 //array holding HRF point data
	 float* HRF = new float [HRF_Samples];
	 for(int x=0;x<HRF_Samples;x++)
	 {
		 v1 = g1[x];
		 v2 = g2[x];
		 v = v1-(c*v2);
		 HRF[x]=v;
		 if(v>max)
		 {
			 max=v;
		 }
		 if(v<min )
		 {
			 min=v;
		 }
	 }

	 delete [] g2;
	 delete [] g1;

	 //# flip HRF for convolution
	 //float* flipHRF=new float [HRF_Samples];
	 vtkFloatArray* flipHRF=vtkFloatArray::New();
	 flipHRF->SetNumberOfTuples(HRF_Samples);

	 int end = HRF_Samples;
	 for(int x=end;x>0;x--)
	 {
		 	v=HRF[x-1];
			flipHRF->SetComponent(HRF_Samples-x,0,v);
			//flipHRF[HRF_Samples-x]=v;
	 }
	 

	//write HRF to a txt file
	char* filename="flipHRF.txt";
	std::ofstream file (filename);
	if(file.is_open())
	{
		file<<"HRF point data"<<std::endl;
		for(int i=0;i<HRF_Samples;i++)
	   {
		 file<<HRF[i]<<std::endl;
	    }
	}
	file.close();
	delete [] HRF;
	return flipHRF;
}
//compute boxcar wave data at single time points
vtkFloatArray* SignalModeling::ComputeBoxCar(vtkFloatArray* onset,vtkFloatArray* duration)
{
	if(onset->GetNumberOfTuples()!=duration->GetNumberOfTuples())
	{

		std::cout<<"error in boxcar"<<std::endl;
	}
	
	//start volume: onset
	//duration is not seconds
	//it is number of TR: duration

	//time of repetition
	float TR=this->TR;
	//time_increment=0.1, as defined in ComputeHRF()
	float time_increment=this->TimeIncrement;
	//number of time point
	int numTimePoint=this->numTimePoints;
    //#--- what second of the signal should this footprint begin on?
    //#--- andhow many seconds should the signal footprintspan?
	
	int vox=int(TR/time_increment);

	int whole	 =  int(numTimePoint*TR/time_increment);
	//output array
	//holding boxcar stimuli signal
	vtkFloatArray* BoxCar=vtkFloatArray::New();
	BoxCar->SetNumberOfComponents(onset->GetNumberOfComponents());
	//zero out the whole array
	for(int j=0;j<onset->GetNumberOfComponents();j++) 
	{
		for(int i=0;i<whole;i++) 
		{BoxCar->InsertComponent(i,j,0.0);}
	}
	std::cout<<BoxCar->GetNumberOfComponents()<<std::endl;
	//#--- compute a boxcar signal footprint and insert into signal list
    //#--- boxcar signal goes from 0.0 to 1.0
    //#--- if event-related or mixed design, footprint might be delta function.
    //#--- In that case, we molde duration as one TimeIncrement.
	for(int i=0;i<onset->GetNumberOfComponents();i++)//number of components = numberof conditions
	{
		for(int j=0;j<onset->GetNumberOfTuples();j++)
		{
			for(int t=onset->GetComponent(j,i)*vox;t<onset->GetComponent(j,i)*vox+duration->GetComponent(j,i)*vox;t++)
			{
				BoxCar->InsertComponent(t,i,1.0);
			}
		}
	}
	fileout(BoxCar,"boxcar.txt");
	return BoxCar;
}
//compute baseline
vtkFloatArray* SignalModeling::Baseline()
{
	float min=0.0;
	float max=1.0;
	int numTimePoints=this->numTimePoints;

	//#--- build a constant basis function to capture baseline.
    //#--- signal first
	//float* baseline=new float[signaldim];
	vtkFloatArray* baseline=vtkFloatArray::New();
	baseline->SetNumberOfTuples(numTimePoints);
	for(int y=0;y<numTimePoints;y++)
	{
		baseline->SetComponent(y,0,1.0);
	}
	fileout(baseline,"baseline.txt");
	return baseline;
}

//function as name described
//invoked by AddDiscreteCosineBasis()
vtkFloatArray* SignalModeling::BuildDiscreteCosineBasis()
{
    vtkFloatArray* cosine=vtkFloatArray::New();

	//#--- cos ( (PI t / 2N) * (2k + 1) where k is the cos.
    int siglen=(this->numTimePoints) * (this->TR) / (this->TimeIncrement);
	float inc=this->TimeIncrement;
	int N    = this->numTimePoints*this->TR;
	int basiscount=this->FindNumberOfCosineBasis();
	cosine->SetNumberOfComponents(basiscount);
	//set  siglen [ llength $::fMRIModelView(Data,Run$r,EV$evnum,Signal) ]
    //set inc $::fMRIModelView(Design,Run$r,TimeIncrement)
    //set N [ expr $::fMRIModelView(Design,Run$r,numTimePoints) * \
     //           $::fMRIModelView(Design,Run$r,TR) ]
    float t = 0.0;
	float v = 0.0;
	for(int k=1;k<=basiscount;k++)
	{
		for ( int y = 0 ;y < siglen ; y++ )
		{
			v = cos (3.14159 * (2.0 * t + 1.0) * k / (2.0 * N));

			if ( t == 0.0 )
			{
				v = v / sqrt (2.0);
			}
			cosine->InsertComponent(y,k-1,v);
			t = t + inc ;
		}
	}
	fileout(cosine,"cosine1.txt");
	vtkFloatArray* cosine2=vtkFloatArray::New();
	cosine2->SetNumberOfComponents(basiscount);
	for(int j=0;j<basiscount;j++)
	{
		for(int i=0;i<this->numTimePoints;i++)
		{
			int t=i*(this->TR)/(this->TimeIncrement);
			float v=cosine->GetComponent(t,j);
			cosine2->InsertComponent(i,j,v);
		}
	}
	cosine->Delete();
	fileout(cosine2,"cosine.txt");
	return cosine2;
}

//function as name described
//invoked by AddDiscreteCosineBasis()
//return number of consine basises
int SignalModeling::FindNumberOfCosineBasis()
{
	//flag for use of DCBasis or not
	//1: use ECBasis, by default
	int UseDCBasis=1;

	if(UseDCBasis==1)
	{
		int     N =this->numTimePoints;
		float 	TC=this->ComputeDefaultHighpassTemporalCutoff();
		float   fc=1.0/TC;

		//#--- think this is right...
        //#--- have cos (PI*u(2t+1) / 2N) as the basic basis functions ...
        //#--- Tcutoff/TR is the cutoff period in samples.
        //#--- 2PI TR / Tcutoff is max cutoff omega;
        //#--- so to find out how many basis functions we need,
        //#--- take (PI*u(2t+1)/2N) and pull out omega,
        //#--- set it equal to cutoff omega, and solve for u.
        //#--- u is the number of frequencies (basis functions) we need.
        //#--- get something like this.
		float k = int(2.0*TR*N/TC);
		//#--- because we don't want the DC term (baseline models this...)
        k=k-1;
		this->NumberOfCosineBasis=k;
		return k;
	}
	else
	{
		this->NumberOfCosineBasis=0;
		return 0;
	}
}
//func as name
//invoked by FindNumberOfCosineBasis
//#--- Here's how the default cutoff frequency is computed:
//#--- Presume T = 1/f_lowest, is the longest epoch spacing in the run.
//#--- fmin is the cutoff frequency of the high-pass filter (lowest frequency
//#--- that we let pass through. Choose to let fmin = 0.666666/T (just less than the
//#--- lowest frequency in paradigm). As recommended in S.M. Smith, "Preparing fMRI
//#--- data for statistical analysis, in 'Functional MRI, an introduction to methods', Jezzard,
//#--- Matthews, Smith Eds. 2002, Oxford University Press.
float SignalModeling::ComputeDefaultHighpassTemporalCutoff()
{
	int numVolume=(this->numTimePoints-1);//=30 for test
	TR=this->TR;
	int longestEpoch=0;
	
	//find longestEpoch
	int diff=0;
	for(int j=0;j<this->onset->GetNumberOfComponents();j++)
	{
		int thisonset=0;
		int lastonset=0;
		for(int i=0;i<this->onset->GetNumberOfTuples();i++)
		{
			thisonset=this->onset->GetComponent(i,j);
			diff     =TR*(thisonset-lastonset);
			if(diff>longestEpoch)
			{
				longestEpoch=diff;
			}
			lastonset=thisonset;
		}
	}

	//#--- T is the number of seconds in the longest epoch
	float T=longestEpoch;
	//#--- set the model parameter, cutoff Period
	float HighpassCutoff=1.5*T;
	return HighpassCutoff;
}

vtkFloatArray* SignalModeling::GaussianDownsampleList(vtkFloatArray* inputList)
{
    vtkFloatArray* Gkernel=ComputeGaussianFilter();
	int numsamps=Gkernel->GetNumberOfTuples();// length of Gaussian kernel
	int numComponent=inputList->GetNumberOfComponents();
	//#--- takes a list in, subsamples it to a new length
    //#--- and returns the new list.
    //#---
    //#--- get or generate filter and find out its length
    //#--- *notice we are expecting inc to be an integer!
    //#--- so we never land between pixels when downsampling.

	int half = (int) (numsamps/2);//int half = floor(numsamps/2);
	
	int olen=inputList->GetNumberOfTuples();////signal length : numTimePoint*TR/Time_Increment
	int nlen=this->numTimePoints;//event length : number of timepoint
	
	int inc = olen/nlen;

	//#--- what if no downsampling is required!
    if ( inc == 1.0 ) 
	{
        return inputList;
    }

	vtkFloatArray* evlist=vtkFloatArray::New();
	evlist->SetNumberOfComponents(numComponent);
	evlist->SetNumberOfTuples(nlen);
	
	for(int x=0;x<numComponent;x++)
	{
		int count=0;
		//#---filter and subsample
		for ( int t = 0 ;  t < olen ; t = t + inc)
		{
			int start,stop=0;

			if ( t < half  )
			{
				start=half - t;
				stop =numsamps - 1;
			} 
			else if ( t > olen - (half + 1)  ) 
			{
				start=0;
				stop =half + (olen - 1 - t ) ;
			} 
			else 
			{
				start=0;
				stop =numsamps - 1;
			}
			float sum = 0.0;
			for (int j=start;j <= stop; j++ ) 
			{
				int k = int( t - half + j ) ;
				j = int (j);
				float v1 = inputList->GetComponent(k,x);
				float v2 = Gkernel->GetComponent(j,0);
				sum = sum + (v1 * v2) ;
			}
			evlist->SetComponent(count++,x,(sum/numsamps));
		}
	}
    return evlist;
}

vtkFloatArray* SignalModeling::ComputeGaussianFilter()
{
	float TR=this->TR;
	float Time_increment=this->TimeIncrement;
	//#--- Computes a gaussian kernel for convolution
    //#--- Define the filter's cutoff frequency fmax = 1/2*TR,
    //#--- or to be more careful, 1/2.01*TR,
    //#--- and wmax = 2pi * fmax = numsigmas*sigma.
    //#--- where numsigmas = 2.0 or 3.0,
    //#--- and numsigmas * sigma defines the half-band
    //#--- in the frequency domain.
    //#---
    //#--- use g(t) = 1/(sqrt(2pi)sigma) * exp ( -t^2 / 2sigma^2)
    //#--- sigma = 2pi * fmax / numsigmas.
    //#--- sigma = 2pi / (2.01*TR*numsigmas).
    //#---
    //#--- Assumes that all explanatory variables within a run
    //#--- have the same TR.
    //#---
    //#--- The signal we're filtering has sampling freq = fs = 1/0.1sec
    //#--- downsampling to a signal with sampling freq 1/TRsec.
    //#--- to signal with fmax = 1/2*TRsec
	float nyquistbuffer=2.001;
	float PI=3.14159265;
	float fmax=( 1.0 / (nyquistbuffer * TR));
	
	//#--- use 3 sigmas out for the kernel size (cutoff) now, 
    //#--- where gaussian approaches zero...
     //float  numsigmas = 3.0;
	 
	 //#--- try setting 2PI fmax = FWHM of the gaussian.
     //#--- so sigma = FWHM/(2*sqrt(2ln(2)))
     float numsigmas= 2.3548;
     float sigma=2.0 * PI * fmax / numsigmas ;

	 //#--- how many samples of the time-domain kernel do
     //#--- we need? 
     float inc = Time_increment;
     float numsamps= TR + inc ;

	 //#--- now compute the gaussian to convolve with.

	 float twoSigmaSq = sigma * sigma * 2.0;
     float twoPI = PI * PI;
     float sqrtTwoPI = sqrt(twoPI) ;

	 vtkFloatArray *kernel=vtkFloatArray::New();
	 //float* kernel;
	 //kernel->SetNumberOfTuples(numsamps);
	 kernel->SetNumberOfComponents(1);
	 
	 int kernel_count=0;
	 for(float t=-numsamps;t<=numsamps;t=t+inc)
	 {
		 float v;
		 v=(1.0/(sigma*sqrtTwoPI))*(exp(-(t*t)/twoSigmaSq));
		 kernel->InsertComponent(kernel_count,0,v);
		// kernel->SetComponent(kernel_count,0,v);
		 kernel_count++;
	 }

	 return kernel;
}

vtkFloatArray* SignalModeling::Rescale(vtkFloatArray* Inputlist)
{
	vtkFloatArray* Outputlist=vtkFloatArray::New();
	Outputlist->SetNumberOfComponents(Inputlist->GetNumberOfComponents());
	for(int c=0;c<Inputlist->GetNumberOfComponents();c++)
	{
		float max = -1000000.0;
		float min =  100000.0;
		int   dim =  Inputlist->GetNumberOfTuples();
		float temp_v=0;
		for (int t=0;t<dim;t++)//{ set t 0 } { $t < $dim } { incr t } 
		{
			temp_v = Inputlist->GetComponent(t,c);//[ lindex $data $t ]
			if ( temp_v > max ) 
			{
				max = temp_v;
			}
			if ( temp_v < min ) 
			{
				min = temp_v;
			}
		}

		//set min [ expr double($min) ]
		//set max [ expr double($max) ]
		float posrange,negrange=0;
		if ( max >= 0.0 ) 
		{
			posrange = max;
		} else 
		{
			posrange = 0.0;
		}
		if ( min <= 0.0 ) 
		{
			negrange = min;
		} 
		else 
		{
			negrange = 0.0;
		}
		float absprange = abs (posrange);
		float absnrange = abs (negrange);
		float range=0;

		if ( absprange > absnrange ) 
		{
			//#--- normalize to positive half
			range = absprange;
		} 
		else if ( absnrange > absprange ) 
		{
			//#--- normalize to negative half
			range = absnrange;
		} 
		else if (absnrange == absprange)
		{
			//#--- either value will work
			range = absprange;
		}
		//#--- normalize
		for ( int i = 0 ;i <dim;i++ )
		{
			float v = Inputlist->GetComponent(i,c);
			float nv=0;
			if (range != 0.0 ) 
			{
				nv = v/range ;
			} else 
			{
				nv = 0.0;
			}
			Outputlist->InsertComponent(i,c,nv);
		}
	}
    return Outputlist;
};

vtkFloatArray* GenerateDesignMatrix(SignalModeling* siglm)
{
	vtkFloatArray* out=vtkFloatArray::New();

	vtkFloatArray* Conv=siglm->ConvolveWithHRF();
	vtkFloatArray* Base=siglm->Baseline();
	vtkFloatArray* DCb =siglm->BuildDiscreteCosineBasis();

	int c=Conv->GetNumberOfComponents();
	int b=Base->GetNumberOfComponents();
	int d=DCb->GetNumberOfComponents();
	int len=Conv->GetNumberOfTuples();
	out->SetNumberOfComponents(c+b+d);
	//write to matrix
	for(int j=0;j<c;j++)
	{
		for(int i=0;i<len;i++)
		{
			out->InsertComponent(i,j,Conv->GetComponent(i,j));
		}
	}
	for(int j=c;j<c+b;j++)
	{
		for(int i=0;i<len;i++)
		{
			out->InsertComponent(i,j,Base->GetComponent(i,j-c));
		}
	}
	for(int j=b+c;j<b+c+d;j++)
	{
		for(int i=0;i<len;i++)
		{
			out->InsertComponent(i,j,DCb->GetComponent(i,j-b-c));
		}
	}
	return out;
}

void fileout(vtkFloatArray* input,char* name)
{
	std::ofstream file (name);
	
	int tuples=input->GetNumberOfTuples();
	int comp  =input->GetNumberOfComponents();
	
	if(file.is_open())
	{
		file<<name<<"point data"<<std::endl;
		file<<"number of tuples:"<<tuples<<endl;
		file<<"number of component:"<<comp<<endl;

		for(int j=0;j<comp;j++)
		{
			file<<"\ncomponent:"<<j<<"\n";
			for(int i=0;i<tuples;i++)
		   {
			 file<<input->GetComponent(i,j)<<"\,";
			}
		}
	}
	file.close();
}



