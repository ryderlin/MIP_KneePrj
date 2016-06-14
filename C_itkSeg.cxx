#include "C_itkSeg.h"

C_itkSeg::C_itkSeg()
{
	numberOfInitialClasses = 4;
	numberOfClasses = 4;
}

C_itkSeg::~C_itkSeg()
{
}

void C_itkSeg::setParameter(int aInteger)
{
	if (aInteger != 0){
		numberOfInitialClasses = aInteger;
		numberOfClasses = aInteger;
	}
}

C_itkSeg::uchar2D::Pointer C_itkSeg::kmeanMethod2D(C_itkSeg::signedShort2D::Pointer  inputImage )
{

	//./kMeansClustering input.jpg output.jpg 1 3 0 100 200
	const unsigned int useNonContiguousLabels = 1;
	const unsigned int argoffset = 5;
	std::vector<double> userMeans;

	for (unsigned k = 0; k < numberOfInitialClasses; k++)
	{
		const double userProvidedInitialMean = (255 / numberOfInitialClasses*k);
		userMeans.push_back(userProvidedInitialMean);
	}

	typedef itk::ScalarImageKmeansImageFilter< signedShort2D > KMeansFilterType;

	KMeansFilterType::Pointer kmeansFilter = KMeansFilterType::New();
	kmeansFilter->SetInput(inputImage);
	kmeansFilter->SetUseNonContiguousLabels(useNonContiguousLabels);

	for (unsigned k = 0; k < numberOfInitialClasses; k++)
	{
		kmeansFilter->AddClassWithInitialMean(userMeans[k]);
	}

	try
	{
		kmeansFilter->Update();
		return kmeansFilter->GetOutput();
	}
	catch (itk::ExceptionObject & excp)
	{
		std::cerr << "Problem encountered while writing ";
		std::cerr << excp << std::endl;
	}




 // //./kMeansClustering input.jpg output.jpg 1 3 0 100 200
 // const unsigned int useNonContiguousLabels = 1;
 // const unsigned int numberOfInitialClasses = 3;
 // const unsigned int argoffset = 5;
 // std::vector<double> userMeans;

 // int iarr[3] = {0, 100, 200};

 // for( unsigned k = 0; k < numberOfInitialClasses; k++ )
 //   {
 //   const double userProvidedInitialMean = iarr[k];
 //   userMeans.push_back(userProvidedInitialMean);
 //   }
 // typedef itk::Image< signed short, 2 > ImageType;
 // typedef itk::ScalarImageKmeansImageFilter< ImageType > KMeansFilterType;
 // KMeansFilterType::Pointer kmeansFilter = KMeansFilterType::New();
 // kmeansFilter->SetInput( inputImage );
 // kmeansFilter->SetUseNonContiguousLabels( useNonContiguousLabels );
 //
 // // initialize using the user input means
 //   for( unsigned k = 0; k < numberOfInitialClasses; k++ )
 //   {
 //   kmeansFilter->AddClassWithInitialMean( userMeans[k] );
 //   }
 //
 // // Create and setup a writer
 //
 // typedef KMeansFilterType::OutputImageType  OutputImageType;
 //
 // typedef itk::ImageFileWriter< OutputImageType > WriterType;

 //
 // // execut the pipeline
 // try
 //   {
 //   kmeansFilter->Update();
	//return kmeansFilter->GetOutput();
 //   }
 // catch( itk::ExceptionObject & excp )
 //   {
 //   std::cerr << "Problem encountered while writing ";
 //   std::cerr << excp << std::endl;
 //   }
	

}

C_itkSeg::uchar2D::Pointer C_itkSeg::markovMethod2D(C_itkSeg::signedShort2D::Pointer  inputImage,C_itkSeg::uchar2D::Pointer labelImage )
{
	/*   MRF   Markov Random Filter          */
	const unsigned int numberOfIterations = 30; //30
	const double       smoothingFactor    = 3;
	const unsigned int numberOfArgumentsBeforeMeans = 7;

	//50 3 3 14.8 91.6 134.9

	typedef signed short        PixelType;
	const unsigned int          Dimension = 2;
	typedef itk::Image<PixelType, Dimension > ImageType;

	typedef unsigned char       LabelPixelType;
	typedef itk::Image<LabelPixelType, Dimension > LabelImageType;

	typedef itk::FixedArray<LabelPixelType,1>  ArrayPixelType;
	typedef itk::Image< ArrayPixelType, Dimension > ArrayImageType;
	typedef itk::ComposeImageFilter<ImageType, ArrayImageType > ScalarToArrayFilterType;

	ScalarToArrayFilterType::Pointer scalarToArrayFilter = ScalarToArrayFilterType::New();
	scalarToArrayFilter->SetInput( inputImage );

	typedef itk::MRFImageFilter< ArrayImageType, LabelImageType > MRFFilterType;
	MRFFilterType::Pointer mrfFilter = MRFFilterType::New();
	mrfFilter->SetInput( scalarToArrayFilter->GetOutput() );

	mrfFilter->SetNumberOfClasses( numberOfClasses );
	mrfFilter->SetMaximumNumberOfIterations( numberOfIterations );
	mrfFilter->SetErrorTolerance( 1e-7 );  //1e-7

	mrfFilter->SetSmoothingFactor( smoothingFactor );

	typedef itk::ImageClassifierBase< ArrayImageType,LabelImageType >   SupervisedClassifierType;
	SupervisedClassifierType::Pointer classifier =	SupervisedClassifierType::New();

	typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;
	DecisionRuleType::Pointer  classifierDecisionRule = DecisionRuleType::New();
	classifier->SetDecisionRule( classifierDecisionRule.GetPointer() );

	typedef itk::Statistics::DistanceToCentroidMembershipFunction<	ArrayPixelType > MembershipFunctionType;
	typedef MembershipFunctionType::Pointer MembershipFunctionPointer;
	double meanDistance = 0;
	MembershipFunctionType::CentroidType centroid(1);

	std::map<int, double> pixelCal;
	std::map <int, double>::iterator m1_Iter;

	pixelCal= getMean2D( inputImage  , labelImage );
	if (pixelCal.size() < numberOfClasses)
	{
		return labelImage;
	}

	for ( m1_Iter = pixelCal.begin( ); m1_Iter != pixelCal.end( ); m1_Iter++ )
	{
		MembershipFunctionPointer membershipFunction =	MembershipFunctionType::New();
		centroid[0] = m1_Iter->second ;
		membershipFunction->SetCentroid( centroid );
		classifier->AddMembershipFunction( membershipFunction );
		meanDistance += static_cast< double > (centroid[0]);

	}

	if (numberOfClasses > 0)
	{
		meanDistance /= numberOfClasses;
	}
	else
	{
		std::cerr << "ERROR: numberOfClasses is 0" << std::endl;
	}

	mrfFilter->SetSmoothingFactor( smoothingFactor );
	mrfFilter->SetNeighborhoodRadius( 1 );

	std::vector< double > weights;
	//weights.push_back(1.5);
	//weights.push_back(2.0);
	//weights.push_back(1.5);
	//weights.push_back(2.0);
	//weights.push_back(0.0); // This is the central pixel
	//weights.push_back(2.0);
	//weights.push_back(1.5);
	//weights.push_back(2.0);
	//weights.push_back(1.5);

	double totalWeight = 0;
	for(std::vector< double >::const_iterator wcIt = weights.begin();
		wcIt != weights.end(); ++wcIt )
	{
		totalWeight += *wcIt;
	}
	for(std::vector< double >::iterator wIt = weights.begin();
		wIt != weights.end(); ++wIt )
	{
		*wIt = static_cast< double > ( (*wIt) * meanDistance / (2 * totalWeight));
	}
	mrfFilter->SetMRFNeighborhoodWeight( weights );

	mrfFilter->SetClassifier( classifier );

	typedef MRFFilterType::OutputImageType  OutputImageType;

	typedef itk::Image< unsigned char, Dimension > RescaledOutputImageType;
	typedef itk::RescaleIntensityImageFilter<
		OutputImageType, RescaledOutputImageType >   RescalerType;
	RescalerType::Pointer intensityRescaler = RescalerType::New();
	intensityRescaler->SetOutputMinimum(   0 );
	intensityRescaler->SetOutputMaximum( 255 );
	intensityRescaler->SetInput( mrfFilter->GetOutput() );

	try
	{
		intensityRescaler->Update();
		return intensityRescaler->GetOutput();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Problem encountered while writing ";
		std::cerr << excp << std::endl;

	}

	std::cout << "Number of Iterations : ";
	std::cout << mrfFilter->GetNumberOfIterations() << std::endl;
	std::cout << "Stop condition: " << std::endl;
	std::cout << "  (1) Maximum number of iterations " << std::endl;
	std::cout << "  (2) Error tolerance:  "  << std::endl;
	std::cout << mrfFilter->GetStopCondition() << std::endl;

}

std::map<int, double> C_itkSeg::getMean2D( C_itkSeg::signedShort2D::Pointer  inputImage,C_itkSeg::uchar2D::Pointer labelImage)
{

	typedef itk::Image<unsigned char, 2>  ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;


	ImageType::RegionType	region = labelImage->GetLargestPossibleRegion();
	ImageType::SizeType		size = region.GetSize();
	//ImageType::PixelType	pixelValue;
	ImageType::IndexType	pixelIndex;

	std::map<int, double> meanCalculate;
	std::map<int, int> count;
	std::map<int ,double >::iterator l_it;;


	for (unsigned int i=0;i<size[0];i++)
	{
		for (unsigned int j=0;j<size[1];j++)
		{
			pixelIndex[0] =i;
			pixelIndex[1] =j;

			int classLabel= (int) labelImage->GetPixel( pixelIndex );
			double pixelValue= (double) inputImage->GetPixel( pixelIndex );
			 
			 l_it=meanCalculate.find(classLabel);
			if(l_it==meanCalculate.end())
			{
				meanCalculate.insert ( std::pair <int, double>  ( classLabel, pixelValue ) );
				count.insert (  std::pair <int, int>  ( classLabel, 1 ) );
			}
			else
			{
				meanCalculate[classLabel]=meanCalculate[classLabel]+pixelValue;
				count[classLabel]=count[classLabel]+1;
			}

		}
	}

	for ( l_it = meanCalculate.begin( ); l_it != meanCalculate.end( ); l_it++ )
	{
		int classLabel = l_it->first ;
		meanCalculate[classLabel]=meanCalculate[classLabel]/count[classLabel];
	}
	return meanCalculate;
}

C_itkSeg::uchar2D::Pointer  C_itkSeg::watershedMethod2D(C_itkSeg::uchar2D::Pointer  inputImage,double waterShedSetting[20] )
{
	typedef itk::Vector< float, 3 >              VectorPixelType;
	typedef itk::Image< VectorPixelType, 2 >     VectorImageType;
	typedef itk::VectorGradientAnisotropicDiffusionImageFilter<VectorImageType, VectorImageType >DiffusionFilterType;
	DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();

	typedef itk::Image< float, 2 >               ScalarImageType;
	typedef itk::WatershedImageFilter< ScalarImageType >WatershedFilterType;
	WatershedFilterType::Pointer watershed = WatershedFilterType::New();

	typedef itk::VectorGradientMagnitudeImageFilter< VectorImageType >GradientMagnitudeFilterType;
	GradientMagnitudeFilterType::Pointer gradient = GradientMagnitudeFilterType::New();


	typedef itk::RGBPixel< unsigned char >       RGBPixelType;
	typedef itk::Image< RGBPixelType, 2 >        RGBImageType;

	typedef itk::Image<unsigned char, 2>  uchar2D;
	typedef itk::CastImageFilter< uchar2D , VectorImageType >CastFilterType;
	CastFilterType::Pointer caster = CastFilterType::New();


	//diffusion->SetConductanceParameter( 2.0 );
	//diffusion->SetNumberOfIterations( 10 );
	//watershed->SetThreshold( 0.001 );
	//watershed->SetLevel( 0.15 );
	//gradient->SetUsePrincipleComponents(0);
	//diffusion->SetTimeStep(0.125);

	diffusion->SetConductanceParameter( double(waterShedSetting[0]) );
	diffusion->SetNumberOfIterations( int(waterShedSetting[1]) );
	watershed->SetThreshold( double(waterShedSetting[2]) );
	watershed->SetLevel( double(waterShedSetting[3]) );
	gradient->SetUsePrincipleComponents(int(waterShedSetting[4]));
	diffusion->SetTimeStep(double(waterShedSetting[5]));

	caster->SetInput( inputImage );

	diffusion->SetInput(caster->GetOutput());
	gradient->SetInput(diffusion->GetOutput());
	watershed->SetInput(gradient->GetOutput());

	typedef itk::Image< itk::IdentifierType, 2 > LabeledImageType;
	typedef itk::CastImageFilter< LabeledImageType, uchar2D >CastToUCharType;
	CastToUCharType::Pointer castToUchar = CastToUCharType::New();

	castToUchar->SetInput(watershed->GetOutput());
	castToUchar->Update();
	return castToUchar->GetOutput();

}

C_itkSeg::uchar2D::Pointer  C_itkSeg::regionGrowMethod2D(C_itkSeg::uchar2D::Pointer  inputImage,double pSignalArray[20])
{

	typedef itk::Image<unsigned char,2>  uchar2D;
	uchar2D::IndexType seed1;

	cout << int(pSignalArray[0]) << endl;  
	cout << int(pSignalArray[1]) << endl;
	cout << int(pSignalArray[2]) << endl;
	cout << int(pSignalArray[3]) << endl;

	seed1[0] = (int)pSignalArray[0];
	seed1[1] = (int)pSignalArray[1];

	typedef itk::ConnectedThresholdImageFilter<uchar2D, uchar2D> ConnectedFilterType;
	ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();

	connectedThreshold->SetLower(pSignalArray[2]);
	connectedThreshold->SetUpper(pSignalArray[3]);
	connectedThreshold->SetReplaceValue(255);

	typedef itk::CurvatureFlowImageFilter<uchar2D, uchar2D> CurvatureFlowImageFilterType;
	CurvatureFlowImageFilterType::Pointer smoothing =CurvatureFlowImageFilterType::New();

	connectedThreshold->SetSeed(seed1);
	smoothing->SetInput(inputImage );
	connectedThreshold->SetInput(inputImage);
	connectedThreshold->Update();
    return  connectedThreshold->GetOutput();
}

C_itkSeg::uchar2D::Pointer  C_itkSeg::anisotropicDiffusion(C_itkSeg::uchar2D::Pointer  inputImage)
{
    typedef unsigned char InputPixelType;
    typedef itk::Image< InputPixelType, 2 > InputImageType;

    typedef float                                     OutputPixelType;
    typedef itk::Image< OutputPixelType, 2 >  OutputImageType;
    typedef itk::GradientAnisotropicDiffusionImageFilter< InputImageType,
            OutputImageType > FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(inputImage);
    filter->SetNumberOfIterations(5.0/*this->ui->leAD_Iteration->text().toDouble()*/);
    filter->SetTimeStep( 0.25 );
    filter->SetConductanceParameter(5.0/*this->ui->leAD_Conductance->text().toDouble()*/);

    //convert output from float to uchar2d
    typedef itk::RescaleIntensityImageFilter<OutputImageType, InputImageType> RescaleType;
    RescaleType::Pointer rescaler = RescaleType::New();
    rescaler->SetInput( filter->GetOutput() );
    rescaler->SetOutputMinimum( itk::NumericTraits< InputPixelType >::min() );
    rescaler->SetOutputMaximum( itk::NumericTraits< InputPixelType >::max() );
    return rescaler->GetOutput();
}



