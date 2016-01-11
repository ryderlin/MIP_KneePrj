#include "C_fileIO.h"

C_fileIO::C_fileIO()
{
}

C_fileIO::~C_fileIO()
{
}

void C_fileIO::readFiletoImages( const char* dataPath )
{
	C_fileIO::FileReaderuchar2D::Pointer reader = C_fileIO::FileReaderuchar2D::New( );
	reader->SetFileName( dataPath);
	reader->Update();
	m_oriImage=reader->GetOutput();
}

void C_fileIO::writeImageToFile( const char* dataPath )
{
	FileWriteruchar2D::Pointer writer = FileWriteruchar2D::New();
	writer->SetInput( m_oriImage );
	writer->SetFileName( "./outSucess.png" );
	writer->Update();
}

void C_fileIO::readFromOtherImage (uchar2D::Pointer image)
{
	typedef itk::ImageDuplicator< uchar2D > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(image);
	duplicator->Update();
	m_oriImage = duplicator->GetOutput();

}

void C_fileIO::readFromOtherOutput(uchar2D::Pointer image)
{
    m_oriImage = image;
}

C_fileIO::uchar2D::Pointer C_fileIO::originalImages( )
{
	return m_oriImage;
}

vtkSmartPointer< vtkRenderer > C_fileIO::vtkRender()
{
	C_fileIO::toVTKuchar2D::Pointer connector = C_fileIO::toVTKuchar2D::New( );
	vtkSmartPointer< vtkImageData> vtkImage = vtkSmartPointer< vtkImageData  >::New();
	vtkSmartPointer< vtkImageActor >  actor = vtkSmartPointer< vtkImageActor >::New();
	vtkSmartPointer< vtkRenderer > renderer = vtkSmartPointer< vtkRenderer   >::New();
	
	connector->GetExporter( )->SetInput( m_oriImage );
	connector->GetImporter( )->Update( );

	vtkImage->Initialize( );
	vtkImage->DeepCopy( connector->GetImporter()->GetOutput() );
	printf("here");
	actor->SetInputData( vtkImage );
	actor->Update();

	actor->RotateX(180);

	renderer->AddActor2D( actor );
	//renderer->AddActor( actor );
	renderer->ResetCamera( );

	return renderer;

}

//cv::Mat C_fileIO::getMat()
//{
//	cv::Mat img;
//	img = itk::OpenCVImageBridge::ITKImageToCVMat< uchar2D >( m_oriImage );
//	return img;
//}

C_fileIO::float2D::Pointer C_fileIO::castfloat2D()
{
  typedef unsigned char                     InputPixelType;
  typedef float							    OutputPixelType;
  typedef itk::Image< InputPixelType, 2 >   InputImageType;
  typedef itk::Image< OutputPixelType, 2 >  OutputImageType;

  typedef itk::RescaleIntensityImageFilter< InputImageType, InputImageType > RescaleType;
  RescaleType::Pointer rescale = RescaleType::New();
  rescale->SetInput( m_oriImage );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( itk::NumericTraits< InputPixelType >::max() );

  //cout<<itk::NumericTraits< InputPixelType >::max()<<endl;

  typedef itk::CastImageFilter< InputImageType, OutputImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( rescale->GetOutput() );

 
  try
    {
    filter->Update();

	return filter->GetOutput();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Error: " << e << std::endl;
    }


}

C_fileIO::signedShort2D::Pointer C_fileIO::castsignedShort2D()
{
  typedef unsigned char                     InputPixelType;
  typedef signed short						OutputPixelType;
  typedef itk::Image< InputPixelType, 2 >   InputImageType;
  typedef itk::Image< OutputPixelType, 2 >  OutputImageType;

  typedef itk::RescaleIntensityImageFilter< InputImageType, InputImageType > RescaleType;
  RescaleType::Pointer rescale = RescaleType::New();
  rescale->SetInput( m_oriImage );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( itk::NumericTraits< OutputPixelType >::max() );

  typedef itk::CastImageFilter< InputImageType, OutputImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( rescale->GetOutput() );

 
  try
    {
    filter->Update();
	return filter->GetOutput();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Error: " << e << std::endl;
    }

}





