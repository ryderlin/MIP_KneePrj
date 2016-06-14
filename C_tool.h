#ifndef _C_TOOL_H_
#define _C_TOOL_H_
//#include "itkImageToVTKImageFilter.h"
#include <stdio.h> //General
#include <iostream>	
#include <vector>
#include <cstdlib>
#include <cstring>
#include <string>
#include <QtWidgets>	//#include <QtGui/QMainWindow>

#include <itkImage.h>	// ITK 
#include <itkRGBPixel.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkImageDuplicator.h"
#include "itkGDCMImageIO.h"	//ITK 3D
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"

#include <vtkSmartPointer.h>	// VTK 
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkEventQtSlotConnect.h>
#include <vtkVersion.h>	// VTK 3D
#include <vtkSphere.h>
#include <vtkSampleFunction.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVolumeProperty.h>
#include <vtkCamera.h>
#include <vtkImageShiftScale.h>
#include <vtkImageData.h>
#include <vtkImageFlip.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkXMLImageDataReader.h>
#include <itkScalarImageKmeansImageFilter.h>	//Kmean
#include <itkComposeImageFilter.h>		//Markov
#include <itkRescaleIntensityImageFilter.h>
#include <itkMRFImageFilter.h>
#include <itkDistanceToCentroidMembershipFunction.h>
#include <itkMinimumDecisionRule.h>
#include <itkVectorGradientAnisotropicDiffusionImageFilter.h>//WaterShed
#include <itkVectorGradientMagnitudeImageFilter.h>
#include <itkWatershedImageFilter.h>
#include <itkVectorCastImageFilter.h>
#include <itkScalarToRGBPixelFunctor.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkNeighborhoodConnectedImageFilter.h>  //RegionGrowing
#include <itkCurvatureFlowImageFilter.h>
//#include "10_include\itkOpenCVImageBridge.h" /*opencv*/
//#include <cv.h>
#include "itkImageToVTKImageFilter.h"
#include <itkExtractImageFilter.h>
/*wide string*/
#include <locale>
//#include <codecvt>
#include <string>

class tickTock{
private:
	time_t tstart, tend;
	std::string name;
public:

	void start()
	{
		tstart = clock();
	}
	void end(std::string aString)
	{
		tend = clock();
		cout << aString << ":" << difftime(tend, tstart) / 1000 << " sec." << endl;
	}
};

class C_tool
{
public:
	typedef signed short signedShort;
	//2D
	static const int  dimension2D= 2;
	typedef itk::Image<float, dimension2D >  float2D;
	typedef itk::Image<unsigned char, dimension2D> uchar2D;
	typedef itk::Image<signed short, dimension2D >  signedShort2D;
	typedef itk::ImageFileReader< uchar2D > FileReaderuchar2D;
	typedef itk::ImageFileWriter< uchar2D > FileWriteruchar2D;
    typedef itk::ImageToVTKImageFilter< uchar2D > toVTKuchar2D;
	//3D
	static const int  dimension= 3;
	typedef itk::Image< float, dimension>  float3D;
	typedef itk::Image< unsigned char, dimension> uchar3D;
	typedef itk::Image< signed short, dimension>  signedShort3D;
	typedef itk::ImageFileReader< signedShort3D > FileReadersignedShort3D;

//	typedef itk::ImageToVTKImageFilter< signedShort3D > toVTKsignedShort3D;
//	typedef itk::ImageToVTKImageFilter< signedShort2D > toVTKsignedShort2D;

	C_tool();
	~C_tool();

protected:
	void showSize(signedShort3D::Pointer inputImage);
	void showOrigin(signedShort3D::Pointer inputImage);
	void showIndexStart(signedShort3D::Pointer inputImage);
	
//	void debug(double tagNumber);

//	template <typename inputImageType >
//	void debug( typename inputImageType::Pointer inputImage)
//	{
//		inputImageType::RegionType	region = inputImage->GetLargestPossibleRegion();
//		inputImageType::SizeType		size = region.GetSize();
//		inputImageType::IndexType	pixelIndex;

//		for (unsigned int i=0;i<size[0];i++)
//		{
//			for (unsigned int j=0;j<size[1];j++)
//			{
//				pixelIndex[0] =i;
//				pixelIndex[1] =j;
//				double pixelValue= (double) inputImage->GetPixel( pixelIndex );
//				cout<< pixelValue <<",";
//			}
//			cout<< endl;
//		}
//		system("PAUSE");

//	}

//	template <typename outputImageType,typename inputImageType >
//	typename outputImageType::Pointer imageTypeTrasnfer3D(typename inputImageType::Pointer inputImage  )
//	{
//		typedef itk::RescaleIntensityImageFilter< inputImageType, inputImageType > RescaleType;
//		RescaleType::Pointer rescale = RescaleType::New();
//		rescale->SetInput( inputImage );
//		rescale->SetOutputMinimum( 0 );

//		//outputImageType::PixelType pixelValue;

//		//if (pixelValue==float)
//		//{
//		//	rescale->SetOutputMaximum( 256 );
//		//}
//		//else
//		//{
//		//	rescale->SetOutputMaximum( itk::NumericTraits<pixelValue>::max() );
//		//}

//		if (std::is_same<outputImageType,signedShort3D>::value)
//			rescale->SetOutputMaximum( itk::NumericTraits<signed short>::max() );
//		else if (std::is_same<outputImageType,uchar3D>::value )
//			rescale->SetOutputMaximum( itk::NumericTraits<unsigned char>::max() );
//		else if (std::is_same<outputImageType,float3D>::value )
//			rescale->SetOutputMaximum( 256 );

//		typedef itk::CastImageFilter< inputImageType, outputImageType > FilterType;
//		FilterType::Pointer filter = FilterType::New();
//		filter->SetInput( rescale->GetOutput() );

//		try
//		{
//			filter->Update();
//			return filter->GetOutput();
//		}
//		catch( itk::ExceptionObject & e )
//		{
//			std::cerr << "Error: " << e << std::endl;
//		}
//	}

	std::string myreplace(std::string &s, const std::string &toReplace, const std::string &replaceWith);
	std::string getFolder(const std::string& str);
	std::string getFileName(const std::string& str);
	std::string getExtension(std::string str);
	bool isDicom(std::string str);
//	__int64 GetFileSize(std::wstring const &path);
	
	
	void makeFolder(std::string aPath);
	void getDir(std::string, std::vector<std::string> & f);

private:

};
#endif
