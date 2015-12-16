#ifndef _C_FILEIO_H_
#define _C_FILEIO_H_
#include "C_tool.h"
class C_fileIO : public C_tool
{
public:
	C_fileIO();
	~C_fileIO();
	void readFiletoImages( const char* dataPath );
	void readFromOtherImage (uchar2D::Pointer);
	void writeImageToFile( const char* dataPath );
	vtkSmartPointer< vtkRenderer > vtkRender();
//	cv::Mat getMat();
	uchar2D::Pointer originalImages( );
	float2D::Pointer castfloat2D();
	signedShort2D::Pointer castsignedShort2D();
private:
	uchar2D::Pointer m_oriImage;
};

#endif
