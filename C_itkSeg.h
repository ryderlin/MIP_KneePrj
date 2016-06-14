#ifndef _C_ITKSEQ_H_
#define _C_ITKSEQ_H_
#include "C_tool.h"
class C_itkSeg: public C_tool
{
public:
	C_itkSeg();
	~C_itkSeg();
	void setParameter(int);

//protected:
	unsigned int numberOfInitialClasses, numberOfClasses;

	uchar2D::Pointer kmeanMethod2D(signedShort2D::Pointer);
	uchar2D::Pointer markovMethod2D(signedShort2D::Pointer, uchar2D::Pointer);
	std::map<int, double> getMean2D(signedShort2D::Pointer, uchar2D::Pointer);
	uchar2D::Pointer watershedMethod2D(uchar2D::Pointer, double*);
	uchar2D::Pointer regionGrowMethod2D(uchar2D::Pointer, double*);
    uchar2D::Pointer anisotropicDiffusion(uchar2D::Pointer  inputImage);
};
#endif
