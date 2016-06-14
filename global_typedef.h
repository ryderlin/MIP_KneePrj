#ifndef GLOBAL_TYPEDEF_H
#define GLOBAL_TYPEDEF_H

//#define DEBUG_DRAW_THICKNESS


/************************************/
/*        pixel type                */
/************************************/
typedef signed short                            ssPixelType;
typedef unsigned short                          usPixelType;
typedef unsigned char                           ucPixelType;
typedef float                                   flPixelType;
typedef itk::RGBPixel <ucPixelType>             rgbPixelType;

/************************************************************************/
/*                          Image Type                                  */
/************************************************************************/
const unsigned int DIMENSION_2D = 2;
const unsigned int DIMENSION_3D = 3;

typedef itk::Image <ssPixelType,  DIMENSION_3D> ss3DImageType;

typedef itk::Image <usPixelType,  DIMENSION_2D> us2DImageType;
typedef itk::Image <ssPixelType,  DIMENSION_2D> ss2DImageType;
typedef itk::Image <ucPixelType,  DIMENSION_2D> uc2DImageType;
typedef itk::Image <flPixelType,  DIMENSION_2D> fl2DImageType;

typedef itk::Image <rgbPixelType>               rgbImageType;
typedef itk::Image <rgbPixelType, DIMENSION_2D> rgb2DImageType;

#endif // GLOBAL_TYPEDEF_H

