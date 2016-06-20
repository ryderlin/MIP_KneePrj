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


/************************************************************************/
/*                          Image dot color define                      */
/************************************************************************/
#define COLOR_LEFT_DOT_1 qRgb(0,255,0)
#define COLOR_LEFT_DOT_2 qRgb(0,255,255)
#define COLOR_CENTER_DOT_1 qRgb(255,0,0)
#define COLOR_CENTER_DOT_2 qRgb(255,0,255)
#define COLOR_RIGHT_DOT_1 qRgb(0,0,255)
#define COLOR_RIGHT_DOT_2 qRgb(255,255,0)


#endif // GLOBAL_TYPEDEF_H

