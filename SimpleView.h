/*=========================================================================
  Program:   Visualization Toolkit
  Module:    SimpleView.h
  Language:  C++
  Copyright 2009 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.
=========================================================================*/
#ifndef SimpleView_H
#define SimpleView_H

#include "vtkSmartPointer.h"    // Required for smart pointer internal ivars.
#include <QMainWindow>
#include "C_fileIO.h"
#include "stdio.h"

/*define*/
//file using
#define OUT_FILE_DIR                "OUTFILE/"
#define FILE_SAVE                   OUT_FILE_DIR"KneeOut.bmp"
#define FILE_SOBEL                  OUT_FILE_DIR"KneeOut_sobel.bmp"
#define FILE_SOBEL_RED              OUT_FILE_DIR"KneeOut_sobel_red.bmp"
#define FILE_MRF_MERGE              OUT_FILE_DIR"KneeOut_merge.bmp"
#define FILE_REGION_GROWING         OUT_FILE_DIR"region_growing.bmp"
#define FILE_SPLINE                 OUT_FILE_DIR"spline_out.bmp"
#define FILE_SPLINE_SAMPLE          OUT_FILE_DIR"spline_sample.bmp"
#define FILE_REGION_GROWING_TOP     OUT_FILE_DIR"region_growing_top.bmp"
#define FILE_REGION_GROWING_BOT     OUT_FILE_DIR"region_growing_bot.bmp"
#define FILE_REMOVE_FRAGMENT        OUT_FILE_DIR"remove_fragment.bmp"
#define FILE_SMOOTH_EDGE            OUT_FILE_DIR"smooth_edge.bmp"
#define FILE_MRF                    OUT_FILE_DIR"KneeOut_MRF.bmp"
#define FILE_THICKNESS              OUT_FILE_DIR"KneeOut_Thickness.bmp"

#if 0//set rgb pixel
typedef itk::RGBPixel<unsigned char> RGBPixelType;
typedef itk::Image<RGBPixelType,2> ImageType;
#else
typedef itk::Image<unsigned char,2> ImageType;
#endif
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileWriter<ImageType> WriterType;
static const int  dimension2D= 2;
typedef itk::Image<float, dimension2D >  float2D;
typedef itk::Image<unsigned char, dimension2D> uchar2D;

/*for MyView using*/
typedef enum
{
    None,
    CalculateDistance,
    RegionGrowing
}ViewMode;

// Forward Qt class declarations
class Ui_SimpleView;

// Forward VTK class declarations
class vtkQtTableView;
class vtkImageViewer2;
class vtkImageData;

class SimpleView : public QMainWindow
{
  Q_OBJECT

public:

  // Constructor/Destructor
  SimpleView();
  ~SimpleView();

public slots:

  virtual void slotOpenFile();
  virtual void slotExit();
  virtual void slotPreProcessMrf16();
  virtual void slotRunMrf();
  virtual void slotRunSig();
  virtual void slotRunAD();
  virtual void slotReset();
  virtual void slotWriteFile();
  virtual void slotTest();
  virtual void slotSpline();
  virtual void slotMerge();
  virtual void slotSobel();
  virtual void slotDeFragment();
  virtual void slotSmoothEdge();
  virtual void displayImage(vtkImageData *image);
  virtual void displayImage2(vtkImageData *image);
  virtual void displayMyView(QImage img, ViewMode view_mode);
  virtual void updatePixInfo(QString pix_info);
protected:

protected slots:

private:

  static void handle_double_click(vtkObject* obj, unsigned long event, void* ClientData, void* CallData);
  bool up_pixel_same(ImageType::Pointer seg_image, ImageType::IndexType pixelIndex, short pixel_value);
  QString getDICOMtagValue(itk::GDCMImageIO::Pointer gdcmIO, QString tagkey);
  vtkSmartPointer<vtkQtTableView> TableView;
  void region_growing(QString image_file, QString out_file, int seed_x, int seed_y, int replaced_pixel);
  void Opening();
  void drawThickness();
  void drawSpline();
  void drawSpline1();
  void drawSpline2();
  void RemoveFragments();
  int getThicknessPoint(int line2_x);
  QString getDistanceInfo(int x1, int y1, int x2, int y2);
  QString getDistanceInfo(int x);
  void debugImage(uchar2D::Pointer inputImage);
  void MrfMerge(int classes);
  int sobelFilter();

  //vtkImageViewer2 *viewer;

  // Designer form
  Ui_SimpleView *ui;
  C_fileIO myImage2D;
  uchar2D::Pointer OutImage;
  QString InputFile;
  int ImgW, ImgH;
  //for recording two red splines' Y cordinates, because X will be sequential from 0 to width.
  int SplineY1[1000], SplineY2[1000];
  int lowestY_x1, lowestY_x2;
//  std::vector<int> SplineY1, SplineY2;

public:
//  static   C_fileIO myImage2D;

};

#endif // SimpleView_H
