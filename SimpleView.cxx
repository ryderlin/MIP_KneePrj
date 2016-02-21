/*
 * Copyright 2007 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */


#include "ui_SimpleView.h"
#include "SimpleView.h"
//QT
#include<QFileDialog>
#include<QDir>
//ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkGDCMImageIO.h"
#include "itkBMPImageIOFactory.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
//VTK

#include <vtkAutoInit.h>
#include <vtkDataObjectToTable.h>
#include <vtkElevationFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkQtTableView.h>
#include <vtkVectorText.h>

#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageViewer2.h>
#include <vtkImageFlip.h>
#include <vtkImageActor.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleImage.h>
#include "itkSigmoidImageFilter.h"

#include "vtkSmartPointer.h"
#include "C_itkSeg.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

typedef itk::Image<unsigned char/*RGBPixelType*/,2> ImageType;
// Constructor
SimpleView::SimpleView()
{
    VTK_MODULE_INIT(vtkRenderingOpenGL);
    VTK_MODULE_INIT(vtkInteractionStyle);
//    VTK_MODULE_INIT(vtkRenderingWindow);sfdfs
    this->ui = new Ui_SimpleView;
    this->ui->setupUi(this);

    // Set up action signals and slots
    connect(this->ui->actionOpenFile, SIGNAL(triggered()), this, SLOT(slotOpenFile()));
    connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));
    connect(this->ui->btnRunMrf, SIGNAL (released()),this, SLOT (slotRunMrf()));
    connect(this->ui->btnRunSig, SIGNAL (released()),this, SLOT (slotRunSig()));
    connect(this->ui->btnRunAD, SIGNAL (released()),this, SLOT (slotRunAD()));
    connect(this->ui->btnReset, SIGNAL (released()),this, SLOT (slotReset()));
    connect(this->ui->btnWriteFile, SIGNAL (released()),this, SLOT (slotWriteFile()));
    this->ui->qvtkWidget_Ori->repaint();
    this->ui->qvtkWidget_Seg->hide();
}

SimpleView::~SimpleView()
{
  // The smart pointers should clean up for up
  delete ui;
}

// Action to be taken upon file open
void SimpleView::slotOpenFile()
{
    itk::BMPImageIOFactory::RegisterOneFactory();
    InputFile = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());

#if 0 //for testing use
#else
//    typedef itk::RGBPixel<unsigned char> RGBPixelType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(InputFile.toLatin1().data());
    reader->Update();

    std::cout << "size is " << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] << " "
                << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] << std::endl;
    ImageType::IndexType pixelIndex;
    for (int c = 0; c < 800;  c++)
    {
        pixelIndex[0] = 580;
        pixelIndex[1] = c;
        reader->GetOutput()->SetPixel(pixelIndex, 255);
    }
//    reader->GetOutput()->GetPixel(pixelIndex, 255);
    myImage2D.readFromOtherOutput(reader->GetOutput());
#endif
    /***** 變數設定 *****/
    // image is grayscale
//    typedef itk::Image<unsigned short,2> ImageType; //要用unsigned short以上型態，dicom pixel value才不會爆掉
    //typedef itk::Image<signed int,2> ImageType;
    // image is RGB

    //typedef itk::ImageFileWriter<ImageType> WriterType;
    typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;


    // setup and connect itk with vtk
    ConnectorType::Pointer connector = ConnectorType::New();
    connector->SetInput(reader->GetOutput());
    connector->Update();

    // flip image in Y axis (翻轉image, 逆時針180旋轉)
    vtkSmartPointer<vtkImageFlip> flipYFilter = vtkSmartPointer<vtkImageFlip>::New();
    flipYFilter->SetFilteredAxis(1); // flip Y axis
    flipYFilter->SetInputData(connector->GetOutput());
    flipYFilter->Update();

    vtkSmartPointer<vtkImageData> vtkimage = vtkImageData::New();
    vtkimage->DeepCopy(flipYFilter->GetOutput());

    this->displayImage(vtkimage);

    reader = NULL;
    connector = NULL;
    flipYFilter = NULL;
    vtkimage = NULL;
}

void SimpleView::slotRunMrf()
{
    C_fileIO kmeanImage;
    C_fileIO mrfImage;
    C_itkSeg mySeg;
    mySeg.setParameter(this->ui->leIterationNum->text().toInt());
    kmeanImage.readFromOtherImage(mySeg.kmeanMethod2D( myImage2D.castsignedShort2D())) ;
    myImage2D.readFromOtherImage(  mySeg.markovMethod2D( myImage2D.castsignedShort2D(), kmeanImage.originalImages()  )  );
    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
    this->ui->qvtkWidget_Seg->repaint();
    this->ui->qvtkWidget_Seg->show();
}

void SimpleView::slotRunAD()
{
    typedef unsigned char InputPixelType;
    typedef itk::Image< InputPixelType, 2 > InputImageType;

//    typedef itk::ImageFileReader< InputImageType >  ReaderType;
//    ReaderType::Pointer reader = ReaderType::New();
//    reader->SetFileName(InputFile.toLatin1().data());

    typedef float                                     OutputPixelType;
    typedef itk::Image< OutputPixelType, 2 >  OutputImageType;
    typedef itk::GradientAnisotropicDiffusionImageFilter< InputImageType,
    OutputImageType > FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput( myImage2D.originalImages());
    filter->SetNumberOfIterations(this->ui->leAD_Iteration->text().toDouble());
    filter->SetTimeStep( 0.125 );
    filter->SetConductanceParameter(this->ui->leAD_Conductance->text().toDouble());

    typedef itk::RescaleIntensityImageFilter<OutputImageType, InputImageType> RescaleType;

    //convert output from float to uchar2d
    RescaleType::Pointer rescaler = RescaleType::New();
    rescaler->SetInput( filter->GetOutput() );
    rescaler->SetOutputMinimum( itk::NumericTraits< InputPixelType >::min() );
    rescaler->SetOutputMaximum( itk::NumericTraits< InputPixelType >::max() );

    myImage2D.readFromOtherOutput(rescaler->GetOutput());
    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
    this->ui->qvtkWidget_Seg->repaint();
    this->ui->qvtkWidget_Seg->show();
}

void SimpleView::slotRunSig()
{
#if 0 //for testing use
#else
    double alpha = this->ui->leSigAlpha->text().toDouble();
    double beta = this->ui->leSigBeta->text().toDouble();
    typedef itk::SigmoidImageFilter <ImageType, ImageType>
      SigmoidImageFilterType;
    SigmoidImageFilterType::Pointer sigmoidFilter
      = SigmoidImageFilterType::New();
    sigmoidFilter->SetInput(myImage2D.originalImages());
    sigmoidFilter->SetOutputMinimum(0);
    sigmoidFilter->SetOutputMaximum(itk::NumericTraits< signed short >::max());
    sigmoidFilter->SetAlpha(alpha);
    sigmoidFilter->SetBeta(beta);
    myImage2D.readFromOtherOutput(sigmoidFilter->GetOutput());
#endif
    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
    this->ui->qvtkWidget_Seg->repaint();
    this->ui->qvtkWidget_Seg->show();
}

void SimpleView::slotReset()
{
    this->ui->qvtkWidget_Seg->hide();
    myImage2D.readFiletoImages(InputFile.toStdString().c_str());
}

void SimpleView::slotWriteFile()
{
myImage2D.writeImageToFile("/home/ryderlin/Documents/KneeOut.bmp");
}

void SimpleView::slotExit() {
  qApp->exit();
}

void SimpleView::displayImage(vtkImageData *image)
{
    // Create image actor
    vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);

    // set actor properties
    actor->SetInputData(image);
    //actor->InterpolateOff();

    renderer->AddActor(actor);
    this->ui->qvtkWidget_Ori->SetRenderWindow(renderWindow);

    // window interactor style for display images
    vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
    // set interactor style to the qvtkWidget Interactor
    this->ui->qvtkWidget_Ori->GetInteractor()->SetInteractorStyle(style);

    this->ui->qvtkWidget_Ori->update();
}
