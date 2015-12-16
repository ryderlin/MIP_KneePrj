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
//VTK

#include <vtkDataObjectToTable.h>
#include <vtkElevationFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkQtTableView.h>
#include <vtkVectorText.h>

#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
//#include <vtkDICOMImageReader.h>
#include <vtkImageViewer2.h>
#include <vtkImageFlip.h>
#include <vtkImageActor.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleImage.h>

#include "vtkSmartPointer.h"
#include "C_itkSeg.h"
#include "C_fileIO.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

// Constructor
SimpleView::SimpleView()
{
  this->ui = new Ui_SimpleView;
  this->ui->setupUi(this);

  // Set up action signals and slots
  connect(this->ui->actionOpenFile, SIGNAL(triggered()), this, SLOT(slotOpenFile()));
  connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));

};

SimpleView::~SimpleView()
{
  // The smart pointers should clean up for up
  delete ui;
}

// Action to be taken upon file open
void SimpleView::slotOpenFile()
{
    C_fileIO myImage2D;

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());

    myImage2D.readFiletoImages("/media/sf_shared_host/1.bmp");
    C_fileIO kmeanImage;
    C_fileIO mrfImage;
    C_itkSeg mySeg;
    mySeg.setParameter(8);

    kmeanImage.readFromOtherImage(mySeg.kmeanMethod2D( myImage2D.castsignedShort2D())) ;
    mrfImage.readFromOtherImage(  mySeg.markovMethod2D( myImage2D.castsignedShort2D(), kmeanImage.originalImages()  )  );
    this->ui->qvtkWidget->GetRenderWindow()->AddRenderer(mrfImage.vtkRender());
    this->ui->qvtkWidget->repaint();
    return;
    /***** 變數設定 *****/
    // image is grayscale
//    typedef itk::Image<unsigned short,2> ImageType; //要用unsigned short以上型態，dicom pixel value才不會爆掉
    //typedef itk::Image<signed int,2> ImageType;
    // image is RGB
    typedef itk::RGBPixel<unsigned char> RGBPixelType;
    typedef itk::Image<RGBPixelType,2> ImageType;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    //typedef itk::ImageFileWriter<ImageType> WriterType;
    typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;

//    this->ui->lbOriImage->back
    /***** DICOM讀檔 *****/
    ReaderType::Pointer reader = ReaderType::New();
    //reader->SetImageIO( gdcmIO );
    //vtkDICOMImageReader *reader=vtkDICOMImageReader::New();
    reader->SetFileName(fileName.toLatin1().data());
    //reader->SetFileName("/Users/Chris/Downloads/ntu.jpg");
    reader->Update();
//    C_itkSeg mySeg;
//    mySeg.setParameter(5);
//    mySeg.kmeanMethod2D(reader->GetOutput());

    /***** DICOM pixel value *****/
    /*
      ImageType::Pointer image = reader->GetOutput();
      ImageType::IndexType pixelIndex;
      pixelIndex[0] = 400;   // x position
      pixelIndex[1] = 1600;   // y position
      std::cout << image->GetPixel(pixelIndex);
    */

    /***** DICOM寫檔 *****/
    /*
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName("/Users/Chris/Downloads/DICOM_IMG/test.dcm");
      writer->SetInput(reader->GetOutput());
      writer->Update();
    */

    // setup and connect itk with vtk
    ConnectorType::Pointer connector = ConnectorType::New();
    //connector->GetExporter()->SetInput(reader->GetOutput());
    //connector->GetImporter()->Update();
    connector->SetInput(reader->GetOutput());
    connector->Update();

    // flip image in Y axis (翻轉image, 逆時針180旋轉)
    vtkSmartPointer<vtkImageFlip> flipYFilter = vtkSmartPointer<vtkImageFlip>::New();
    flipYFilter->SetFilteredAxis(1); // flip Y axis
    flipYFilter->SetInputData(connector->GetOutput());
    flipYFilter->Update();

    vtkSmartPointer<vtkImageData> vtkimage = vtkImageData::New();
    //vtkimage->DeepCopy(connector->GetOutput());
    vtkimage->DeepCopy(flipYFilter->GetOutput());

    this->displayImage(vtkimage);

    reader = NULL;
    connector = NULL;
    flipYFilter = NULL;
    vtkimage = NULL;

    /* display image in vtkViewer
        vtkSmartPointer<vtkImageViewer2> viewer = vtkImageViewer2::New();
        //viewer = vtkImageViewer2::New();

        //set VTK Viewer to QVTKWidget in Qt's UI
        this->ui->qvtkWidget->SetRenderWindow(viewer->GetRenderWindow());
        viewer->SetupInteractor(this->ui->qvtkWidget->GetRenderWindow()->GetInteractor());

        //Set input image to VTK viewer
        viewer->SetInputData(vtkimage);
        //viewer->SetInputData(connector->GetOutput());
        viewer->GetRenderer()->ResetCamera();
        viewer->Render();
        //viewer->SetColorWindow(255);
        //viewer->SetColorLevel(128);
        this->ui->qvtkWidget->update();
    */
}

void SimpleView::slotExit() {
  qApp->exit();
}

void SimpleView::displayImage(vtkImageData *image)
{
    /*
    int *dim= image->GetDimensions();
    double *spacing = image->GetSpacing();
    double *origin = image->GetOrigin();

    float Cx = (dim[0] * spacing[0])/2. + origin[0];
    float Cy = (dim[1] * spacing[1])/2. + origin[1];
    */

    // Create image actor
    vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();
    // Create camera
    //vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
    // Create renderer and render window
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);

    /*
    camera->ParallelProjectionOn();
    camera->SetFocalPoint(Cx,Cy,0);
    camera->SetPosition(Cx,Cy,1);
    */
    //
    //    // to flip de image
    //    camera->SetViewUp (0, 1, 0);
    //
    // set actor properties
    actor->SetInputData(image);
    //actor->InterpolateOff();

    renderer->AddActor(actor);
    //renderer->SetActiveCamera(camera);
    //renderer->ResetCamera();

    this->ui->qvtkWidget->SetRenderWindow(renderWindow);

    // window interactor style for display images
    vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
    // set interactor style to the qvtkWidget Interactor
    this->ui->qvtkWidget->GetInteractor()->SetInteractorStyle(style);

    this->ui->qvtkWidget->update();
}
