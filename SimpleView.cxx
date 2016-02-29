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
#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkSigmoidImageFilter.h"
//VTK

#include <vtkAutoInit.h>
#include <vtkDataObjectToTable.h>
#include <vtkElevationFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkQtTableView.h>
#include <vtkVectorText.h>
#include <vtkPropPicker.h>
#include <vtkCornerAnnotation.h>
#include <vtkAssemblyPath.h>
#include <vtkMath.h>
#include <vtkTextProperty.h>

#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageViewer2.h>
#include <vtkImageFlip.h>
#include <vtkImageActor.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleImage.h>

#include "vtkSmartPointer.h"
#include "C_itkSeg.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

// Template for image value reading
template<typename T>
void vtkValueMessageTemplate(vtkImageData* image, int* position,
                             std::string& message)
{
  T* tuple = ((T*)image->GetScalarPointer(position));
  int components = image->GetNumberOfScalarComponents();
  for (int c = 0; c < components; ++c)
    {
    message += vtkVariant(tuple[c]).ToString();
    if (c != (components - 1))
      {
      message += ", ";
      }
    }
  message += " )";
}

// The mouse motion callback, to pick the image and recover pixel values
class vtkImageInteractionCallback : public vtkCommand
{
public:
  static vtkImageInteractionCallback *New()
    {
    return new vtkImageInteractionCallback;
    }

  vtkImageInteractionCallback()
    {
    this->Viewer     = NULL;
    this->Picker     = NULL;
    this->Annotation = NULL;
    }

  ~vtkImageInteractionCallback()
    {
    this->Viewer     = NULL;
    this->Picker     = NULL;
    this->Annotation = NULL;
    }

  void SetPicker(vtkPropPicker *picker)
    {
    this->Picker = picker;
    }

  void SetAnnotation(vtkCornerAnnotation *annotation)
    {
    this->Annotation = annotation;
    }

  void SetViewer(vtkImageViewer2 *viewer)
    {
    this->Viewer = viewer;
    }

  virtual void Execute(vtkObject *, unsigned long vtkNotUsed(event), void *)
    {
    vtkRenderWindowInteractor *interactor =
      this->Viewer->GetRenderWindow()->GetInteractor();
    vtkRenderer* renderer = this->Viewer->GetRenderer();
    vtkImageActor* actor = this->Viewer->GetImageActor();
    vtkImageData* image = this->Viewer->GetInput();
    vtkInteractorStyle *style = vtkInteractorStyle::SafeDownCast(
      interactor->GetInteractorStyle());

#if VTK_MAJOR_VERSION <= 5
    image->Update();
#endif

    // Pick at the mouse location provided by the interactor
    this->Picker->Pick(interactor->GetEventPosition()[0],
                       interactor->GetEventPosition()[1],
                       0.0, renderer);

    // There could be other props assigned to this picker, so
    // make sure we picked the image actor
    vtkAssemblyPath* path = this->Picker->GetPath();
    bool validPick = false;

    if (path)
      {
      vtkCollectionSimpleIterator sit;
      path->InitTraversal(sit);
      vtkAssemblyNode *node;
      for (int i = 0; i < path->GetNumberOfItems() && !validPick; ++i)
        {
        node = path->GetNextNode(sit);
        if (actor == vtkImageActor::SafeDownCast(node->GetViewProp()))
          {
          validPick = true;
          }
        }
      }

    if (!validPick)
      {
      this->Annotation->SetText(0, "Off Image");
      interactor->Render();
      // Pass the event further on
      style->OnMouseMove();
      return;
      }

    // Get the world coordinates of the pick
    double pos[3];
    this->Picker->GetPickPosition(pos);

    int image_coordinate[3];

    int axis = this->Viewer->GetSliceOrientation();
    switch (axis)
      {
      case vtkImageViewer2::SLICE_ORIENTATION_XZ:
        image_coordinate[0] = vtkMath::Round(pos[0]);
        image_coordinate[1] = this->Viewer->GetSlice();
        image_coordinate[2] = vtkMath::Round(pos[2]);
        break;
      case vtkImageViewer2::SLICE_ORIENTATION_YZ:
        image_coordinate[0] = this->Viewer->GetSlice();
        image_coordinate[1] = vtkMath::Round(pos[0]);
        image_coordinate[2] = vtkMath::Round(pos[1]);
        break;
      default:  // vtkImageViewer2::SLICE_ORIENTATION_XY
        image_coordinate[0] = vtkMath::Round(pos[0]);
        image_coordinate[1] = vtkMath::Round(pos[1]);
        image_coordinate[2] = this->Viewer->GetSlice();
        break;
      }

    std::string message = "Location: ( ";
    message += vtkVariant(image_coordinate[0]).ToString();
    message += ", ";
    message += vtkVariant(image_coordinate[1]).ToString();
    message += ", ";
    message += vtkVariant(image_coordinate[2]).ToString();
    message += " )\nValue: ( ";

    switch (image->GetScalarType())
      {
      vtkTemplateMacro((vtkValueMessageTemplate<VTK_TT>(image,
                                                        image_coordinate,
                                                        message)));

      default:
        return;
      }

    this->Annotation->SetText( 0, message.c_str() );
    interactor->Render();
    style->OnMouseMove();
    }

private:
  vtkImageViewer2*      Viewer;      // Pointer to the viewer
  vtkPropPicker*        Picker;      // Pointer to the picker
  vtkCornerAnnotation*  Annotation;  // Pointer to the annotation
};

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
    connect(this->ui->btnSeg, SIGNAL (released()),this, SLOT (slotSegmentation()));
    connect(this->ui->btnTest, SIGNAL (released()),this, SLOT (slotTest()));
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
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(InputFile.toLatin1().data());
    reader->Update();

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

    this->displayImage2(vtkimage);

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

    typedef float                                     OutputPixelType;
    typedef itk::Image< OutputPixelType, 2 >  OutputImageType;
    typedef itk::GradientAnisotropicDiffusionImageFilter< InputImageType,
    OutputImageType > FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput( myImage2D.originalImages());
    filter->SetNumberOfIterations(this->ui->leAD_Iteration->text().toDouble());
    filter->SetTimeStep( 0.125 );
    filter->SetConductanceParameter(this->ui->leAD_Conductance->text().toDouble());

    //convert output from float to uchar2d
    typedef itk::RescaleIntensityImageFilter<OutputImageType, InputImageType> RescaleType;
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

bool SimpleView::up_pixel_same(ImageType::Pointer seg_image, ImageType::IndexType pixelIndex, short pixel_value)
{
    signed short p_value;
    ImageType::IndexType p_index;
    p_index[0] = pixelIndex[0];
    for (int i = 1; i <= 2; i++)
    {
        p_index[1] = pixelIndex[1] - i;
        p_value = seg_image->GetPixel(p_index);
        if(abs(p_value - pixel_value) > 20)
            return false;
    }
    return true;
}

void SimpleView::slotTest()
{
    //temp for sobel
    typedef itk::Image<float, 2>          FloatImageType;
    typedef itk::SobelEdgeDetectionImageFilter <ImageType, FloatImageType>
            SobelEdgeDetectionImageFilterType;
    SobelEdgeDetectionImageFilterType::Pointer sobelFilter
            = SobelEdgeDetectionImageFilterType::New();
    sobelFilter->SetInput(myImage2D.originalImages());
    //convert output from float to uchar2d
    typedef itk::RescaleIntensityImageFilter<FloatImageType, ImageType> RescaleType;
    RescaleType::Pointer rescaler = RescaleType::New();
    rescaler->SetInput( sobelFilter->GetOutput() );
    rescaler->SetOutputMinimum( itk::NumericTraits< unsigned char >::min() );
    rescaler->SetOutputMaximum( itk::NumericTraits< unsigned char >::max() );

    myImage2D.readFromOtherOutput(rescaler->GetOutput());
    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
    this->ui->qvtkWidget_Seg->repaint();
    this->ui->qvtkWidget_Seg->show();
}

void SimpleView::slotSegmentation()
{
#if 1   //temp test for finding an edge from image center to up
    ImageType::Pointer seg_image = myImage2D.originalImages();
    int w = seg_image->GetLargestPossibleRegion().GetSize()[0];
    int h = seg_image->GetLargestPossibleRegion().GetSize()[1];
    ImageType::IndexType pixelIndex;
    signed short pixel_value;
    bool b1 = false;
//    int paint_cnt = 0;

    std::cout << "size is " << w << " x " << h << endl;
    for (int x = 200; x < w; x++)
    {
//        paint_cnt = 0;
        for (int y = h/2; y > 0; y--)
        {
            pixelIndex[0] = x;
            pixelIndex[1] = y;
            pixel_value = seg_image->GetPixel(pixelIndex);
            if(pixel_value > 100)
                b1 = true;
            else
                b1 = false;

            if(b1 && up_pixel_same(seg_image, pixelIndex, pixel_value))
            {
                for (int i = 0; i < 5; i++)
                {
                    pixelIndex[1] = pixelIndex[1] - i;
                    seg_image->SetPixel(pixelIndex, 255);
                }
                y = 0;
            }
        }
    }
//    for (int x = 0; x < w; x++)
//    {
//        pixelIndex[0] = x;
//        pixelIndex[1] = h/5;
//        seg_image->SetPixel(pixelIndex, 255);
//    }
    myImage2D.readFromOtherImage(seg_image);
    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
    this->ui->qvtkWidget_Seg->repaint();
    this->ui->qvtkWidget_Seg->show();
//    ImageType::IndexType pixelIndex;
//    signed short pixel_value;
//    for (int c = 0; c < 800;  c++)
//    {
//        pixelIndex[0] = 580;
//        pixelIndex[1] = c;
//        pixel_value = reader->GetOutput()->GetPixel(pixelIndex);
//        std::cout << "pixel before is : " << pixel_value;
//        reader->GetOutput()->SetPixel(pixelIndex, 255);
//        pixel_value = reader->GetOutput()->GetPixel(pixelIndex);
//        std::cout << "pixel after  is : " << pixel_value;
//    }
#endif
}

void SimpleView::slotExit() {
  qApp->exit();
}

void SimpleView::displayImage2(vtkImageData *image)
{
#if 1 //for picker
    vtkSmartPointer<vtkImageViewer2> imageViewer =
      vtkSmartPointer<vtkImageViewer2>::New();
    imageViewer->SetInputData(image);
    // Picker to pick pixels
    vtkSmartPointer<vtkPropPicker> propPicker =
      vtkSmartPointer<vtkPropPicker>::New();
    propPicker->PickFromListOn();

    // Give the picker a prop to pick
    vtkImageActor* imageActor = imageViewer->GetImageActor();
    propPicker->AddPickList(imageActor);

    // disable interpolation, so we can see each pixel
    imageActor->InterpolateOff();

    // Visualize
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//      vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    imageViewer->SetupInteractor(renderWindowInteractor);
//    imageViewer->SetSize(600, 600);

    vtkRenderer* renderer = imageViewer->GetRenderer();
    renderer->ResetCamera();
    renderer->GradientBackgroundOn();
    renderer->SetBackground(0.6, 0.6, 0.5);
    renderer->SetBackground2(0.3, 0.3, 0.2);

    // Annotate the image with window/level and mouse over pixel
    // information
    vtkSmartPointer<vtkCornerAnnotation> cornerAnnotation =
      vtkSmartPointer<vtkCornerAnnotation>::New();
    cornerAnnotation->SetLinearFontScaleFactor(2);
    cornerAnnotation->SetNonlinearFontScaleFactor(1);
    cornerAnnotation->SetMaximumFontSize(20);
    cornerAnnotation->SetText(0, "Off Image");
    cornerAnnotation->SetText(3, "<window>\n<level>");
    cornerAnnotation->GetTextProperty()->SetColor(1, 0, 0);

    imageViewer->GetRenderer()->AddViewProp(cornerAnnotation);

    // Callback listens to MouseMoveEvents invoked by the interactor's style
    vtkSmartPointer<vtkImageInteractionCallback> callback =
      vtkSmartPointer<vtkImageInteractionCallback>::New();
    callback->SetViewer(imageViewer);
    callback->SetAnnotation(cornerAnnotation);
    callback->SetPicker(propPicker);

    // InteractorStyleImage allows for the following controls:
    // 1) middle mouse + move = camera pan
    // 2) left mouse + move = window/level
    // 3) right mouse + move = camera zoom
    // 4) middle mouse wheel scroll = zoom
    // 5) 'r' = reset window/level
    // 6) shift + 'r' = reset camera
    vtkInteractorStyleImage* imageStyle =
      imageViewer->GetInteractorStyle();
    imageStyle->AddObserver(vtkCommand::MouseMoveEvent, callback);
    this->ui->qvtkWidget_Ori->GetInteractor()->SetInteractorStyle(imageStyle);

//    renderWindowInteractor->Initialize();
//    renderWindowInteractor->Start();
    this->ui->qvtkWidget_Ori->SetRenderWindow(imageViewer->GetRenderWindow());
//    imageViewer->SetupInteractor(this->ui->qvtkWidget_Ori->GetRenderWindow()->GetInteractor());
//    imageViewer->Render();
    this->ui->qvtkWidget_Ori->update();
#endif
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
