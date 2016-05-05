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
#include <QGraphicsView>
//ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkGDCMImageIO.h"
#include "itkBMPImageIOFactory.h"
#include "itkJPEGImageIOFactory.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkRegionGrowImageFilter.h"
//#include "itkKLMRegionGrowImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkNeighborhoodConnectedImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
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
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleImage.h>
// double click
#include <vtkCallbackCommand.h>
#include <vtkCoordinate.h>
#include <vtkRendererCollection.h>
// color point
#include <vtkImageTracerWidget.h>
#include <vtkProperty.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>

#include "vtkSmartPointer.h"
#include "C_itkSeg.h"
#include "spline.h"
#define VTK_CREATE(type, name) \
    vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

using namespace std;
//static C_fileIO myImage2D;

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

    void SetViewer2(QVTKWidget *qvtkw)
    {
        this->qvtkw = qvtkw;
    }

    void SetActor(vtkImageActor* actor)
    {
        this->Actor = actor;
    }

    virtual void Execute(vtkObject *, unsigned long vtkNotUsed(event), void *)
    {
        //    vtkRenderWindowInteractor *interactor = this->Viewer->GetRenderWindow()->GetInteractor();
        //      vtkRenderer* renderer = this->Viewer->GetRenderer();
        //      vtkImageActor* actor = this->Viewer->GetImageActor();
        //      vtkImageData* image = this->Viewer->GetInput();
        //      vtkInteractorStyle *style = vtkInteractorStyle::SafeDownCast(
        //        interactor->GetInteractorStyle());
        vtkRenderWindowInteractor *interactor = this->qvtkw->GetRenderWindow()->GetInteractor();
        vtkRenderer* renderer = this->qvtkw->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
        vtkImageActor* actor = this->Actor;//(vtkImageActor*)(renderer->GetActors()->GetLastActor());
        vtkImageData* image = actor->GetInput();
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
        image_coordinate[0] = vtkMath::Round(pos[0]);
        image_coordinate[1] = vtkMath::Round(pos[1]);
        image_coordinate[2] = 0;

        //    int axis = this->Viewer->GetSliceOrientation();
        //    switch (axis)
        //      {
        //      case vtkImageViewer2::SLICE_ORIENTATION_XZ:
        //        image_coordinate[0] = vtkMath::Round(pos[0]);
        //        image_coordinate[1] = this->Viewer->GetSlice();
        //        image_coordinate[2] = vtkMath::Round(pos[2]);
        //        break;
        //      case vtkImageViewer2::SLICE_ORIENTATION_YZ:
        //        image_coordinate[0] = this->Viewer->GetSlice();
        //        image_coordinate[1] = vtkMath::Round(pos[0]);
        //        image_coordinate[2] = vtkMath::Round(pos[1]);
        //        break;
        //      default:  // vtkImageViewer2::SLICE_ORIENTATION_XY
        //        image_coordinate[0] = vtkMath::Round(pos[0]);
        //        image_coordinate[1] = vtkMath::Round(pos[1]);
        //        image_coordinate[2] = this->Viewer->GetSlice();
        //        break;
        //      }

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
        printf("gogogogo!!!\n");
        printf("%s\n", message.c_str());
        interactor->Render();
        style->OnMouseMove();
    }

private:
    vtkImageViewer2*      Viewer;      // Pointer to the viewer
    vtkPropPicker*        Picker;      // Pointer to the picker
    vtkCornerAnnotation*  Annotation;  // Pointer to the annotation
    QVTKWidget*           qvtkw;
    vtkImageActor*        Actor;
};

class MyView:public QGraphicsView
{
public:
    MyView( SimpleView *parent=0 ) :  QGraphicsView(parent = 0)
    {
        this->Parent = parent;
        qlPixelInfo = new QLabel(this);
        qlPixelInfo->setStyleSheet("QLabel { color : red; }");
        this->PointCount = 0;
    }

    void SetImage(QImage img)
    {
        this->Qimage = img.convertToFormat(QImage::Format_RGB888);
    }

    void Display()
    {
        setScene( new QGraphicsScene( this ) );
        this->scene()->addPixmap(QPixmap::fromImage(this->Qimage)); // this->image is a QImage
        this->show();
    }

    void setViewMode(ViewMode mode)
    {
        this->Mode = mode;
    }

public slots:
    void mousePressEvent( QMouseEvent* event )
    {
        QPointF pos =  mapToScene( event->pos() );
        int x = ( int )pos.x();
        int y = ( int )pos.y();
        if (this->Mode == None || this->Mode == CalculateDistance)
        {
            ClickPointX[PointCount] = x;
            ClickPointY[PointCount] = y;
            PointCount++;
            QRgb rgb = this->Qimage.pixel( x, y );
            drawGreenDot(x,y);
            if (this->Mode == CalculateDistance && PointCount % 2 == 0)
            {
                QPainter pt(&Qimage);
                pt.setPen(Qt::green);
                pt.drawLine(ClickPointX[PointCount-1],ClickPointY[PointCount-1],
                        ClickPointX[PointCount-2],ClickPointY[PointCount-2] );
                pt.end();
                QLabel *qlDistance = new QLabel(this);
                qlDistance->setStyleSheet("QLabel { color : red; }");
                qlDistance->setGeometry(0, 0+(PointCount/2)*30, 300, 30);
                QString line_info;
                double distance = sqrt(pow(ClickPointX[PointCount-1] - ClickPointX[PointCount-2], 2) +
                        pow(ClickPointY[PointCount-1] - ClickPointY[PointCount-2], 2));
                line_info.sprintf("(%d,%d) to (%d,%d) is %lf",
                                  ClickPointX[PointCount-1],ClickPointY[PointCount-1],
                        ClickPointX[PointCount-2],ClickPointY[PointCount-2], distance);
                qlDistance->setText(line_info);
                qlDistance->show();
            }

            this->Display();
            QString info;
            info.sprintf("(%d,%d)=(%d,%d,%d)", x, y, qRed(rgb), qGreen(rgb), qBlue(rgb));

            qlPixelInfo->hide();
            qlPixelInfo->setGeometry(0, 0, 300, 30);
            qlPixelInfo->setText(info);
            qlPixelInfo->show();
        }
        else if (this->Mode == RegionGrowing)
        {
            this->Qimage.save("region_growing.bmp");
            region_growing("region_growing.bmp", x, y);
            this->SetImage(QImage("region_growing.bmp"));
            this->Display();
        }
    }

private slots:
    QImage region_growing(QString image_file, int seed_x, int seed_y)
    {
        ITKImageReader = ReaderType::New();
        ITKImageReader->SetFileName(image_file.toLatin1().data());
        typedef itk::NeighborhoodConnectedImageFilter<ImageType, ImageType> RegionGrowImageFilterType;
        RegionGrowImageFilterType::Pointer regionGrow = RegionGrowImageFilterType::New();
        float lower = 0.0;
        float upper = 50.0;
        regionGrow->SetInput(ITKImageReader->GetOutput());
        regionGrow->SetLower(lower);
        regionGrow->SetUpper(upper);

        regionGrow->SetReplaceValue(255);

        // Seed 1: (25, 35)
        ImageType::IndexType seed1;
        seed1[0] = seed_x;
        seed1[1] = seed_y;
        regionGrow->SetSeed(seed1);
        regionGrow->Update();

        //opening
        typedef itk::BinaryBallStructuringElement<ImageType::PixelType, ImageType::ImageDimension>
                    StructuringElementType;
        StructuringElementType structuringElement_opening;
        structuringElement_opening.SetRadius(5);
        structuringElement_opening.CreateStructuringElement();
        typedef itk::BinaryMorphologicalOpeningImageFilter <ImageType, ImageType, StructuringElementType>
                BinaryMorphologicalOpeningImageFilterType;
        BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter
                = BinaryMorphologicalOpeningImageFilterType::New();
        openingFilter->SetInput(regionGrow->GetOutput());
        openingFilter->SetKernel(structuringElement_opening);
        openingFilter->Update();
        C_fileIO tmp_image;
        tmp_image.readFromOtherOutput(regionGrow->GetOutput());
        tmp_image.writeImageToFile("region_growing.bmp");
    }

    void drawGreenDot(int x, int y)
    {
        for (int r = x-2; r < x+3; r++)
        {
            for (int c = y-2; c < y+3; c++)
            {
                this->Qimage.setPixel(r, c, qRgb(0,255,0));
            }
        }
    }

private:
    QImage Qimage;
    ReaderType::Pointer ITKImageReader;
    WriterType::Pointer ITKImageWriter;
    QLabel *qlPixelInfo;
    SimpleView *Parent;
    int PointCount;
    //to recording 6 clicked point
    int ClickPointX[6], ClickPointY[6];
    ViewMode Mode;
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
    connect(this->ui->btnTest, SIGNAL (released()),this, SLOT (slotTest()));
    connect(this->ui->btnSpline, SIGNAL (released()),this, SLOT (slotSpline()));
    connect(this->ui->btnMerge, SIGNAL (released()),this, SLOT (slotMerge()));
    connect(this->ui->btnSobel, SIGNAL (released()),this, SLOT (slotSobel()));
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
    itk::JPEGImageIOFactory::RegisterOneFactory();
    InputFile = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());

#if 0 //for testing use
#else
    // get DICOM info
    typedef itk::GDCMImageIO ImageIOType;
    ImageIOType::Pointer gdcmIO = ImageIOType::New();
    gdcmIO->UseCompressionOn();
    gdcmIO->SetCompressionType(itk::GDCMImageIO::JPEG);
    gdcmIO->UseCompressionOff();

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(InputFile.toLatin1().data());
    //    reader->SetImageIO(gdcmIO);
    reader->Update();

    /***** get DICOM tag value *****/
    QString dicomTagValue;
    int height, width;
    //int maxPixelValue, minPixelValue;
    int winCenter, winWidth;
    double pixelSpacing_row, pixelSpacing_col;

    height = (getDICOMtagValue(gdcmIO, "0028|0010")).toInt(); //height
    width = (getDICOMtagValue(gdcmIO, "0028|0011")).toInt(); //width
    //maxPixelValue = stoi(getDICOMtagValue(gdcmIO, "0028|0107")); //not found
    //minPixelValue = stoi(getDICOMtagValue(gdcmIO, "0028|0106")); //not found
    winCenter = (getDICOMtagValue(gdcmIO, "0028|1050")).toInt(); //window center
    winWidth = (getDICOMtagValue(gdcmIO, "0028|1051")).toInt(); //window width
    // pixelSpacing
    dicomTagValue = getDICOMtagValue(gdcmIO, "0018|1164") + ".";
    //    pixelSpacing_row = (dicomTagValue.substr(0, dicomTagValue.find("\\")));
    //    pixelSpacing_col = (dicomTagValue.substr(dicomTagValue.find("\\")+1, dicomTagValue.find(".")-dicomTagValue.find("\\")-1));

#if 0 //Normalize
    /***** Normalize DICOM pixel value (0~255) *****/
    ImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;

    ImageType::SizeType size;
    size[0] = width;
    size[1] = height;

    // Pixel data is allocated
    ImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    ImageType::Pointer image = reader->GetOutput();
    image->SetRegions(region);
    image->Allocate();

    int outMin=0, outMax=255; // Normalize to 0~255
    ImageType::IndexType pixelIndex;
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            pixelIndex[0] = i;
            pixelIndex[1] = j;
            ImageType::PixelType pixelVal;
            pixelVal = image->GetPixel(pixelIndex); // origin value
            // Normalization
            if (pixelVal <= winCenter - 0.5 - (winWidth-1)/2) {
                pixelVal = outMin;
            }
            else if (pixelVal > winCenter - 0.5 + (winWidth-1)/2) {
                pixelVal = outMax;
            }
            else {
                pixelVal = ((pixelVal - (winCenter - 0.5)) / (winWidth-1) + 0.5) * (outMax - outMin )+ outMin;
            }
            image->SetPixel(pixelIndex, pixelVal); // normalized value
        }
    }
#endif
    myImage2D.readFromOtherOutput(reader->GetOutput());
    myImage2D.writeImageToFile("dicom_test.jpg");
#endif
    /***** 變數設定 *****/
    // image is grayscale
    //    typedef itk::Image<unsigned short,2> ImageType; //要用unsigned short以上型態，dicom pixel value才不會爆掉
    //typedef itk::Image<signed int,2> ImageType;
    // image is RGB

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

    QImage img(InputFile);
    this->displayMyView(img, RegionGrowing);
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
    myImage2D.writeImageToFile("KneeOut.bmp");
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

void SimpleView::slotSobel()
{
    myImage2D.readFiletoImages("region_growing.bmp");
    //sobel
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
    myImage2D.writeImageToFile("KneeOut_sobel.bmp");

    //convert sobel to red line
    int row, col;
    QImage inImg("KneeOut_sobel.bmp");
    QImage outImg = QImage(InputFile.toLatin1().data());
    QImage outImgC = outImg.convertToFormat(QImage::Format_RGB888);
    for (col = 0; col<outImgC.width(); col++)
    {
        for(row = 0; row < outImgC.height(); row++)
        {
            QRgb rgb = inImg.pixel(col, row);
            if(qRed(rgb) > this->ui->leTestValue->text().toInt())
            {
                outImgC.setPixel(col, row, qRgb(255, 0, 0));
            }
            else
            {
//                outImgC.setPixel(col, row, qRgb(0, 0, 0));
            }
        }
    }
    this->displayMyView(outImgC, None);
    outImgC.save("KneeOut_sobel_red.bmp");
#if 0//test code
    // Create an image data
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->DeepCopy(myImage2D.vtkImage());
    int* dims = myImage2D.vtkImage()->GetDimensions();
    imageData->AllocateScalars(VTK_DOUBLE,1);
    imageData->SetDimensions(dims);
    // int dims[3]; // can't do this

    std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;

    std::cout << "Number of points: " << imageData->GetNumberOfPoints() << std::endl;
    std::cout << "Number of cells: " << imageData->GetNumberOfCells() << std::endl;
    for (int z = 0; z < dims[2]; z++)
    {
        for (int y = 0; y < dims[1]; y++)
        {
            for (int x = 0; x < dims[0]; x++)
            {
                int* pixel = static_cast<int*>(myImage2D.vtkImage()->GetScalarPointer(x,y,z));
                //                if (pixel[0] != 0)
                //                {
                //                    imageData->SetScalarComponentFromDouble(x,y,z,0, 255);
                //                    imageData->SetScalarComponentFromDouble(x,y,z,1, 0);
                //                    imageData->SetScalarComponentFromDouble(x,y,z,2, 0);
                //                }
                //                if (myImage2D.vtkImage()->getscGetScalarComponentAsDouble(x,y,z,0) != 0.0)
                //                {
                //                    imageData->SetScalarComponentFromDouble(x,y,z,0, 255);
                //                    imageData->SetScalarComponentFromDouble(x,y,z,1, 0);
                //                    imageData->SetScalarComponentFromDouble(x,y,z,2, 0);
                //                }
                //                else
                //                {
                //                    imageData->SetScalarComponentFromDouble(x,y,z,0, 0);
                //                    imageData->SetScalarComponentFromDouble(x,y,z,1, 0);
                //                    imageData->SetScalarComponentFromDouble(x,y,z,2, 0);
                //                }
            }
        }
    }
    displayImage(imageData);

#endif
}

void SimpleView::slotMerge()
{
    myImage2D.writeImageToFile("KneeOut_merge.bmp");
    QImage inImg("KneeOut_merge.bmp");
    int Original_Image_Height = inImg.height();
    int Original_Image_Width = inImg.width();
    float ClassNumber;                                 //the number of classes used by MRF
    float SelectingClass;                              //how many classes would be merged?  1+2? 1+2+3? 1+2+3+4?...
    int ThresholdingValue;


    ClassNumber = (float)(this->ui->leIterationNum->text().toInt());
    SelectingClass = (float)(this->ui->leMergedCls->text().toInt());
    ThresholdingValue = (int)(((255.0 / (ClassNumber - 1.0)) * (SelectingClass - 1.0)) + 2.0);
    for(int j = 0;j < Original_Image_Height;j++)
    {
        for(int i = 0;i < Original_Image_Width;i++)
        {
            QRgb rgb = inImg.pixel(i, j);
            if(qRed(rgb) <= ThresholdingValue)
                inImg.setPixel(i, j, 0);
            else
                inImg.setPixel(i, j, 255);
        }
    }
    this->displayMyView(inImg, RegionGrowing);
    inImg.save("KneeOut_merge.bmp");
}

void SimpleView::region_growing(int seed_x, int seed_y)
{
    typedef itk::ConnectedThresholdImageFilter<ImageType, ImageType> RegionGrowImageFilterType;
    RegionGrowImageFilterType::Pointer regionGrow = RegionGrowImageFilterType::New();
    float lower = 0.0;
    float upper = 50.0;
    regionGrow->SetInput(myImage2D.originalImages());
    regionGrow->SetLower(lower);
    regionGrow->SetUpper(upper);

    regionGrow->SetReplaceValue(100);

    // Seed 1: (25, 35)
    ImageType::IndexType seed1;
    seed1[0] = 0;
    seed1[1] = 0;
    ImageType::IndexType seed2;
    seed2[0] = 0;
    seed2[1] = 175;
    regionGrow->SetSeed(seed1);
    regionGrow->SetSeed(seed2);
    regionGrow->Update();

    myImage2D.readFromOtherOutput(regionGrow->GetOutput());
}

void SimpleView::slotTest()
{
    //merge other fragments to the cartilage
    ImageType::IndexType index;
    int width = myImage2D.originalImages()->GetLargestPossibleRegion().GetSize()[0];
    int height = myImage2D.originalImages()->GetLargestPossibleRegion().GetSize()[1];
    cout << "w is: "<<width<<endl;
    cout << "h is: "<<height<<endl;
    bool b_find_region = false;
    for(int y = 0; y < height; y++)
    {
        for(int x = 0; x < width; x++)
        {
            index[0] = x;
            index[1] = y;
            if (myImage2D.originalImages()->GetPixel(index) == 0)
            {
                region_growing(x, y);
                b_find_region = true;
                cout << "x1 is: "<<x<<endl;
                cout << "y1 is: "<<y<<endl;
                break;
            }
        }
        if(b_find_region) break;
    }
//    b_find_region = false;
//    for(int y = height - 1; y >=0; y--)
//    {
//        for(int x = 0; x < width; x++)
//        {
//            index[0] = x;
//            index[1] = y;
//            if (myImage2D.originalImages()->GetPixel(index) == 0)
//            {
//                region_growing(x, y);
//                b_find_region = true;
//                cout << "x2 is: "<<x<<endl;
//                cout << "y2 is: "<<y<<endl;
//                break;
//            }
//        }
//        if(b_find_region) break;
//    }
    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
    this->ui->qvtkWidget_Seg->repaint();
    this->ui->qvtkWidget_Seg->show();

}

void SimpleView::Opening()
{
#if 1 //openinig & closing test
    typedef itk::BinaryBallStructuringElement<ImageType::PixelType, ImageType::ImageDimension>
                StructuringElementType;

    //closing
    StructuringElementType structuringElement_closing;
    structuringElement_closing.SetRadius(5);
    structuringElement_closing.CreateStructuringElement();
    typedef itk::BinaryMorphologicalClosingImageFilter <ImageType, ImageType, StructuringElementType>
            BinaryMorphologicalClosingImageFilterType;
    BinaryMorphologicalClosingImageFilterType::Pointer closingFilter
            = BinaryMorphologicalClosingImageFilterType::New();
    closingFilter->SetInput(myImage2D.originalImages());
    closingFilter->SetKernel(structuringElement_closing);
    closingFilter->Update();
    myImage2D.readFromOtherOutput(closingFilter->GetOutput());

    //opening
//    StructuringElementType structuringElement_opening;
//    structuringElement_opening.SetRadius(5);
//    structuringElement_opening.CreateStructuringElement();
//    typedef itk::BinaryMorphologicalOpeningImageFilter <ImageType, ImageType, StructuringElementType>
//            BinaryMorphologicalOpeningImageFilterType;
//    BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter
//            = BinaryMorphologicalOpeningImageFilterType::New();
//    openingFilter->SetInput(myImage2D.originalImages());
//    openingFilter->SetKernel(structuringElement_opening);
//    openingFilter->Update();
//    myImage2D.readFromOtherOutput(openingFilter->GetOutput());

    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
    this->ui->qvtkWidget_Seg->repaint();
    this->ui->qvtkWidget_Seg->show();
#endif
}

void drawColorDot(QImage* img, int rgb, int x, int y)
{
    for (int r = x-2; r < x+3; r++)
    {
        for (int c = y-2; c < y+3; c++)
        {
            img->setPixel(r, c, rgb);
        }
    }
}


void SimpleView::slotSpline()
{
#if 1 //test
    QImage inImg = QImage("KneeOut_sobel_red.bmp");
    QImage outImg = QImage(InputFile.toLatin1().data());
    QImage outImgC = outImg.convertToFormat(QImage::Format_RGB888);
    std::vector<double> x1, y1, x2, y2;
    tk::spline s1, s2;
    int last_red_y1 = -1, last_red_x1 = -1;
    int last_red_y2 = -1, last_red_x2 = -1;
    float slope = 0.0;
    for(int x = 0; x < inImg.width(); x += 15)
    {
        for(int y = 0; y < inImg.height(); y++)
        {
            QRgb rgb = inImg.pixel(x, y);
            if(qRed(rgb) >=255)
            {
                slope = (float)abs(y-last_red_y1) / (float)abs(x-last_red_x1);
                /*must sample first point*/
                cout<<"(x,y) = ("<<x<<","<<y<<")";
                cout<<"(last_red_x1, last_red_y1) = ("<<last_red_x1<<","<<last_red_y1<<")";
                cout<<"slope = "<<slope<<endl;
                if(last_red_y1 == -1 || slope < 0.6)
                {
                    x1.push_back(x);
                    y1.push_back(y);
                    drawColorDot(&inImg,qRgb(0,255,0), x, y);
                    last_red_x1 = x;
                    last_red_y1 = y;
//                    cout<<"in"<<endl;
//                    cout<<"(x,y) = ("<<x<<","<<y<<")";
//                    cout<<"slope = "<<slope<<endl;
                }
                break;
            }
        }
        for(int y = inImg.height()-1; y >= 0; y--)
        {
            QRgb rgb = inImg.pixel(x, y);
            if(qRed(rgb) >=255)
            {
                slope = (float)abs(y-last_red_y2) / (float)abs(x-last_red_x2);
                if(last_red_y2 == -1 || slope < 0.6)
                {
                    x2.push_back(x);
                    y2.push_back(y);
                    drawColorDot(&inImg,qRgb(0,0,255), x, y);
                    last_red_x2 = x;
                    last_red_y2 = y;
                }
//                cout<<"(x,y) = ("<<x<<","<<y<<")"<<endl;
                break;
            }
        }
    }
    inImg.save("spline_test.bmp");
    s1.set_points(x1,y1);
    s2.set_points(x2,y2);
    for (int col = 0; col<outImgC.width(); col++)
    {
        outImgC.setPixel(col, (int)s1((double)col), qRgb(255, 0, 0));
        outImgC.setPixel(col, (int)s2((double)col), qRgb(255, 0, 0));
    }
    this->displayMyView(outImgC, None);
#endif
#if 0//test for addimage filter
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName("/home/ryderlin/Documents/001.bmp");
    reader->Update();

    ReaderType::Pointer reader1 = ReaderType::New();
    reader1->SetFileName("KneeOut_sobel_red.bmp");
    reader1->Update();
    typedef itk::AddImageFilter <ImageType, ImageType >
            AddImageFilterType;

    AddImageFilterType::Pointer addFilter
            = AddImageFilterType::New ();
    addFilter->SetInput1(reader->GetOutput());
    addFilter->SetInput2(reader1->GetOutput());
    addFilter->Update();
    myImage2D.readFromOtherImage(addFilter->GetOutput());
    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
    this->ui->qvtkWidget_Seg->repaint();
    this->ui->qvtkWidget_Seg->show();
#endif
#if 0   //temp test for finding an edge from image center to up
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

//temp test for pick pixel on QVTKWidget
void SimpleView::displayImage2(vtkImageData *image)
{
    // Picker to pick pixels
    vtkSmartPointer<vtkPropPicker> propPicker =
            vtkSmartPointer<vtkPropPicker>::New();
    propPicker->PickFromListOn();
    // Create image actor
    vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(800, 600);

    // set actor properties
    actor->SetInputData(image);
    // Give the picker a prop to pick
    propPicker->AddPickList(actor);

    // disable interpolation, so we can see each pixel
    actor->InterpolateOff();

    renderer->AddActor(actor);
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

    renderer->AddViewProp(cornerAnnotation);
    // Callback listens to MouseMoveEvents invoked by the interactor's style
    vtkSmartPointer<vtkImageInteractionCallback> callback =
            vtkSmartPointer<vtkImageInteractionCallback>::New();
    //    callback->SetViewer(imageViewer);
    callback->SetAnnotation(cornerAnnotation);
    callback->SetPicker(propPicker);
    callback->SetViewer2(this->ui->qvtkWidget_Ori);
    callback->SetActor(actor);

    this->ui->qvtkWidget_Ori->GetRenderWindow()->GetInteractor()->Initialize();
    this->ui->qvtkWidget_Ori->GetRenderWindow()->GetInteractor()->Start();
    this->ui->qvtkWidget_Ori->SetRenderWindow(renderWindow);

    // window interactor style for display images
    vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
    // set interactor style to the qvtkWidget Interactor
    this->ui->qvtkWidget_Ori->GetInteractor()->SetInteractorStyle(style);
    this->ui->qvtkWidget_Ori->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::MouseMoveEvent, callback);
    this->ui->qvtkWidget_Ori->update();
}

//using QImage to show original image
void SimpleView::displayMyView(QImage img, ViewMode view_mode)
{
    MyView *mv = new(std::nothrow) MyView(this);
    mv->setGeometry(100, 200, 500, 500);
    mv->SetImage(img);
    mv->setViewMode(view_mode);
    mv->Display();
}

void SimpleView::updatePixInfo(QString pix_info)
{

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
#if 1   //test from chris, get mouse event and put a pixel on it
    // get double click events
    vtkCallbackCommand *callback = vtkCallbackCommand::New();
    callback->SetCallback(SimpleView::handle_double_click);
    callback->SetClientData(this);
    this->ui->qvtkWidget_Ori->GetInteractor()->AddObserver(vtkCommand::LeftButtonPressEvent, callback, 1.0);
#endif
}

void SimpleView::handle_double_click(vtkObject* obj, unsigned long event, void* ClientData, void* CallData) {
    SimpleView* self = reinterpret_cast<SimpleView*>(ClientData);
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(obj);

    int currPos[2];
    iren->GetEventPosition(currPos);
    cout << currPos[0] <<", " << currPos[1] << endl;

    int x, y;
    x = currPos[0];
    y = currPos[1];
    vtkSmartPointer<vtkCoordinate> coordinate = vtkSmartPointer<vtkCoordinate>::New();
    coordinate->SetCoordinateSystemToDisplay();
    coordinate->SetValue(x,y,0);

    // This doesn't produce the right value if the sphere is zoomed in???
    double *world = coordinate->GetComputedWorldValue(iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
    cout << "World coordinate: " << world[0] << ", " << world[1] << ", " << world[2] << endl;

    // Draw colored point
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    //points->InsertNextPoint (0.0, 0.0, 0.0);
    //points->InsertNextPoint (x, y, 0.0);
    points->InsertNextPoint (world[0], world[1], world[2]);

    vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
    pointsPolydata->SetPoints(points);

    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilter->SetInputData(pointsPolydata);
    vertexFilter->Update();

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->ShallowCopy(vertexFilter->GetOutput());

    // Setup colors
    //unsigned char red[3] = {255, 0, 0};
    unsigned char green[3] = {0, 255, 0};
    //unsigned char blue[3] = {0, 0, 255};

    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName ("Colors");
    //colors->InsertNextTupleValue(red);
    colors->InsertNextTupleValue(green);

    polydata->GetPointData()->SetScalars(colors);

    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polydata);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(5);

    // store camera original value
    double *tmp_p = self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->GetPosition();
    double *tmp_fp = self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->GetFocalPoint();
    double p[3] = {tmp_p[0], tmp_p[1], tmp_p[2]};
    double fp[3] = {tmp_fp[0], tmp_fp[1], tmp_fp[2]};
    // draw color point
    self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actor);
    // reset camera
    self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->ResetCamera();
    self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->Zoom(1.5);
    // restore camera with original value
    self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->SetPosition(p);
    self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->SetFocalPoint(fp);
    self->ui->qvtkWidget_Ori->GetRenderWindow()->Render();
    self->ui->qvtkWidget_Ori->update();
}

QString SimpleView::getDICOMtagValue(itk::GDCMImageIO::Pointer gdcmIO, QString tagkey) {
    string labelId;
    if( itk::GDCMImageIO::GetLabelFromTag( tagkey.toStdString(), labelId ) ) {
        string value;
        cout << labelId << " (" << tagkey.toStdString() << "): ";
        if( gdcmIO->GetValueFromTag(tagkey.toStdString(), value) ) {
            cout << value << endl;
        }
        else {
            cout << "(No Value Found in File)" << endl;
        }
        QString ret_value = QString::fromStdString(value);
        return ret_value;
    }
    else {
        cerr << "Trying to access inexistant DICOM tag." << endl;
        return "";
    }
}
