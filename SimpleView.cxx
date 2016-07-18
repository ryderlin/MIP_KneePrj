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
#include "global_typedef.h"
//QT
#include<QFileDialog>
#include<QDir>
#include <QGraphicsView>
#include "QVTKWidget.h"
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
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkRegionGrowImageFilter.h"
//#include "itkKLMRegionGrowImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkNeighborhoodConnectedImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
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


//static utility
static QImage region_growing(QString image_file, QString out_file, int seed_x, int seed_y, int replaced_pixel)
{
    ReaderType::Pointer ITKImageReader;
    ITKImageReader = ReaderType::New();
    ITKImageReader->SetFileName(image_file.toLatin1().data());
    typedef itk::NeighborhoodConnectedImageFilter<ImageType, ImageType> RegionGrowImageFilterType;
    RegionGrowImageFilterType::Pointer regionGrow = RegionGrowImageFilterType::New();
    float lower = 0.0;
    float upper = 50.0;
    regionGrow->SetInput(ITKImageReader->GetOutput());
    regionGrow->SetLower(lower);
    regionGrow->SetUpper(upper);

    regionGrow->SetReplaceValue(replaced_pixel);

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
    tmp_image.writeImageToFile(out_file.toLatin1().data());
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
            drawColorDot(&Qimage, qRgb(0,255,0), x, y);
//            drawGreenDot(x,y);
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
            this->Qimage.save(FILE_REGION_GROWING);
            region_growing(FILE_REGION_GROWING, FILE_REGION_GROWING, x, y, 255);
            this->SetImage(QImage(FILE_REGION_GROWING));
            this->Display();
        }
    }

private slots:
//    QImage region_growing(QString image_file, int seed_x, int seed_y)
//    {
//        ITKImageReader = ReaderType::New();
//        ITKImageReader->SetFileName(image_file.toLatin1().data());
//        typedef itk::NeighborhoodConnectedImageFilter<ImageType, ImageType> RegionGrowImageFilterType;
//        RegionGrowImageFilterType::Pointer regionGrow = RegionGrowImageFilterType::New();
//        float lower = 0.0;
//        float upper = 50.0;
//        regionGrow->SetInput(ITKImageReader->GetOutput());
//        regionGrow->SetLower(lower);
//        regionGrow->SetUpper(upper);

//        regionGrow->SetReplaceValue(255);

//        // Seed 1: (25, 35)
//        ImageType::IndexType seed1;
//        seed1[0] = seed_x;
//        seed1[1] = seed_y;
//        regionGrow->SetSeed(seed1);
//        regionGrow->Update();

//        //opening
//        typedef itk::BinaryBallStructuringElement<ImageType::PixelType, ImageType::ImageDimension>
//                    StructuringElementType;
//        StructuringElementType structuringElement_opening;
//        structuringElement_opening.SetRadius(5);
//        structuringElement_opening.CreateStructuringElement();
//        typedef itk::BinaryMorphologicalOpeningImageFilter <ImageType, ImageType, StructuringElementType>
//                BinaryMorphologicalOpeningImageFilterType;
//        BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter
//                = BinaryMorphologicalOpeningImageFilterType::New();
//        openingFilter->SetInput(regionGrow->GetOutput());
//        openingFilter->SetKernel(structuringElement_opening);
//        openingFilter->Update();
//        C_fileIO tmp_image;
//        tmp_image.readFromOtherOutput(regionGrow->GetOutput());
//        tmp_image.writeImageToFile(FILE_REGION_GROWING);
//    }

//    void drawGreenDot(int x, int y)
//    {
//        for (int r = x-2; r < x+3; r++)
//        {
//            for (int c = y-2; c < y+3; c++)
//            {
//                this->Qimage.setPixel(r, c, qRgb(0,255,0));
//            }
//        }
//    }

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

    //create output file dir
    if(!QDir(OUT_FILE_DIR).exists()) QDir().mkdir(OUT_FILE_DIR);

    // Set up action signals and slots
    connect(this->ui->actionOpenFile, SIGNAL(triggered()), this, SLOT(slotOpenFile()));
    connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));
    connect(this->ui->btnRunMrf, SIGNAL (released()),this, SLOT (slotRunMrf()));
    connect(this->ui->btnRunSig, SIGNAL (released()),this, SLOT (slotRunSig()));
    connect(this->ui->btnPreProcessMrf, SIGNAL (released()),this, SLOT (slotPreProcessMrf16()));
    connect(this->ui->btnRunAD, SIGNAL (released()),this, SLOT (slotRunAD()));
    connect(this->ui->btnReset, SIGNAL (released()),this, SLOT (slotReset()));
    connect(this->ui->btnWriteFile, SIGNAL (released()),this, SLOT (slotWriteFile()));
    connect(this->ui->btnTest, SIGNAL (released()),this, SLOT (slotTest()));
    connect(this->ui->btnSpline, SIGNAL (released()),this, SLOT (slotSpline()));
    connect(this->ui->btnMerge, SIGNAL (released()),this, SLOT (slotMerge()));
    connect(this->ui->btnSobel, SIGNAL (released()),this, SLOT (slotSobel()));
    connect(this->ui->btnDeFrangment, SIGNAL (released()),this, SLOT (slotDeFragment()));
    connect(this->ui->btnSmooth, SIGNAL (released()),this, SLOT (slotSmoothEdge()));
    connect(this->ui->btnOpenDocImg, SIGNAL (released()),this, SLOT (slotOpenDocImg()));
    connect(this->ui->btnAutoCompare, SIGNAL (released()),this, SLOT (slotAutoCompare()));

    tblCmp = new QTableWidget(3, 3, this);
    tblCmp->setGeometry(20, 600, 1000, 150);
    tblCmp->setWindowTitle("compare result");
//    tblCmp->resize(350, 200);  //设置表格
    QStringList h_header, v_header;
    h_header<<"Center Thickness"<<"Left Thickness"<<"Right Thickness";
    v_header<<"Computer"<<"Doctor"<<"error %";
    tblCmp->setHorizontalHeaderLabels(h_header);
    tblCmp->setVerticalHeaderLabels(v_header);
    tblCmp->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
//    tblCmp->setItem(0,0,new QTableWidgetItem("Jan"));
//    tblCmp->setItem(1,0,new QTableWidgetItem("Feb"));
//    tblCmp->setItem(2,0,new QTableWidgetItem("Mar"));

//    tblCmp->setItem(0,1,new QTableWidgetItem("Jan's month"));
//    tblCmp->setItem(1,1,new QTableWidgetItem("Feb's month"));
//    tblCmp->setItem(2,1,new QTableWidgetItem("Mar's month"));
    tblCmp->show();
    connect(this->ui->btnSelect, SIGNAL (released()),this, SLOT (slotSelectCmpResult()));
//    this->ui->qvtkWidget_Ori->repaint();
//    this->ui->qvtkWidget_Seg->hide();
}

SimpleView::~SimpleView()
{
    // The smart pointers should clean up for up
    delete ui;
}

void SimpleView::slotSelectCmpResult()
{
    QAbstractItemModel * model = tblCmp->model();
    QItemSelectionModel * selection = tblCmp->selectionModel();
    QModelIndexList indexes = selection->selectedIndexes();
    QString selected_text;
    // You need a pair of indexes to find the row changes
    QModelIndex previous = indexes.first();
//    indexes.removeFirst();
    foreach(QModelIndex current, indexes)
    {
        QVariant data = model->data(current);
        QString text = data.toString();
        // If you are at the start of the row the row number of the previous index
        // isn't the same.  Text is followed by a row separator, which is a newline.
        if (current.row() != previous.row())
        {
            selected_text.append('\n');
        }
        // Otherwise it's the same row, so append a column separator, which is a tab.
        else if (current != previous)
        {
            selected_text.append('\t');
        }
        // At this point `text` contains the text in one cell
        selected_text.append(text);
        previous = current;
    }
    QClipboard *clipboard = QApplication::clipboard();
    clipboard->setText(selected_text);
}

void SimpleView::slotAutoCompare()
{
    itk::BMPImageIOFactory::RegisterOneFactory();
    QString auto_compare_path = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());
    QFileInfo fi(auto_compare_path);
#if 1 //test
    QDir dir(fi.absolutePath());
    QStringList nameFilter;
    nameFilter << "*.bmp";
    QFileInfoList list = dir.entryInfoList(nameFilter, QDir::Files, QDir::Name);
    int count = 0;
    qDebug() << "start auto compare "<<endl;
    foreach(QFileInfo finfo, list) {
        count++;
        InputFile = finfo.absoluteFilePath();
        qDebug() << InputFile;
        QImage img(InputFile);
        ImgW = img.width();
        ImgH = img.height();
        cout << "W x H : " << ImgW <<","<<ImgH<<endl;

        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(InputFile.toLatin1().data());
        reader->Update();
        myImage2D.readFromOtherOutput(reader->GetOutput());

        slotPreProcessMrf16();
        AutoOpenDocImage();
    }
    qDebug() << endl<<"count = "<<count;
#else
    QDirIterator it(fi.absolutePath(), QStringList() << "*.bmp", QDir::Files, QDirIterator::NoIteratorFlags);
    qDebug() << "start auto compare "<<endl;
    int count = 0;
    while (it.hasNext()) {
        count++;
        InputFile = it.next();
        qDebug() << InputFile;
        QImage img(InputFile);
        ImgW = img.width();
        ImgH = img.height();
        cout << "W x H : " << ImgW <<","<<ImgH<<endl;

        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(InputFile.toLatin1().data());
        reader->Update();
        myImage2D.readFromOtherOutput(reader->GetOutput());

        slotPreProcessMrf16();
        AutoOpenDocImage();
    }
    qDebug() << endl<<"count = "<<count;
#endif
}

// Action to be taken upon file open
void SimpleView::slotOpenFile()
{
    itk::BMPImageIOFactory::RegisterOneFactory();
//    itk::JPEGImageIOFactory::RegisterOneFactory();
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
//    slotRunAD();
//    myImage2D.writeImageToFile("dicom_test.jpg");
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
    ImgW = img.width();
    ImgH = img.height();
    cout << "W x H : " << ImgW <<","<<ImgH<<endl;
//    this->displayMyView(img, None);
    showComputerSegImage(img);
    this->displayImage(vtkimage);

    reader = NULL;
    connector = NULL;
    flipYFilter = NULL;
    vtkimage = NULL;
}

bool SimpleView::is5PixelDot(QImage in_img, int x, int y, int color)
{
    for(int i = x; i < x + 5; i++)
    {
        for (int j = y; j < y + 5; j++)
        {
            if (in_img.pixel(i,j) != color)
                return false;
        }
    }
    return true;
}

bool SimpleView::is5PixelDot_RightThickness(QImage in_img, int x, int y, int color)
{
    for(int i = x; i < x + 5; i++)
    {
        for (int j = y; j < y + 5; j++)
        {
            if (qBlue(in_img.pixel(i,j)) != qBlue(color))
                return false;
        }
    }
    return true;
}


void SimpleView::slotOpenDocImg()
{
    DoctorInputFile = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());
    showDoctorSegImage(QImage(DoctorInputFile));
    //load back the image that processed by doctor
    QImage info_img(DoctorInputFile.remove("_out.bmp") + "_out_load.bmp");
    //clear previous record
    DCenterX.clear(); DCenterY.clear(); DLeftX.clear(); DLeftY.clear(); DRightX.clear(); DRightY.clear();
    for (int x = 0; x < info_img.width(); x++)
    {
        for (int y = 0; y < info_img.height(); y++)
        {
            //get center point
            if(is5PixelDot(info_img, x, y, COLOR_CENTER_DOT_1))
            {
                qDebug() << "center point = (" << x+2 << "," << y+2 << ")" << endl;
                DCenterX.push_back(x+2);
                DCenterY.push_back(y+2);
            }
            //get left point
            if(is5PixelDot(info_img, x, y, COLOR_LEFT_DOT_1))
            {
                qDebug() << "left point = (" << x+2 << "," << y+2 << ")" << endl;
                DLeftX.push_back(x+2);
                DLeftY.push_back(y+2);
            }
            //get right point
            if(is5PixelDot_RightThickness(info_img, x, y, COLOR_RIGHT_DOT_1))
            {
                qDebug() << "right point = (" << x+2 << "," << y+2 << ")" << endl;
                DRightX.push_back(x+2);
                DRightY.push_back(y+2);
            }
        }
    }
    //show information of doctor image
    QString distance_info;
    if(DCenterX.size() >= 2)
    {
        distance_info = getDistanceInfo(DCenterX[0], DCenterY[0], DCenterX[1], DCenterY[1], &DCenterThickness);
        ui->lbDCenterThick->setText(distance_info);
        tblCmp->setItem(1,0,new QTableWidgetItem(distance_info));
    }
    if(DLeftX.size() >= 2)
    {
        distance_info = getDistanceInfo(DLeftX[0], DLeftY[0], DLeftX[1], DLeftY[1], &DLeftThickness);
        ui->lbDLeftThick->setText(distance_info);
        tblCmp->setItem(1,1,new QTableWidgetItem(distance_info));
    }
    if(DRightX.size() >= 2)
    {
        distance_info = getDistanceInfo(DRightX[0], DRightY[0], DRightX[1], DRightY[1], &DRightThickness);
        ui->lbDRightThick->setText(distance_info);
        tblCmp->setItem(1,2,new QTableWidgetItem(distance_info));
    }
    //show compare information
    double difference = (DCenterThickness - CCenterThickness)/CCenterThickness;
    QString s_diff = "%" + QString::number(difference*100);
    ui->lbCmpCenter->setText(s_diff);
    tblCmp->setItem(2,0,new QTableWidgetItem(s_diff));

    difference = (DLeftThickness - CLeftThickness)/CLeftThickness;
    s_diff = "%" + QString::number(difference*100);
    ui->lbCmpLeft->setText(s_diff);
    tblCmp->setItem(2,1,new QTableWidgetItem(s_diff));

    difference = (DRightThickness - CRightThickness)/CRightThickness;
    s_diff = "%" + QString::number(difference*100);
    ui->lbCmpRight->setText(s_diff);
    tblCmp->setItem(2,2,new QTableWidgetItem(s_diff));
}

void SimpleView::AutoOpenDocImage()
{
    QFileInfo fi(InputFile);
    //load back the image that processed by doctor
    QImage info_img(fi.absolutePath() + "/KneeOutput/" + fi.baseName().remove(".bmp") + "_out_load.bmp");
    //clear previous record
    DCenterX.clear(); DCenterY.clear(); DLeftX.clear(); DLeftY.clear(); DRightX.clear(); DRightY.clear();
    for (int x = 0; x < info_img.width(); x++)
    {
        for (int y = 0; y < info_img.height(); y++)
        {
            //get center point
            if(is5PixelDot(info_img, x, y, COLOR_CENTER_DOT_1))
            {
                qDebug() << "center point = (" << x+2 << "," << y+2 << ")" << endl;
                DCenterX.push_back(x+2);
                DCenterY.push_back(y+2);
            }
            //get left point
            if(is5PixelDot(info_img, x, y, COLOR_LEFT_DOT_1))
            {
                qDebug() << "left point = (" << x+2 << "," << y+2 << ")" << endl;
                DLeftX.push_back(x+2);
                DLeftY.push_back(y+2);
            }
            //get right point
            if(is5PixelDot_RightThickness(info_img, x, y, COLOR_RIGHT_DOT_1))
            {
                qDebug() << "right point = (" << x+2 << "," << y+2 << ")" << endl;
                DRightX.push_back(x+2);
                DRightY.push_back(y+2);
            }
        }
    }
    QFile statisticFile(FILE_STATISTIC);
    statisticFile.open(QIODevice::Append | QIODevice::Text);
    if(!statisticFile.isOpen()){
        qDebug() << "- Error, unable to open" << FILE_STATISTIC << "for output";
        return;
    }
    QTextStream outStream(&statisticFile);

    //show information of doctor image
    QString distance_info;
    if(DCenterX.size() >= 2)
    {
        distance_info = getDistanceInfo(DCenterX[0], DCenterY[0], DCenterX[1], DCenterY[1], &DCenterThickness);
        ui->lbDCenterThick->setText(distance_info);
        tblCmp->setItem(1,0,new QTableWidgetItem(distance_info));
        outStream << distance_info.mid(distance_info.indexOf("=")+1) + ";";
    }
    if(DLeftX.size() >= 2)
    {
        distance_info = getDistanceInfo(DLeftX[0], DLeftY[0], DLeftX[1], DLeftY[1], &DLeftThickness);
        ui->lbDLeftThick->setText(distance_info);
        tblCmp->setItem(1,1,new QTableWidgetItem(distance_info));
        outStream << distance_info.mid(distance_info.indexOf("=")+1) + ";";
    }
    if(DRightX.size() >= 2)
    {
        distance_info = getDistanceInfo(DRightX[0], DRightY[0], DRightX[1], DRightY[1], &DRightThickness);
        ui->lbDRightThick->setText(distance_info);
        tblCmp->setItem(1,2,new QTableWidgetItem(distance_info));
        outStream << distance_info.mid(distance_info.indexOf("=")+1) <<endl;
    }
    //show compare information
    double denominator = DCenterThickness > CCenterThickness ? DCenterThickness : CCenterThickness;
    double diff = abs(DCenterThickness - CCenterThickness)/denominator;
    QString s_diffC = "%" + QString::number(diff*100);
    double difference = abs(DCenterThickness - CCenterThickness);
    ui->lbCmpCenter->setText(s_diffC);
    tblCmp->setItem(2,0,new QTableWidgetItem(s_diffC));
    outStream << QString::number(difference) + ";";

    denominator = DLeftThickness > CLeftThickness ? DLeftThickness : CLeftThickness;
    diff = abs(DLeftThickness - CLeftThickness)/denominator;
    QString s_diffL = "%" + QString::number(diff*100);
    difference = abs(DLeftThickness - CLeftThickness);
    ui->lbCmpLeft->setText(s_diffL);
    tblCmp->setItem(2,1,new QTableWidgetItem(s_diffL));
    outStream << QString::number(difference) + ";";

    denominator = DRightThickness > CRightThickness ? DRightThickness : CRightThickness;
    diff = abs(DRightThickness - CRightThickness)/CRightThickness;
    QString s_diffR = "%" + QString::number(diff*100);
    difference = abs(DRightThickness - CRightThickness);
    ui->lbCmpRight->setText(s_diffR);
    tblCmp->setItem(2,2,new QTableWidgetItem(s_diffR));
    outStream << QString::number(difference) <<endl;

    outStream << s_diffC + ";" << s_diffL + ";" << s_diffR + ";" << endl;

    statisticFile.close();
}

void SimpleView::showComputerSegImage(QImage img)
{
    int w = ui->gbComputerSegImg->width(), h = ui->gbComputerSegImg->height();
    ui->lbComputerSegImg->setGeometry(0,10+ui->lbComputerSegImgFileName->height(),w,h);
    ui->lbComputerSegImg->setPixmap(QPixmap::fromImage(img).scaled(w,h,Qt::KeepAspectRatio));
    ui->lbComputerSegImgFileName->setText(InputFile);
    ui->lbComputerSegImgFileName->show();
    ui->lbComputerSegImg->show();
}

void SimpleView::showDoctorSegImage(QImage img)
{
    int w = ui->gbDoctorSegImg->width(), h = ui->gbDoctorSegImg->height();
    ui->lbDoctorSegImg->setGeometry(0,10+ui->lbDoctorSegImgFileName->height(),w,h);
    ui->lbDoctorSegImg->setPixmap(QPixmap::fromImage(img).scaled(w,h,Qt::KeepAspectRatio));
    ui->lbDoctorSegImgFileName->setText(DoctorInputFile);
    ui->lbDoctorSegImgFileName->show();
    ui->lbDoctorSegImg->show();
}

void SimpleView::slotPreProcessMrf16()
{
    typedef unsigned char InputPixelType;
    typedef itk::Image< InputPixelType, 2 > InputImageType;

    typedef float                                     OutputPixelType;
    typedef itk::Image< OutputPixelType, 2 >  OutputImageType;
    typedef itk::GradientAnisotropicDiffusionImageFilter< InputImageType,
            OutputImageType > FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput( myImage2D.originalImages());
    filter->SetNumberOfIterations(5.0/*this->ui->leAD_Iteration->text().toDouble()*/);
    filter->SetTimeStep( 0.25 );
    filter->SetConductanceParameter(5.0/*this->ui->leAD_Conductance->text().toDouble()*/);

    //convert output from float to uchar2d
    typedef itk::RescaleIntensityImageFilter<OutputImageType, InputImageType> RescaleType;
    RescaleType::Pointer rescaler = RescaleType::New();
    rescaler->SetInput( filter->GetOutput() );
    rescaler->SetOutputMinimum( itk::NumericTraits< InputPixelType >::min() );
    rescaler->SetOutputMaximum( itk::NumericTraits< InputPixelType >::max() );

    double alpha = this->ui->leSigAlpha->text().toDouble();
    double beta = this->ui->leSigBeta->text().toDouble();
    typedef itk::SigmoidImageFilter <ImageType, ImageType>
            SigmoidImageFilterType;
    SigmoidImageFilterType::Pointer sigmoidFilter
            = SigmoidImageFilterType::New();
    sigmoidFilter->SetInput(rescaler->GetOutput());
    sigmoidFilter->SetOutputMinimum(0);
    sigmoidFilter->SetOutputMaximum(itk::NumericTraits< signed short >::max());
    sigmoidFilter->SetAlpha(alpha);
    sigmoidFilter->SetBeta(beta);
    myImage2D.readFromOtherOutput(sigmoidFilter->GetOutput());
//    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
//    this->ui->qvtkWidget_Seg->repaint();
//    this->ui->qvtkWidget_Seg->show();
    slotRunMrf();

    //merge and then see whick one has largest edges
    int edge_pixel_count[this->ui->leIterationNum->text().toInt()];
    float edge_count_th_low = ImgW * 6;
    float edge_count_th_hi = ImgW * 6 * 1.4;
    cout << "ImgW : "<<ImgW <<" edge count threshold low : "<<edge_count_th_low<<" edge count threshold high: "<<edge_count_th_hi<<endl;
    for (int i = 2; i < this->ui->leIterationNum->text().toInt()-1; i++)
    {
        MrfMerge(i);
        region_growing(FILE_MRF_MERGE, FILE_REGION_GROWING, ImgW/2, ImgH-10, 255);
        edge_pixel_count[i] = sobelFilter();
        cout<<"merge:"<<i<<", edge_pixel_count:"<<edge_pixel_count[i]<<endl;
//        if (edge_pixel_count > th_edge_count) break;
    }
    for (int i = 2; i < this->ui->leIterationNum->text().toInt()-1; i++)
    {
        if (edge_pixel_count[i] > edge_count_th_low)
        {
            cout << "i:"<<i;
            for(int j = i; j < this->ui->leIterationNum->text().toInt()-1; j++)
            {
                if (edge_pixel_count[j] > edge_count_th_hi)
                {
                    cout << " j:"<<j;
                    int merge_no = (i+j-1)/2;
                    if (merge_no < i) merge_no = i;
                    MrfMerge(merge_no);
                    region_growing(FILE_MRF_MERGE, FILE_REGION_GROWING, ImgW/2, ImgH-10, 255);
                    edge_pixel_count[i] = sobelFilter();
                    cout<<" do merge:"<<merge_no<<endl;
                    break;
                }
            }
            break;
        }
    }
    drawSpline3();
    drawThickness();
}

void SimpleView::slotRunMrf()
{
    C_fileIO kmeanImage;
    C_fileIO mrfImage;
    C_itkSeg mySeg;
    mySeg.setParameter(this->ui->leIterationNum->text().toInt());
    kmeanImage.readFromOtherImage(mySeg.kmeanMethod2D( myImage2D.castsignedShort2D())) ;
    myImage2D.readFromOtherImage(  mySeg.markovMethod2D( myImage2D.castsignedShort2D(), kmeanImage.originalImages()  )  );
    myImage2D.writeImageToFile(FILE_MRF);
//    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
//    this->ui->qvtkWidget_Seg->repaint();
//    this->ui->qvtkWidget_Seg->show();
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
    filter->SetNumberOfIterations(5.0/*this->ui->leAD_Iteration->text().toDouble()*/);
    filter->SetTimeStep( 0.25 );
    filter->SetConductanceParameter(5.0/*this->ui->leAD_Conductance->text().toDouble()*/);

    //convert output from float to uchar2d
    typedef itk::RescaleIntensityImageFilter<OutputImageType, InputImageType> RescaleType;
    RescaleType::Pointer rescaler = RescaleType::New();
    rescaler->SetInput( filter->GetOutput() );
    rescaler->SetOutputMinimum( itk::NumericTraits< InputPixelType >::min() );
    rescaler->SetOutputMaximum( itk::NumericTraits< InputPixelType >::max() );

    myImage2D.readFromOtherOutput(rescaler->GetOutput());
    myImage2D.writeImageToFile(FILE_ANISOTROPIC_DIFUSSION);
//    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
//    this->ui->qvtkWidget_Seg->repaint();
//    this->ui->qvtkWidget_Seg->show();
}

void SimpleView::slotRunSig()
{
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
    myImage2D.writeImageToFile(FILE_SIGMOID);
//    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
//    this->ui->qvtkWidget_Seg->repaint();
//    this->ui->qvtkWidget_Seg->show();
}

void SimpleView::slotReset()
{
//    this->ui->qvtkWidget_Seg->hide();
    myImage2D.readFiletoImages(InputFile.toStdString().c_str());
}

void SimpleView::slotWriteFile()
{
    myImage2D.writeImageToFile(FILE_SAVE);
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

int SimpleView::sobelFilter()
{
#if 1 //escape the remove fragment step
    myImage2D.readFiletoImages(FILE_REGION_GROWING);
#else
    myImage2D.readFiletoImages(FILE_SMOOTH_EDGE);
#endif
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
    myImage2D.writeImageToFile(FILE_SOBEL);

    //convert sobel to red line
    int row, col, edge_pixel_count = 0;
    QImage inImg(FILE_SOBEL);
    QImage outImg = QImage(InputFile.toLatin1().data());
    QImage outImgC = outImg.convertToFormat(QImage::Format_RGB888);
    for (col = 0; col<outImgC.width(); col++)
    {
        for(row = 0; row < outImgC.height(); row++)
        {
            QRgb rgb = inImg.pixel(col, row);
            if(qRed(rgb) > 0/*this->ui->leTestValue->text().toInt()*/)
            {
                edge_pixel_count ++;
                outImgC.setPixel(col, row, qRgb(255, 0, 0));
            }
        }
    }
//    this->displayMyView(outImgC, None);
    outImgC.save(FILE_SOBEL_RED);
    return edge_pixel_count;
}

void SimpleView::slotSobel()
{
    sobelFilter();
}

void SimpleView::MrfMerge(int classes)
{
    QImage inImg(FILE_MRF);
    int Original_Image_Height = inImg.height();
    int Original_Image_Width = inImg.width();
    float ClassNumber;                                 //the number of classes used by MRF
    float SelectingClass;                              //how many classes would be merged?  1+2? 1+2+3? 1+2+3+4?...
    int ThresholdingValue;


    ClassNumber = this->ui->leIterationNum->text().toInt();
    SelectingClass = (float)(classes);
    ThresholdingValue = (int)(((255.0 / (ClassNumber - 1.0)) * (SelectingClass - 1.0)) + 2.0);
//    cout << "ThresholdingValue is : " << ThresholdingValue << endl;
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

//    this->displayMyView(inImg, RegionGrowing);
    QString saved_file_name = FILE_MRF_MERGE;
    inImg.save(saved_file_name.remove(".bmp") + QString::number(classes) + ".bmp");
    inImg.save(FILE_MRF_MERGE);
}

void SimpleView::slotMerge()
{
    int edge_pixel_count = 0;
    int i = this->ui->leMergedCls->text().toInt();
    MrfMerge(i);
    region_growing(FILE_MRF_MERGE, FILE_REGION_GROWING, ImgW/2, ImgH-10, 255);
//    myImage2D.readFiletoImages(FILE_MRF_MERGE);
//    //count image object
//    typedef itk::ConnectedComponentImageFilter <ImageType, ImageType >
//      ConnectedComponentImageFilterType;

//    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
//    connected->SetInput(myImage2D.originalImages());
//    connected->Update();

//    std::cout << "Number of objects: " << connected->GetObjectCount() << std::endl;

    edge_pixel_count = sobelFilter();
    cout<<"merge:"<<i<<", edge_pixel_count:"<<edge_pixel_count<<endl;
}

void SimpleView::region_growing(QString image_file, QString out_file, int seed_x, int seed_y, int replaced_pixel)
{
    ReaderType::Pointer ITKImageReader;
    ITKImageReader = ReaderType::New();
    ITKImageReader->SetFileName(image_file.toLatin1().data());
    typedef itk::ConnectedThresholdImageFilter<ImageType, ImageType> RegionGrowImageFilterType;
    RegionGrowImageFilterType::Pointer regionGrow = RegionGrowImageFilterType::New();
    float lower = 0.0;
    float upper = 50.0;
    regionGrow->SetInput(ITKImageReader->GetOutput());
    regionGrow->SetLower(lower);
    regionGrow->SetUpper(upper);

    regionGrow->SetReplaceValue(replaced_pixel);

    ImageType::IndexType seed1;
    seed1[0] = seed_x;
    seed1[1] = seed_y;
    regionGrow->SetSeed(seed1);
    regionGrow->Update();

    C_fileIO tmp_image;
    tmp_image.readFromOtherOutput(regionGrow->GetOutput());
    tmp_image.writeImageToFile(out_file.toLatin1().data());
}

void SimpleView::slotDeFragment()
{
    RemoveFragments();
}

int SimpleView::getThicknessPoint(int line2_x)
{
    //the upper spline point, for calculating the cartilage thickness
    int line1x;
    float slope1, slope2, diff = 40;
    //find the center thickness
    slope1 = (float)(SplineY2[line2_x+15] - SplineY2[line2_x-15]) / 30.0;
#ifdef DEBUG_DRAW_THICKNESS
    cout << "slope1 = " << slope1 << endl;
#endif
    if (slope1 == 0.0)
    {
        line1x = line2_x;
    }
    else
    {
        for(int x = 0; x < ImgW; x++)
        {
            slope2 = (float)(SplineY1[x] - SplineY2[line2_x]) / (float)(x - line2_x);
#ifdef DEBUG_DRAW_THICKNESS
                cout << "slope2 = " << slope2 << endl;
                cout << "diff is : " << diff << endl;
#endif
            if ((float)abs(slope2 - ((-1.0)/slope1)) < diff)
            {
                diff = (float)abs(slope2 - ((-1.0)/slope1));
                line1x = x;
            }
        }
    }
#ifdef DEBUG_DRAW_THICKNESS
    cout << "line1xy is :" << line1x << "," << SplineY1[line1x]<<endl;
#endif
    return line1x;
}

void SimpleView::slotTest()
{
//    debugImage(myImage2D.originalImages());
    drawThickness();
}

QString SimpleView::getDistanceInfo(int x1, int y1, int x2, int y2, double *dis)
{
    QString distance_info;
    double distance = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
    *dis = distance;
    distance_info.sprintf("(%d,%d)to(%d,%d)=%lf",x1,y1,x2,y2, distance);
    return distance_info;
}

QString SimpleView::getDistanceInfo(int x, double *dis)
{
    QString distance_info;
    //find 5 distances around the given point, then calculate the mean distance
    double distance = 0;
    int x1[5], x2[5], offset = 1;
    for(int i = 0; i < 5; i++)
    {
        x2[i] = x+(i-2)*offset;
        x1[i] = getThicknessPoint(x2[i]);
        cout << __func__<<"(x1,y1) = (" << x1[i] << "," << SplineY1[x1[i]] << "); " ;
        cout << "(x2,y2) = (" << x2[i] << "," << SplineY2[x2[i]] << ")" << endl;
        distance += sqrt(pow(x2[i] - x1[i], 2) + pow(SplineY2[x2[i]] - SplineY1[x1[i]], 2));
    }
    distance = distance / 5.0;
    *dis = distance;
    distance_info.sprintf("(%d,%d)to(%d,%d)=%lf",x1[2],SplineY1[x1[2]],x2[2],SplineY2[x2[2]], distance);
    return distance_info;
}

void SimpleView::drawThickness()
{
    cout << "spline bottom1 is :" << lowestY_x1 << "," << SplineY1[lowestY_x1]<<endl;
    cout << "spline bottom2 is :" << lowestY_x2 << "," << SplineY2[lowestY_x2]<<endl;
    int line1_center_x, line1_left_x, line1_right_x, line2_left_x, line2_right_x;
    line2_left_x = lowestY_x2/3;
    line2_right_x = lowestY_x2 + (ImgW-lowestY_x2)*2/3;
    line1_center_x = getThicknessPoint(lowestY_x2);
    line1_left_x = getThicknessPoint(line2_left_x);
    line1_right_x = getThicknessPoint(line2_right_x);
    QImage img(FILE_SPLINE);
    QPainter pt(&img);
    int x1, y1, x2, y2;
    QString distance_info;
    pt.setPen(Qt::green);

    QFile statisticFile(FILE_STATISTIC);
    statisticFile.open(QIODevice::Append | QIODevice::Text);
    if(!statisticFile.isOpen()){
        qDebug() << "- Error, unable to open" << FILE_STATISTIC << "for output";
        return;
    }
    QTextStream outStream(&statisticFile);
    outStream << QFileInfo(InputFile).baseName() + "; ; " << endl;

    //draw center and distance info
    x1 = line1_center_x;    y1 = SplineY1[line1_center_x];
    x2 = lowestY_x2;        y2 = SplineY2[lowestY_x2];
    distance_info = getDistanceInfo(x2, &CCenterThickness);
    ui->lbCCenterThick->setText(distance_info);
    tblCmp->setItem(0,0,new QTableWidgetItem(distance_info));
    outStream << distance_info.mid(distance_info.indexOf("=")+1) + ";";
    pt.drawLine(x1,y1, x2,y2);
    pt.setPen(Qt::yellow);
    pt.drawText(QRect(0, 0, 500, 20),Qt::AlignLeft,distance_info);
    pt.setPen(Qt::green);

    //draw left and distance info
    x1 = line1_left_x;      y1 = SplineY1[line1_left_x];
    x2 = line2_left_x;      y2 = SplineY2[line2_left_x];
    distance_info = getDistanceInfo(x2, &CLeftThickness);
    ui->lbCLeftThick->setText(distance_info);
    tblCmp->setItem(0,1,new QTableWidgetItem(distance_info));
    outStream << distance_info.mid(distance_info.indexOf("=")+1) + ";";
    pt.drawLine(x1,y1, x2,y2);
    pt.setPen(Qt::yellow);
    pt.drawText(QRect(0, 20, 300, 20),Qt::AlignLeft,distance_info);
    pt.setPen(Qt::green);

    //draw right and distance info
    pt.drawLine(line1_right_x,SplineY1[line1_right_x], line2_right_x,SplineY2[line2_right_x]);
    x1 = line1_right_x;      y1 = SplineY1[line1_right_x];
    x2 = line2_right_x;      y2 = SplineY2[line2_right_x];
    distance_info = getDistanceInfo(x2, &CRightThickness);
    ui->lbCRightThick->setText(distance_info);
    tblCmp->setItem(0,2,new QTableWidgetItem(distance_info));
    outStream << distance_info.mid(distance_info.indexOf("=")+1) <<endl;
    pt.drawLine(x1,y1, x2,y2);
    pt.setPen(Qt::yellow);
    pt.drawText(QRect(0, 40, 300, 30),Qt::AlignLeft,distance_info);
    pt.setPen(Qt::green);

    pt.end();
//    displayMyView(img, CalculateDistance);
    showComputerSegImage(img);
    img.save(FILE_THICKNESS);
    statisticFile.close();
}

void SimpleView::RemoveFragments()
{
    //merge other fragments to the cartilage
    QImage ori_image(FILE_REGION_GROWING);
    int width = ori_image.width();
    int height = ori_image.height();
    cout << "w is: "<<width<<endl;
    cout << "h is: "<<height<<endl;
    //generate the top part of the region growing image
    bool b_find_region = false;
    for(int y = 0; y < height; y++)
    {
        for(int x = 0; x < width; x++)
        {
            QRgb pixel = ori_image.pixel(x, y);
            if (qRed(pixel) == 0)
            {
                region_growing(FILE_REGION_GROWING, FILE_REGION_GROWING_TOP, x, y, 150);
                b_find_region = true;
                cout << "x1 is: "<<x<<endl;
                cout << "y1 is: "<<y<<endl;
                break;
            }
        }
        if(b_find_region) break;
    }
    //generate the bottom part of the region growing image
    b_find_region = false;
    for(int y = height - 1; y >=0; y--)
    {
        for(int x = 0; x < width; x++)
        {
            QRgb pixel = ori_image.pixel(x, y);
            if (qRed(pixel) == 0)
            {
                region_growing(FILE_REGION_GROWING, FILE_REGION_GROWING_BOT, x, y, 150);
                b_find_region = true;
                cout << "x2 is: "<<x<<endl;
                cout << "y2 is: "<<y<<endl;
                break;
            }
        }
        if(b_find_region) break;
    }

    //merge top & bottom part to a new image as to remove fragments
    QImage top_img(FILE_REGION_GROWING_TOP);
    QImage bot_img(FILE_REGION_GROWING_BOT);
    QImage out_img(width, height, QImage::Format_RGB888);
    for(int x = 0; x < width; x++)
    {
        for(int y = 0; y < height; y++)
        {
            if (qRed(top_img.pixel(x, y)) == 150 || qRed(bot_img.pixel(x, y)))
            {
                out_img.setPixel(x,y,qRgb(0,0,0));
            }
            else
            {
                out_img.setPixel(x,y,qRgb(255,255,255));
            }
        }
    }
    out_img.save(FILE_REMOVE_FRAGMENT);
    displayMyView(out_img, None);
//    Opening();

//    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
//    this->ui->qvtkWidget_Seg->repaint();
//    this->ui->qvtkWidget_Seg->show();

}

void SimpleView::slotSmoothEdge()
{
    ReaderType::Pointer ITKImageReader;
    ITKImageReader = ReaderType::New();
#if 1 //escape the remove fragment step
    ITKImageReader->SetFileName(FILE_REGION_GROWING);
#else
    ITKImageReader->SetFileName(FILE_REMOVE_FRAGMENT);
#endif
    typedef itk::BinaryBallStructuringElement<ImageType::PixelType, ImageType::ImageDimension>
                StructuringElementType;

    StructuringElementType structuringElement;

    typedef itk::BinaryErodeImageFilter<ImageType, ImageType, StructuringElementType>
            BinaryErodeImageFilterType;
    BinaryErodeImageFilterType::Pointer erodeFilter
            = BinaryErodeImageFilterType::New();
    typedef itk::BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType>
            BinaryDilateImageFilterType;
    BinaryDilateImageFilterType::Pointer dilateFilter
            = BinaryDilateImageFilterType::New();
    //erode 2
    structuringElement.SetRadius(this->ui->leSmoothArg1->text().toInt());
    structuringElement.CreateStructuringElement();
    erodeFilter->SetInput(ITKImageReader->GetOutput());
    erodeFilter->SetKernel(structuringElement);
    erodeFilter->Update();
    //dilate 5
    structuringElement.SetRadius(this->ui->leSmoothArg2->text().toInt());
    structuringElement.CreateStructuringElement();
    dilateFilter->SetInput(erodeFilter->GetOutput());
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->Update();
    //erode 3
    structuringElement.SetRadius(this->ui->leSmoothArg3->text().toInt());
    structuringElement.CreateStructuringElement();
    erodeFilter->SetInput(dilateFilter->GetOutput());
    erodeFilter->SetKernel(structuringElement);
    erodeFilter->Update();
    myImage2D.readFromOtherOutput(erodeFilter->GetOutput());
    myImage2D.writeImageToFile(FILE_SMOOTH_EDGE);
    QImage img(FILE_SMOOTH_EDGE);
    displayMyView(img, None);

//    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
//    this->ui->qvtkWidget_Seg->repaint();
//    this->ui->qvtkWidget_Seg->show();
}

void SimpleView::Opening()
{
#if 1 //openinig & closing test
    ReaderType::Pointer ITKImageReader;
    ITKImageReader = ReaderType::New();
    ITKImageReader->SetFileName(FILE_REMOVE_FRAGMENT);
    typedef itk::BinaryBallStructuringElement<ImageType::PixelType, ImageType::ImageDimension>
                StructuringElementType;

    //closing
//    StructuringElementType structuringElement_closing;
//    structuringElement_closing.SetRadius(5);
//    structuringElement_closing.CreateStructuringElement();
//    typedef itk::BinaryMorphologicalClosingImageFilter <ImageType, ImageType, StructuringElementType>
//            BinaryMorphologicalClosingImageFilterType;
//    BinaryMorphologicalClosingImageFilterType::Pointer closingFilter
//            = BinaryMorphologicalClosingImageFilterType::New();
//    closingFilter->SetInput(ITKImageReader->GetOutput());
//    closingFilter->SetKernel(structuringElement_closing);
//    closingFilter->Update();
//    myImage2D.readFromOtherOutput(closingFilter->GetOutput());

    //opening
    StructuringElementType structuringElement_opening;
    structuringElement_opening.SetRadius(6);
    structuringElement_opening.CreateStructuringElement();
    typedef itk::BinaryMorphologicalOpeningImageFilter <ImageType, ImageType, StructuringElementType>
            BinaryMorphologicalOpeningImageFilterType;
    BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter
            = BinaryMorphologicalOpeningImageFilterType::New();
    openingFilter->SetInput(ITKImageReader->GetOutput());
    openingFilter->SetKernel(structuringElement_opening);
    openingFilter->Update();
    myImage2D.readFromOtherOutput(openingFilter->GetOutput());
    myImage2D.writeImageToFile(FILE_REMOVE_FRAGMENT);

//    this->ui->qvtkWidget_Seg->GetRenderWindow()->AddRenderer(myImage2D.vtkRender());
//    this->ui->qvtkWidget_Seg->repaint();
//    this->ui->qvtkWidget_Seg->show();
#endif
}

void SimpleView::slotSpline()
{
#if 1 //the difference of 2 slope must < a constant
    drawSpline3();
#else
    drawSpline();
#endif
}

//first version of drawSpline, just only check the slope between the neighor
void SimpleView::drawSpline()
{
#if 1 //test
    QImage inImg = QImage(FILE_SOBEL_RED);
    QImage outImg = QImage(InputFile.toLatin1().data());
    QImage outImgC = outImg.convertToFormat(QImage::Format_RGB888);
    std::vector<double> x1, y1, x2, y2;
    tk::spline s1, s2;
    int last_red_y1 = -1, last_red_x1 = -1;
    int last_red_y2 = -1, last_red_x2 = -1;
    float slope = 0.0;
    for(int x = 0; x < inImg.width(); x += 20)
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
    inImg.save(FILE_SPLINE_SAMPLE);
    s1.set_points(x1,y1);
    s2.set_points(x2,y2);
    memset(SplineY1, 0, sizeof(SplineY1));
    memset(SplineY2, 0, sizeof(SplineY2));
    lowestY_x1 = lowestY_x2 = 0;
    int max_y1 = 0, max_y2 = 0;
    for (int col = 0; col<outImgC.width(); col++)
    {
        outImgC.setPixel(col, (int)s1((double)col), qRgb(255, 0, 0));
        outImgC.setPixel(col, (int)s2((double)col), qRgb(255, 0, 0));
        //recording the 2 spline cordinates, and lowestY_x1, lowestY_x2
        SplineY1[col] = (int)s1((double)col);
        SplineY2[col] = (int)s2((double)col);
        if (SplineY1[col] > max_y1)
        {
            max_y1 = SplineY1[col];
            lowestY_x1 = col;
        }
        if (SplineY2[col] > max_y2)
        {
            max_y2 = SplineY2[col];
            lowestY_x2 = col;
        }
    }
    drawColorDot(&outImgC, qRgb(0,0,255),lowestY_x1, SplineY1[lowestY_x1]);
    drawColorDot(&outImgC, qRgb(0,255,0),lowestY_x2, SplineY2[lowestY_x2]);
    this->displayMyView(outImgC, None);
    outImgC.save(FILE_SPLINE);
#endif
}

void SimpleView::drawSpline1()
{
    QImage inImg = QImage(FILE_SOBEL_RED);
    QImage sampleImg = QImage(FILE_SOBEL_RED);
    QImage outImg = QImage(InputFile.toLatin1().data());
    QImage outImgC = outImg.convertToFormat(QImage::Format_RGB888);
    std::vector<double> sp_x1, sp_y1, sp_x2, sp_y2;
    tk::spline s1, s2;
    int last_red_y1 = -1, last_red_x1 = -1;
    int last_red_y2 = -1, last_red_x2 = -1;
    float slope = 0.0, old_slope1 = 0.0, old_slope2 = 0.0;
    for(int x = 0; x < inImg.width(); x += this->ui->leSplineW->text().toInt())
    {
        for(int y = 0; y < inImg.height(); y++)
        {
            QRgb rgb = inImg.pixel(x, y);
            if(qRed(rgb) >=255)
            {
                slope = (float)(y-last_red_y1) / (float)(x-last_red_x1);
                cout<<"############Line1\n";
                cout<<"(x,y) = ("<<x<<","<<y<<")";
                cout<<", (last_red_x1, last_red_y1) = ("<<last_red_x1<<","<<last_red_y1<<")";
                cout<<", slope = "<<slope;
                cout<<", old_slope1 = "<<old_slope1<<endl;

                //the skip point condition
                    //if slope is going up and before width/2, skip
                if (x < (inImg.height()/2 + 10) && slope < 0)
                    break;

                if(last_red_y1 == -1 || //must sample first point
//                   last_red_x1 == 0  || //must sample second point, let the old_slope1 be correct
                   //the diffefenct of 2 slope must < 0.6 && itself must < 0.6
                   fabs(slope-old_slope1) < 0.6 && fabs(slope) < 0.65 ||
                   //slope begin to going up && can not going up too much (< 0.6)
                   (slope < 0 && old_slope1 > 0) && fabs(slope) < 0.55 /*||
                   //
                   x > inImg.height()/5 && slope > 0*/)
                {
                    cout<<"sample taken!!"<<endl;
                    sp_x1.push_back(x);
                    sp_y1.push_back(y);
                    drawColorDot(&sampleImg,qRgb(0,0,255), x, y);
                    last_red_x1 = x;
                    last_red_y1 = y;
                    old_slope1 = slope;
                }
                //get first old_slope1
                if(last_red_x1 == 0) old_slope1 = slope;
                break;
            }
        }
        //find top line of the bottom edge
        int y2 = 0, y_use = 0, y_cnt = 0, y_tmp = 0;
        bool every_white = false;
        for(int y = inImg.height()-1; y >= 0; y--)
        {
            QRgb rgb = inImg.pixel(x, y);
            if(qRed(rgb) >=255)
            {
                y_cnt ++;
                cout<<"############Line2\n";
                cout<<"(x,y) = ("<<x<<","<<y<<")";
                for (int i = y; i >=0; i--)
                {
                    QRgb t_rgb = inImg.pixel(x, i);
                    if (qGreen(t_rgb) > 0) every_white = true;
                    if(qRed(t_rgb) >=255 && every_white)
                    {
                        every_white = false;
                        y_cnt ++;
                        if (y_cnt == 2) y2 = i; //remember the line2's Y cordinate
                        cout<<", (x,i) = ("<<x<<","<<i<<")";
                        cout<<", y_cnt = "<<y_cnt<<endl;
                    }
                    else if (qRed(t_rgb) >=255 && !every_white)
                    {
                        //to remember the top y of the continous red y points
                        if (y_tmp == 0 || i == y_tmp - 1)
                            y_tmp = i;
                        cout<<", y_tmp = "<<y_tmp<<endl;
                    }
                }
                if (y_cnt >= 3) //means total >= 3 lines, but I needs line2
                {
                    y_use = y2;
                }
                else if (y_cnt == 2) //means total has 2 lines, but I needs line1
                {
                    y_use = y_tmp;
                }
                cout<<", y_cnt = "<<y_cnt<<endl;

                slope = (float)(y_use-last_red_y2) / (float)(x-last_red_x2);
                cout<<"(x,y_use) = ("<<x<<","<<y_use<<")";
                cout<<", (last_red_x2, last_red_y2) = ("<<last_red_x2<<","<<last_red_y2<<")";
                cout<<", slope = "<<slope;
                cout<<", old_slope2 = "<<old_slope2<<endl;
                if(last_red_y2 == -1 ||
                   last_red_x2 == 0  ||
                   //the diffefenct of 2 slope must < 0.6 && itself must < 0.6
                   fabs(slope-old_slope2) < 0.6  && fabs(slope) < 0.65 ||
                   //slope begin to going up && can not going up too much (< 0.6)
                   (slope < 0 && old_slope2 > 0) && fabs(slope) < 0.6)
                {
                    cout<<"sample taken!!"<<endl;
                    sp_x2.push_back(x);
                    sp_y2.push_back(y_use);
                    drawColorDot(&sampleImg,qRgb(0,255,0), x, y_use);
                    last_red_x2 = x;
                    last_red_y2 = y_use;
                    old_slope2 = slope;
                }
                break;
            }
        }
    }
    s1.set_points(sp_x1,sp_y1);
    s2.set_points(sp_x2,sp_y2);
    memset(SplineY1, 0, sizeof(SplineY1));
    memset(SplineY2, 0, sizeof(SplineY2));
    lowestY_x1 = lowestY_x2 = 0;
    int max_y1 = 0, max_y2 = 0;
    for (int col = 0; col<outImgC.width(); col++)
    {
        outImgC.setPixel(col, (int)s1((double)col), qRgb(255, 0, 0));
        outImgC.setPixel(col, (int)s2((double)col), qRgb(255, 0, 0));
        drawColorDot(&inImg, qRgb(0,0,255), col, (int)s1((double)col));
        drawColorDot(&inImg, qRgb(0,255,0), col, (int)s2((double)col));
        //recording the 2 spline cordinates, and lowestY_x1, lowestY_x2
        SplineY1[col] = (int)s1((double)col);
        SplineY2[col] = (int)s2((double)col);
        if (SplineY1[col] > max_y1)
        {
            max_y1 = SplineY1[col];
            lowestY_x1 = col;
        }
        if (SplineY2[col] > max_y2)
        {
            max_y2 = SplineY2[col];
            lowestY_x2 = col;
        }
    }
    sampleImg.save(FILE_SPLINE_SAMPLE);
    drawColorDot(&outImgC, qRgb(0,0,255),lowestY_x1, SplineY1[lowestY_x1]);
    drawColorDot(&outImgC, qRgb(0,255,0),lowestY_x2, SplineY2[lowestY_x2]);
    this->displayMyView(outImgC, None);
    outImgC.save(FILE_SPLINE);
}

//clone from drawSpline1, add reverse check
void SimpleView::drawSpline3()
{
    QImage inImg = QImage(FILE_SOBEL_RED);
    QImage sampleImg = QImage(FILE_SOBEL_RED);
    QImage outImg = QImage(InputFile.toLatin1().data());
    QImage outImgC = outImg.convertToFormat(QImage::Format_RGB888);
    std::vector<double> sp_x1, sp_y1, sp_x2, sp_y2;
    std::vector<double> sp_all_x1, sp_all_y1, sp_all_x2, sp_all_y2;
    tk::spline s1, s2;
    int last_red_y1 = -1, last_red_x1 = -1;
    int last_red_y2 = -1, last_red_x2 = -1;
    float slope = 0.0, old_slope1 = 0.0, old_slope2 = 0.0;
    //1.get the all sample point first
    for(int x = 0; x < inImg.width(); x += this->ui->leSplineW->text().toInt())
    {
        //line1
        for(int y = 0; y < inImg.height(); y++)
        {
            QRgb rgb = inImg.pixel(x, y);
            if(qRed(rgb) >=255)
            {
                sp_all_x1.push_back(x);
                sp_all_y1.push_back(y);
                break;
            }
        }

        //line2
        //find top line of the bottom edge
        int y2 = 0, y_use = 0, y_cnt = 0, y_tmp = 0;
        bool every_white = false;
        for(int y = inImg.height()-1; y >= 0; y--)
        {
            QRgb rgb = inImg.pixel(x, y);
            if(qRed(rgb) >=255)
            {
                y_cnt ++;
#ifdef DEBUG_DRAW_THICKNESS
                cout<<"############Line2\n";
                cout<<"(x,y) = ("<<x<<","<<y<<")";
#endif
                //find the top point of line2
                for (int i = y; i >=0; i--)
                {
                    QRgb t_rgb = inImg.pixel(x, i);
                    if (qGreen(t_rgb) > 0) every_white = true;
                    if(qRed(t_rgb) >=255 && every_white)
                    {
                        every_white = false;
                        y_cnt ++;
                        if (y_cnt == 2) y2 = i; //remember the line2's Y cordinate
#ifdef DEBUG_DRAW_THICKNESS
                        cout<<", (x,i) = ("<<x<<","<<i<<")";
                        cout<<", y_cnt = "<<y_cnt<<endl;
#endif
                    }
                    else if (qRed(t_rgb) >=255 && !every_white)
                    {
                        //to remember the top y of the continous red y points
                        if (y_tmp == 0 || i == y_tmp - 1)
                            y_tmp = i;
#ifdef DEBUG_DRAW_THICKNESS
                        cout<<", y_tmp = "<<y_tmp<<endl;
#endif
                    }
                }
                if (y_cnt >= 3) //means total >= 3 lines, but I needs line2
                {
                    y_use = y2;
                }
                else if (y_cnt == 2) //means total has 2 lines, but I needs line1
                {
                    y_use = y_tmp;
                }
#ifdef DEBUG_DRAW_THICKNESS
                cout<<", y_cnt = "<<y_cnt<<endl;
#endif
                sp_all_x2.push_back(x);
                sp_all_y2.push_back(y_use);
                break;
            }
        }
    }
    //2.to check to dispose first point or not
    float slope_first_point = (sp_all_y1[1] - sp_all_y1[0])/(sp_all_x1[1] - sp_all_x1[0]);
#ifdef DEBUG_DRAW_THICKNESS
    cout<<"point0 :"<<sp_all_x1[0]<<","<<sp_all_y1[0]<<endl;
    cout<<"point1 :"<<sp_all_x1[1]<<","<<sp_all_y1[1]<<endl;
    cout<<"slope_first_point :"<<slope_first_point<<endl;
#endif
    //if first slope too large, delete the first point then interpolate it
    if (slope_first_point > 0.65)
    {
#ifdef DEBUG_DRAW_THICKNESS
        cout<<"######after interpolate first point:"<<endl;
#endif
        int interpolate_y = 0;
        for (int i = 1; i <=5; i++)
        {
            interpolate_y += sp_all_y1[i+1] - sp_all_y1[i];
        }
        //average
        interpolate_y = interpolate_y/5;
        sp_all_y1[0] = sp_all_y1[1] - interpolate_y;
        slope_first_point = (sp_all_y1[1] - sp_all_y1[0])/(sp_all_x1[1] - sp_all_x1[0]);
#ifdef DEBUG_DRAW_THICKNESS
        cout<<"point0 :"<<sp_all_x1[0]<<","<<sp_all_y1[0]<<endl;
        cout<<"point1 :"<<sp_all_x1[1]<<","<<sp_all_y1[1]<<endl;
        cout<<"slope_first_point :"<<slope_first_point<<endl;
#endif
    }

    //3.remove some unsuitable point of line1 & line2
    for (int i = 0; i < sp_all_x1.size(); i++)
    {
        /*****line1*****/
        int x = sp_all_x1[i], y = sp_all_y1[i];
        slope = (float)(y-last_red_y1) / (float)(x-last_red_x1);
        //get first old_slope1
        if(last_red_x1 == 0) old_slope1 = slope;
#ifdef DEBUG_DRAW_THICKNESS
        cout<<"############Line1\n";
        cout<<"(x,y) = ("<<x<<","<<y<<")";
        cout<<", (last_red_x1, last_red_y1) = ("<<last_red_x1<<","<<last_red_y1<<")";
        cout<<", slope = "<<slope;
        cout<<", old_slope1 = "<<old_slope1<<endl;
#endif
        //the skip point condition
            //if slope is going up and before width/2, skip
        if (x < (inImg.height()/2 + 10) && slope < 0)
            continue;

        if(last_red_y1 == -1 || //must sample first point
//                   last_red_x1 == 0  || //must sample second point, let the old_slope1 be correct
            //slope begin to going up && can not going up too much (< 0.6)
            (slope < 0 && old_slope1 > 0) && fabs(slope) < 0.55 ||
           //the diffefenct of 2 slope must < 0.6 && itself must < 0.6
           fabs(slope-old_slope1) < 0.6 && fabs(slope) < 0.65
                )
        {
#ifdef DEBUG_DRAW_THICKNESS
            cout<<"sample taken!!"<<endl;
#endif
            sp_x1.push_back(x);
            sp_y1.push_back(y);
            drawColorDot(&sampleImg,qRgb(0,0,255), x, y);
            last_red_x1 = x;
            last_red_y1 = y;
            old_slope1 = slope;
        }
    }

    for (int i = 0; i < sp_all_x2.size(); i++)
    {
        /*****line2*****/
        int x = sp_all_x2[i], y = sp_all_y2[i];
        slope = (float)(y-last_red_y2) / (float)(x-last_red_x2);
#ifdef DEBUG_DRAW_THICKNESS
        cout<<"(x,y) = ("<<x<<","<<y<<")";
        cout<<", (last_red_x2, last_red_y2) = ("<<last_red_x2<<","<<last_red_y2<<")";
        cout<<", slope = "<<slope;
        cout<<", old_slope2 = "<<old_slope2<<endl;
#endif
        //if slope is going up and before width/2, skip
        if (x < (inImg.height()/2 + 10) && slope < 0)
            continue;
        if(last_red_y2 == -1 ||
           last_red_x2 == 0  ||
            //slope begin to going up && can not going up too much (< 0.6)
            (slope < 0 && old_slope2 > 0) && fabs(slope) < 0.6 ||
           //the diffefenct of 2 slope must < 0.6 && itself must < 0.6
           fabs(slope-old_slope2) < 0.6  && fabs(slope) < 0.65
                )
        {
#ifdef DEBUG_DRAW_THICKNESS
            cout<<"sample taken!!"<<endl;
#endif
            sp_x2.push_back(x);
            sp_y2.push_back(y);
            drawColorDot(&sampleImg,qRgb(0,255,0), x, y);
            last_red_x2 = x;
            last_red_y2 = y;
            old_slope2 = slope;
        }
    }

    s1.set_points(sp_x1,sp_y1);
    s2.set_points(sp_x2,sp_y2);
    memset(SplineY1, 0, sizeof(SplineY1));
    memset(SplineY2, 0, sizeof(SplineY2));
    lowestY_x1 = lowestY_x2 = 0;
    int max_y1 = 0, max_y2 = 0;
    for (int col = 0; col<outImgC.width(); col++)
    {
        outImgC.setPixel(col, (int)s1((double)col), qRgb(255, 0, 0));
        outImgC.setPixel(col, (int)s2((double)col), qRgb(255, 0, 0));
        drawColorDot(&inImg, qRgb(0,0,255), col, (int)s1((double)col));
        drawColorDot(&inImg, qRgb(0,255,0), col, (int)s2((double)col));
        //recording the 2 spline cordinates, and lowestY_x1, lowestY_x2
        SplineY1[col] = (int)s1((double)col);
        SplineY2[col] = (int)s2((double)col);
        if (SplineY1[col] > max_y1)
        {
            max_y1 = SplineY1[col];
            lowestY_x1 = col;
        }
        if (SplineY2[col] > max_y2)
        {
            max_y2 = SplineY2[col];
            lowestY_x2 = col;
        }
    }
    sampleImg.save(FILE_SPLINE_SAMPLE);
    //draw lowest point
    drawColorDot(&outImgC, qRgb(0,0,255),lowestY_x1, SplineY1[lowestY_x1]);
    drawColorDot(&outImgC, qRgb(0,255,0),lowestY_x2, SplineY2[lowestY_x2]);
//    this->displayMyView(outImgC, None);
    outImgC.save(FILE_SPLINE);
}

//clone from drawSpline3, add reverse check
void SimpleView::drawSpline4()
{
    QImage inImg = QImage(FILE_SOBEL_RED);
    QImage sampleImg = QImage(FILE_SOBEL_RED);
    QImage outImg = QImage(InputFile.toLatin1().data());
    QImage outImgC = outImg.convertToFormat(QImage::Format_RGB888);
    std::vector<double> sp_x1, sp_y1, sp_x2, sp_y2;
    std::vector<double> sp_all_x1, sp_all_y1, sp_all_x2, sp_all_y2;
    tk::spline s1, s2;
    int last_red_y1 = -1, last_red_x1 = -1;
    int last_red_y2 = -1, last_red_x2 = -1;
    float slope = 0.0, old_slope1 = 0.0, old_slope2 = 0.0;
    //1.get the all sample point first
    for(int x = 0; x < inImg.width(); x += this->ui->leSplineW->text().toInt())
    {
        //line1
        for(int y = 0; y < inImg.height(); y++)
        {
            QRgb rgb = inImg.pixel(x, y);
            if(qRed(rgb) >=255)
            {
                sp_all_x1.push_back(x);
                sp_all_y1.push_back(y);
                break;
            }
        }

        //line2
        //find top line of the bottom edge
        int y2 = 0, y_use = 0, y_cnt = 0, y_tmp = 0;
        bool every_white = false;
        for(int y = inImg.height()-1; y >= 0; y--)
        {
            QRgb rgb = inImg.pixel(x, y);
            if(qRed(rgb) >=255)
            {
                y_cnt ++;
#ifdef DEBUG_DRAW_THICKNESS
                cout<<"############Line2\n";
                cout<<"(x,y) = ("<<x<<","<<y<<")";
#endif
                //find the top point of line2
                for (int i = y; i >=0; i--)
                {
                    QRgb t_rgb = inImg.pixel(x, i);
                    if (qGreen(t_rgb) > 0) every_white = true;
                    if(qRed(t_rgb) >=255 && every_white)
                    {
                        every_white = false;
                        y_cnt ++;
                        if (y_cnt == 2) y2 = i; //remember the line2's Y cordinate
#ifdef DEBUG_DRAW_THICKNESS
                        cout<<", (x,i) = ("<<x<<","<<i<<")";
                        cout<<", y_cnt = "<<y_cnt<<endl;
#endif
                    }
                    else if (qRed(t_rgb) >=255 && !every_white)
                    {
                        //to remember the top y of the continous red y points
                        if (y_tmp == 0 || i == y_tmp - 1)
                            y_tmp = i;
#ifdef DEBUG_DRAW_THICKNESS
                        cout<<", y_tmp = "<<y_tmp<<endl;
#endif
                    }
                }
                if (y_cnt >= 3) //means total >= 3 lines, but I needs line2
                {
                    y_use = y2;
                }
                else if (y_cnt == 2) //means total has 2 lines, but I needs line1
                {
                    y_use = y_tmp;
                }
#ifdef DEBUG_DRAW_THICKNESS
                cout<<", y_cnt = "<<y_cnt<<endl;
#endif
                sp_all_x2.push_back(x);
                sp_all_y2.push_back(y_use);
                break;
            }
        }
    }
    //2.to check to dispose first point or not
    float slope_first_point = (sp_all_y1[1] - sp_all_y1[0])/(sp_all_x1[1] - sp_all_x1[0]);
#ifdef DEBUG_DRAW_THICKNESS
    cout<<"point0 :"<<sp_all_x1[0]<<","<<sp_all_y1[0]<<endl;
    cout<<"point1 :"<<sp_all_x1[1]<<","<<sp_all_y1[1]<<endl;
    cout<<"slope_first_point :"<<slope_first_point<<endl;
#endif
    //if first slope too large, delete the first point then interpolate it
    if (slope_first_point > 0.65)
    {
#ifdef DEBUG_DRAW_THICKNESS
        cout<<"######after interpolate first point:"<<endl;
#endif
        int interpolate_y = 0;
        for (int i = 1; i <=5; i++)
        {
            interpolate_y += sp_all_y1[i+1] - sp_all_y1[i];
        }
        //average
        interpolate_y = interpolate_y/5;
        sp_all_y1[0] = sp_all_y1[1] - interpolate_y;
        slope_first_point = (sp_all_y1[1] - sp_all_y1[0])/(sp_all_x1[1] - sp_all_x1[0]);
#ifdef DEBUG_DRAW_THICKNESS
        cout<<"point0 :"<<sp_all_x1[0]<<","<<sp_all_y1[0]<<endl;
        cout<<"point1 :"<<sp_all_x1[1]<<","<<sp_all_y1[1]<<endl;
        cout<<"slope_first_point :"<<slope_first_point<<endl;
#endif
    }

    //3.remove some unsuitable point of line1 & line2
    for (int i = 0; i < sp_all_x1.size(); i++)
    {
        /*****line1*****/
        int x = sp_all_x1[i], y = sp_all_y1[i];
        slope = (float)(y-last_red_y1) / (float)(x-last_red_x1);
        //get first old_slope1
        if(last_red_x1 == 0) old_slope1 = slope;
#ifdef DEBUG_DRAW_THICKNESS
        cout<<"############Line1\n";
        cout<<"(x,y) = ("<<x<<","<<y<<")";
        cout<<", (last_red_x1, last_red_y1) = ("<<last_red_x1<<","<<last_red_y1<<")";
        cout<<", slope = "<<slope;
        cout<<", old_slope1 = "<<old_slope1<<endl;
#endif
        //the skip point condition
            //if slope is going up and before width/2, skip
        if (x < (inImg.height()/2 + 10) && slope < 0)
            continue;

        if(last_red_y1 == -1 || //must sample first point
//                   last_red_x1 == 0  || //must sample second point, let the old_slope1 be correct
            //slope begin to going up && can not going up too much (< 0.6)
            (slope < 0 && old_slope1 > 0) && fabs(slope) < 0.55 ||
           //the diffefenct of 2 slope must < 0.6 && itself must < 0.6
           fabs(slope-old_slope1) < 0.6 && fabs(slope) < 0.65
                )
        {
#ifdef DEBUG_DRAW_THICKNESS
            cout<<"sample taken!!"<<endl;
#endif
            sp_x1.push_back(x);
            sp_y1.push_back(y);
            drawColorDot(&sampleImg,qRgb(0,0,255), x, y);
            last_red_x1 = x;
            last_red_y1 = y;
            old_slope1 = slope;
        }
    }

    for (int i = 0; i < sp_all_x2.size(); i++)
    {
        /*****line2*****/
        int x = sp_all_x2[i], y = sp_all_y2[i];
        slope = (float)(y-last_red_y2) / (float)(x-last_red_x2);
#ifdef DEBUG_DRAW_THICKNESS
        cout<<"(x,y) = ("<<x<<","<<y<<")";
        cout<<", (last_red_x2, last_red_y2) = ("<<last_red_x2<<","<<last_red_y2<<")";
        cout<<", slope = "<<slope;
        cout<<", old_slope2 = "<<old_slope2<<endl;
#endif
        //if slope is going up and before width/2, skip
        if (x < (inImg.height()/2 + 10) && slope < 0)
            continue;
        if(last_red_y2 == -1 ||
           last_red_x2 == 0  ||
            //slope begin to going up && can not going up too much (< 0.6)
            (slope < 0 && old_slope2 > 0) && fabs(slope) < 0.6 ||
           //the diffefenct of 2 slope must < 0.6 && itself must < 0.6
           fabs(slope-old_slope2) < 0.6  && fabs(slope) < 0.65
                )
        {
#ifdef DEBUG_DRAW_THICKNESS
            cout<<"sample taken!!"<<endl;
#endif
            sp_x2.push_back(x);
            sp_y2.push_back(y);
            drawColorDot(&sampleImg,qRgb(0,255,0), x, y);
            last_red_x2 = x;
            last_red_y2 = y;
            old_slope2 = slope;
        }
    }

    s1.set_points(sp_x1,sp_y1);
    s2.set_points(sp_x2,sp_y2);
    memset(SplineY1, 0, sizeof(SplineY1));
    memset(SplineY2, 0, sizeof(SplineY2));
    lowestY_x1 = lowestY_x2 = 0;
    int max_y1 = 0, max_y2 = 0;
    for (int col = 0; col<outImgC.width(); col++)
    {
        outImgC.setPixel(col, (int)s1((double)col), qRgb(255, 0, 0));
        outImgC.setPixel(col, (int)s2((double)col), qRgb(255, 0, 0));
        drawColorDot(&inImg, qRgb(0,0,255), col, (int)s1((double)col));
        drawColorDot(&inImg, qRgb(0,255,0), col, (int)s2((double)col));
        //recording the 2 spline cordinates, and lowestY_x1, lowestY_x2
        SplineY1[col] = (int)s1((double)col);
        SplineY2[col] = (int)s2((double)col);
        if (SplineY1[col] > max_y1)
        {
            max_y1 = SplineY1[col];
            lowestY_x1 = col;
        }
        if (SplineY2[col] > max_y2)
        {
            max_y2 = SplineY2[col];
            lowestY_x2 = col;
        }
    }
    sampleImg.save(FILE_SPLINE_SAMPLE);
    //draw lowest point
    drawColorDot(&outImgC, qRgb(0,0,255),lowestY_x1, SplineY1[lowestY_x1]);
    drawColorDot(&outImgC, qRgb(0,255,0),lowestY_x2, SplineY2[lowestY_x2]);
//    this->displayMyView(outImgC, None);
    outImgC.save(FILE_SPLINE);
}

void SimpleView::debugImage(uchar2D::Pointer inputImage) {
    uchar2D::RegionType region = inputImage->GetLargestPossibleRegion();
    uchar2D::SizeType size;
    size = region.GetSize();

    const int width = size[0];
    const int height = size[1];

    uchar2D::IndexType index;

    std::vector<unsigned char> valueList;
    for (int w = 0 ; w < width ; w++) {
        for (int h = 0 ; h < height ; h++) {
            index[0] = w;
            index[1] = h;
            const unsigned char value = inputImage->GetPixel(index);
            if (std::find(valueList.begin(), valueList.end(), value) == valueList.end()) {
                valueList.push_back(value);
            }
        }
    }

    std::sort(valueList.begin(), valueList.end());

    std::vector<unsigned char>::iterator it;
    for (it = valueList.begin() ; it < valueList.end() ; ++it) {
        cout << (int)*it << " ";
    }
    cout << endl << "size=" << valueList.size() << endl;
}
void SimpleView::drawSpline2()
{
    QImage inImg = QImage(FILE_SOBEL_RED);
    QImage outImg = QImage(InputFile.toLatin1().data());
    QImage outImgC = outImg.convertToFormat(QImage::Format_RGB888);
    std::vector<double> sp_x1, sp_y1, sp_x2, sp_y2;
    tk::spline s1, s2;
    int last_red_y1 = -1, last_red_x1 = -1;
    int last_red_y2 = -1, last_red_x2 = -1;
    float slope = 0.0, old_slope = 0.0;
    for(int x = 0; x < inImg.width(); x += 18)
    {
        for(int y = 0; y < inImg.height(); y++)
        {
            QRgb rgb = inImg.pixel(x, y);
            if(qRed(rgb) >=255)
            {
                slope = (float)(y-last_red_y1) / (float)(x-last_red_x1);
                cout<<"(x,y) = ("<<x<<","<<y<<")";
                cout<<", (last_red_x1, last_red_y1) = ("<<last_red_x1<<","<<last_red_y1<<")";
                cout<<", slope = "<<slope;
                cout<<", old_slope = "<<old_slope<<endl;

                if(last_red_y1 == -1 || //must sample first point
                   last_red_x1 == 0  || //must sample second point, let the old_slope be correct
                   fabs(slope-old_slope) < 0.55)
                {
                    sp_x1.push_back(x);
                    sp_y1.push_back(y);
                    drawColorDot(&inImg,qRgb(255,255,0), x, y);
                    last_red_x1 = x;
                    last_red_y1 = y;
                    old_slope = slope;
                }
                break;
            }
        }
        //find top line of the bottom edge
        int y2 = 0, y_use = 0, y_cnt = 0;
        bool every_white = false;
        for(int y = inImg.height()-1; y >= 0; y--)
        {
            QRgb rgb = inImg.pixel(x, y);
            if(qRed(rgb) >=255)
            {
                y_cnt ++;
                cout<<"(x,y) = ("<<x<<","<<y<<")";
                for (int i = y; i >=0; i-=2)
                {
                    QRgb t_rgb = inImg.pixel(x, i);
                    if (qGreen(t_rgb) > 0) every_white = true;
                    if(qRed(t_rgb) >=255 && every_white)
                    {
                        y_cnt ++;
                        if (y_cnt == 2) y2 = i; //remember the line2's Y cordinate
                        cout<<", (x,i) = ("<<x<<","<<i<<")";
                        cout<<", y_cnt = "<<y_cnt<<endl;
                    }
                }
                if (y_cnt >= 3) //means total >= 3 lines, but I needs line2
                {
                    y_use = y2;
                }
                else if (y_cnt == 2) //means total has 2 lines, but I needs line1
                {
                    y_use = y;
                }
                cout<<"(x,y_use) = ("<<x<<","<<y_use<<")";
                cout<<", y_cnt = "<<y_cnt<<endl;

                slope = (float)(y_use-last_red_y2) / (float)(x-last_red_x2);
                if(last_red_y2 == -1 ||
                   last_red_x2 == 0  ||
                   fabs(slope-old_slope) < 0.55)
                {
                    sp_x2.push_back(x);
                    sp_y2.push_back(y_use);
                    drawColorDot(&inImg,qRgb(255,0,255), x, y_use);
                    last_red_x2 = x;
                    last_red_y2 = y_use;
                    old_slope = slope;
                }
                break;
            }
        }
    }
    inImg.save(FILE_SPLINE_SAMPLE);
    s1.set_points(sp_x1,sp_y1);
    s2.set_points(sp_x2,sp_y2);
    memset(SplineY1, 0, sizeof(SplineY1));
    memset(SplineY2, 0, sizeof(SplineY2));
    lowestY_x1 = lowestY_x2 = 0;
    int max_y1 = 0, max_y2 = 0;
    for (int col = 0; col<outImgC.width(); col++)
    {
        outImgC.setPixel(col, (int)s1((double)col), qRgb(255, 0, 0));
        outImgC.setPixel(col, (int)s2((double)col), qRgb(255, 0, 0));
        //recording the 2 spline cordinates, and lowestY_x1, lowestY_x2
        SplineY1[col] = (int)s1((double)col);
        SplineY2[col] = (int)s2((double)col);
        if (SplineY1[col] > max_y1)
        {
            max_y1 = SplineY1[col];
            lowestY_x1 = col;
        }
        if (SplineY2[col] > max_y2)
        {
            max_y2 = SplineY2[col];
            lowestY_x2 = col;
        }
    }
    drawColorDot(&outImgC, qRgb(0,0,255),lowestY_x1, SplineY1[lowestY_x1]);
    drawColorDot(&outImgC, qRgb(0,255,0),lowestY_x2, SplineY2[lowestY_x2]);
    this->displayMyView(outImgC, None);
    outImgC.save(FILE_SPLINE);
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
//    callback->SetViewer2(this->ui->qvtkWidget_Ori);
    callback->SetActor(actor);

//    this->ui->qvtkWidget_Ori->GetRenderWindow()->GetInteractor()->Initialize();
//    this->ui->qvtkWidget_Ori->GetRenderWindow()->GetInteractor()->Start();
//    this->ui->qvtkWidget_Ori->SetRenderWindow(renderWindow);

    // window interactor style for display images
    vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
    // set interactor style to the qvtkWidget Interactor
//    this->ui->qvtkWidget_Ori->GetInteractor()->SetInteractorStyle(style);
//    this->ui->qvtkWidget_Ori->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::MouseMoveEvent, callback);
//    this->ui->qvtkWidget_Ori->update();
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
//    this->ui->qvtkWidget_Ori->SetRenderWindow(renderWindow);

    // window interactor style for display images
    vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
    // set interactor style to the qvtkWidget Interactor
//    this->ui->qvtkWidget_Ori->GetInteractor()->SetInteractorStyle(style);

//    this->ui->qvtkWidget_Ori->update();
#if 1   //test from chris, get mouse event and put a pixel on it
    // get double click events
    vtkCallbackCommand *callback = vtkCallbackCommand::New();
    callback->SetCallback(SimpleView::handle_double_click);
    callback->SetClientData(this);
//    this->ui->qvtkWidget_Ori->GetInteractor()->AddObserver(vtkCommand::LeftButtonPressEvent, callback, 1.0);
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
//    double *tmp_p = self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->GetPosition();
//    double *tmp_fp = self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->GetFocalPoint();
//    double p[3] = {tmp_p[0], tmp_p[1], tmp_p[2]};
//    double fp[3] = {tmp_fp[0], tmp_fp[1], tmp_fp[2]};
    // draw color point
//    self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actor);
    // reset camera
//    self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->ResetCamera();
//    self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->Zoom(1.5);
    // restore camera with original value
//    self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->SetPosition(p);
//    self->ui->qvtkWidget_Ori->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->SetFocalPoint(fp);
//    self->ui->qvtkWidget_Ori->GetRenderWindow()->Render();
//    self->ui->qvtkWidget_Ori->update();
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
