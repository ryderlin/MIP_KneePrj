#-------------------------------------------------
#
# Project created by QtCreator 2015-10-19T16:13:56
#
#-------------------------------------------------

QT       += core

#QT       -= gui

QT       += widgets
TARGET = Knee_Project
#CONFIG   += console
#CONFIG   -= app_bundle
#QMAKE_CXXFLAGS += -std=c++11
TEMPLATE = app

#-------------------------------------------------
# include path
#-------------------------------------------------
INCLUDEPATH += /home/ryderlin/VTK/install_6.3.0/include/vtk-6.3 \
               /home/ryderlin/ITK/install_4.8.0/include/ITK-4.8 \
               ./

#-------------------------------------------------
# lib path
#-------------------------------------------------
LIBS += -L/home/ryderlin/VTK/install_6.3.0/lib \ # VTK
        -lvtkGUISupportQt-6.3 \
        -lvtkIOImage-6.3 \
        -lvtkInteractionImage-6.3 \
        -lvtkRenderingCore-6.3 \
        -lvtkCommonExecutionModel-6.3 \
        -lvtkCommonCore-6.3 \
        -lvtkRenderingOpenGL-6.3 \
        -lvtkalglib-6.3 \
        -lvtkChartsCore-6.3 \
        -lvtkCommonColor-6.3 \
        -lvtkCommonDataModel-6.3 \
        -lvtkCommonMath-6.3 \
        -lvtkCommonCore-6.3 \
        -lvtksys-6.3 \
        -lvtkCommonMisc-6.3 \
        -lvtkCommonSystem-6.3 \
        -lvtkCommonTransforms-6.3 \
        -lvtkInfovisCore-6.3 \
        -lvtkFiltersExtraction-6.3 \
        -lvtkCommonExecutionModel-6.3 \
        -lvtkFiltersCore-6.3 \
        -lvtkFiltersGeneral-6.3 \
        -lvtkCommonComputationalGeometry-6.3 \
        -lvtkFiltersStatistics-6.3 \
        -lvtkImagingFourier-6.3 \
        -lvtkImagingCore-6.3 \
        -lvtkRenderingContext2D-6.3 \
        -lvtkRenderingCore-6.3 \
        -lvtkFiltersGeometry-6.3 \
        -lvtkFiltersSources-6.3 \
        -lvtkRenderingFreeType-6.3 \
        -lvtkfreetype-6.3 \
        -lvtkzlib-6.3 \
        -lvtkftgl-6.3 \
        -lvtkRenderingOpenGL-6.3 \
        -lvtkImagingHybrid-6.3 \
        -lvtkIOImage-6.3 \
        -lvtkDICOMParser-6.3 \
        -lvtkIOCore-6.3 \
        -lvtkmetaio-6.3 \
        -lvtkjpeg-6.3 \
        -lvtkpng-6.3 \
        -lvtktiff-6.3 \
        -lvtkDomainsChemistry-6.3 \
        -lvtkIOXML-6.3 \
        -lvtkIOGeometry-6.3 \
        -lvtkjsoncpp-6.3 \
        -lvtkIOXMLParser-6.3 \
        -lvtkexpat-6.3 \
        -lvtkexoIIc-6.3 \
        -lvtkNetCDF-6.3 \
        -lvtkNetCDF_cxx-6.3 \
        -lvtkFiltersAMR-6.3 \
        -lvtkParallelCore-6.3 \
        -lvtkIOLegacy-6.3 \
        -lvtkFiltersFlowPaths-6.3 \
        -lvtkFiltersGeneric-6.3 \
        -lvtkFiltersHybrid-6.3 \
        -lvtkImagingSources-6.3 \
        -lvtkFiltersHyperTree-6.3 \
        -lvtkFiltersImaging-6.3 \
        -lvtkImagingGeneral-6.3 \
        -lvtkFiltersModeling-6.3 \
        -lvtkFiltersParallel-6.3 \
        -lvtkFiltersParallelImaging-6.3 \
        -lvtkFiltersProgrammable-6.3 \
        -lvtkFiltersSelection-6.3 \
        -lvtkFiltersSMP-6.3 \
        -lvtkFiltersTexture-6.3 \
        -lvtkFiltersVerdict-6.3 \
        -lvtkverdict-6.3 \
        -lvtkGeovisCore-6.3 \
        -lvtkInfovisLayout-6.3 \
        -lvtkInteractionStyle-6.3 \
        -lvtkInteractionWidgets-6.3 \
        -lvtkRenderingAnnotation-6.3 \
        -lvtkImagingColor-6.3 \
        -lvtkRenderingVolume-6.3 \
        -lvtkViewsCore-6.3 \
        -lvtkproj4-6.3 \
        -lvtkgl2ps-6.3 \
        -lvtkGUISupportQt-6.3 \
        -lvtkGUISupportQtOpenGL-6.3 \
        -lvtkGUISupportQtSQL-6.3 \
        -lvtkIOSQL-6.3 \
        -lvtksqlite-6.3 \
        -lvtkGUISupportQtWebkit-6.3 \
        -lvtkViewsQt-6.3 \
        -lvtkViewsInfovis-6.3 \
        -lvtkRenderingLabel-6.3 \
        -lvtkImagingMath-6.3 \
        -lvtkImagingMorphological-6.3 \
        -lvtkImagingStatistics-6.3 \
        -lvtkImagingStencil-6.3 \
        -lvtkInteractionImage-6.3 \
        -lvtkIOAMR-6.3 \
        -lvtkIOEnSight-6.3 \
        -lvtkIOExodus-6.3 \
        -lvtkIOExport-6.3 \
        -lvtkRenderingGL2PS-6.3 \
        -lvtkIOImport-6.3 \
        -lvtkIOInfovis-6.3 \
        -lvtklibxml2-6.3 \
        -lvtkIOLSDyna-6.3 \
        -lvtkIOMINC-6.3 \
        -lvtkIOMovie-6.3 \
        -lvtkoggtheora-6.3 \
        -lvtkIONetCDF-6.3 \
        -lvtkIOParallel-6.3 \
        -lvtkIOPLY-6.3 \
        -lvtkIOVideo-6.3 \
#        -lvtkRenderingFreeTypeOpenGL-6.3 \
        -lvtkRenderingImage-6.3 \
        -lvtkRenderingLIC-6.3 \
        -lvtkRenderingLOD-6.3 \
        -lvtkRenderingQt-6.3 \
#        -lvtkRenderingVolumeAMR-6.3 \
        -lvtkRenderingVolumeOpenGL-6.3 \
        -lvtkViewsContext2D-6.3 \
#        -lvtkViewsGeovis-6.3 \
        -lvtkgl2ps-6.3 \
        -lvtkexoIIc-6.3 \
        -lvtkFiltersParallel-6.3 \
        -lvtkIONetCDF-6.3 \
        -lvtkNetCDF_cxx-6.3 \
        -lvtkNetCDF-6.3 \
        -lvtkFiltersTexture-6.3 \
        -lvtkGUISupportQt-6.3 \
        -lvtkFiltersAMR-6.3 \
        -lvtkParallelCore-6.3 \
        -lvtkIOLegacy-6.3 \
        -lvtkGeovisCore-6.3 \
        -lvtkIOXML-6.3 \
        -lvtkIOGeometry-6.3 \
        -lvtkjsoncpp-6.3 \
        -lvtkIOXMLParser-6.3 \
        -lvtkexpat-6.3 \
        -lvtkproj4-6.3 \
        -lvtkViewsInfovis-6.3 \
        -lvtkChartsCore-6.3 \
        -lvtkCommonColor-6.3 \
        -lvtkRenderingContext2D-6.3 \
        -lvtkRenderingOpenGL-6.3 \
        -lvtkFiltersImaging-6.3 \
        -lvtkInfovisLayout-6.3 \
        -lvtkInfovisCore-6.3 \
        -lvtkViewsCore-6.3 \
        -lvtkInteractionWidgets-6.3 \
        -lvtkImagingHybrid-6.3 \
        -lvtkIOImage-6.3 \
        -lvtkDICOMParser-6.3 \
        -lvtkIOCore-6.3 \
        -lvtkmetaio-6.3 \
        -lvtkpng-6.3 \
        -lvtktiff-6.3 \
        -lvtkjpeg-6.3 \
        -lvtkFiltersHybrid-6.3 \
        -lvtkImagingGeneral-6.3 \
        -lvtkImagingSources-6.3 \
        -lvtkFiltersModeling-6.3 \
        -lvtkInteractionStyle-6.3 \
        -lvtkRenderingAnnotation-6.3 \
        -lvtkImagingColor-6.3 \
        -lvtkRenderingVolume-6.3 \
        -lvtkRenderingLabel-6.3 \
        -lvtkRenderingFreeType-6.3 \
        -lvtkRenderingCore-6.3 \
        -lvtkFiltersExtraction-6.3 \
        -lvtkFiltersStatistics-6.3 \
        -lvtkalglib-6.3 \
        -lvtkImagingFourier-6.3 \
        -lvtkImagingCore-6.3 \
        -lvtkFiltersGeometry-6.3 \
        -lvtkFiltersSources-6.3 \
        -lvtkFiltersGeneral-6.3 \
        -lvtkFiltersCore-6.3 \
        -lvtkCommonExecutionModel-6.3 \
        -lvtkCommonComputationalGeometry-6.3 \
        -lvtkCommonDataModel-6.3 \
        -lvtkCommonMisc-6.3 \
        -lvtkCommonTransforms-6.3 \
        -lvtkCommonMath-6.3 \
        -lvtkCommonSystem-6.3 \
        -lvtkCommonCore-6.3 \
        -lvtksys-6.3 \
        -lvtkftgl-6.3 \
        -lvtkfreetype-6.3 \
        -lvtkzlib-6.3 \
        -L/home/ryderlin/ITK/install_4.8.0/lib \ # ITK
        -litkdouble-conversion-4.8 \
        -litksys-4.8 \
        -litkvnl_algo-4.8 \
        -litkvnl-4.8 \
        -lITKCommon-4.8 \
        -lITKStatistics-4.8 \
        -lITKIOImageBase-4.8 \
        -lITKMesh-4.8 \
        -lITKMetaIO-4.8 \
        -lITKSpatialObjects-4.8 \
        -lITKPath-4.8 \
        -lITKLabelMap-4.8 \
        -lITKQuadEdgeMesh-4.8 \
        -lITKOptimizers-4.8 \
        -lITKPolynomials-4.8 \
        -lITKBiasCorrection-4.8 \
        -lITKBioCell-4.8 \
        -lITKDICOMParser-4.8 \
        -lITKEXPAT-4.8 \
        -lITKIOXML-4.8 \
        -lITKIOSpatialObjects-4.8 \
        -lITKFEM-4.8 \
        -litkopenjpeg-4.8 \
        -litkgdcmDICT-4.8 \
        -litkgdcmMSFF-4.8 \
        -lITKniftiio-4.8 \
        -lITKgiftiio-4.8 \
        -litkhdf5_cpp-4.8 \
        -litkhdf5-4.8 \
        -lITKIOBMP-4.8 \
        -lITKIOBioRad-4.8 \
        -lITKIOCSV-4.8 \
        -lITKIOGDCM-4.8 \
        -lITKIOIPL-4.8 \
        -lITKIOGE-4.8 \
        -lITKIOGIPL-4.8 \
        -lITKIOHDF5-4.8 \
        -litkjpeg-4.8 \
        -lITKIOJPEG-4.8 \
        -litktiff-4.8 \
        -lITKIOTIFF-4.8 \
        -lITKIOLSM-4.8 \
        -lITKIOMRC-4.8 \
        -lITKIOMesh-4.8 \
        -lITKIOMeta-4.8 \
        -lITKIONIFTI-4.8 \
        -lITKNrrdIO-4.8 \
        -lITKIONRRD-4.8 \
        -litkpng-4.8 \
        -lITKIOPNG-4.8 \
        -lITKIOSiemens-4.8 \
        -lITKIOStimulate-4.8 \
        -lITKIOTransformBase-4.8 \
        -lITKIOTransformHDF5-4.8 \
        -lITKIOTransformInsightLegacy-4.8 \
        -lITKIOTransformMatlab-4.8 \
        -lITKIOVTK-4.8 \
        -lITKKLMRegionGrowing-4.8 \
        -lITKOptimizersv4-4.8 \
        -lITKVTK-4.8 \
        -lITKVideoCore-4.8 \
        -lITKVideoIO-4.8 \
        -lITKWatersheds-4.8 \
        -lITKPolynomials-4.8 \
        -lITKIOXML-4.8 \
        -litkgdcmMSFF-4.8 \
        -litkopenjpeg-4.8 \
        -litkgdcmDICT-4.8 \
        -litkgdcmIOD-4.8 \
        -litkgdcmDSED-4.8 \
        -litkgdcmCommon-4.8 \
        -lITKIOTIFF-4.8 \
        -litktiff-4.8 \
        -litkjpeg-4.8 \
        -lITKgiftiio-4.8 \
        -lITKEXPAT-4.8 \
        -lITKMetaIO-4.8 \
        -lITKniftiio-4.8 \
        -lITKznz-4.8 \
        -lITKNrrdIO-4.8 \
        -litkpng-4.8 \
        -lITKIOGE-4.8 \
        -lITKIOIPL-4.8 \
        -litkhdf5_cpp-4.8 \
        -litkhdf5-4.8 \
        -lITKIOTransformBase-4.8 \
        -lITKOptimizers-4.8 \
        -lITKVideoCore-4.8 \
        -lITKSpatialObjects-4.8 \
        -lITKIOImageBase-4.8 \
        -lITKMesh-4.8 \
        -lITKPath-4.8 \
        -lITKStatistics-4.8 \
        -lITKCommon-4.8 \
        -litkdouble-conversion-4.8 \
        -litksys-4.8 \
        -lITKVNLInstantiation-4.8 \
        -litkvnl_algo-4.8 \
        -litkv3p_lsqr-4.8 \
        -litkvnl-4.8 \
        -litkvcl-4.8

SOURCES += \
    main.cxx \
    SimpleView.cxx \
    C_fileIO.cxx \
    C_itkSeg.cxx \
    C_tool.cxx

HEADERS += \
#    itkInclude.h \
#    itkTypeDef.h \
    SimpleView.h \
    C_fileIO.h \
    C_itkSeg.h \
    C_tool.h \
    spline.h \
    itkImageToVTKImageFilter.h

FORMS += \
    SimpleView.ui

RESOURCES += \
    Icons/icons.qrc

DISTFILES += \
    CMakeLists.txt.user \
    Knee_Project.pro.user \
    CMakeLists.txt
