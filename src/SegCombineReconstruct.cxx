#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesWriter.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkMetaDataObject.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkBasicDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkCastImageFilter.h"
#include "four-neighbour.h"
#include "itkOffset.h"
#include "itkGDCMImageIO.h"

#include "itkImageToVTKImageFilter.h"

#include "itkExceptionObject.h"
#include <time.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMarchingCubes.h>
#include <vtkVoxelModeller.h>
#include <vtkSphereSource.h>
#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageReader2.h>
 
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include "vtkAutoInit.h"

VTK_MODULE_INIT(vtkRenderingOpenGL2); 
VTK_MODULE_INIT(vtkInteractionStyle);
typedef int16_t VolumePixelType;
typedef itk::Image<VolumePixelType, 3> VolumeImageType;
typedef itk::ImageRegionIterator<VolumeImageType> VolumeIteratorType;

inline bool IsROITag(int tag,const std::vector<int>&ROITag){
    bool flag = false;
    for(int i=0;i<ROITag.size();i++){
           if(ROITag[i] == tag){
            flag = true;
            break;
        }
    }
    return flag;
}

int SegmentSingleDCM(int sliceNumber,std::string inputFile, VolumeIteratorType& voliter){
    clock_t startTime,finishTime;
    startTime = clock();
    typedef double InputPixelType;
    typedef double OutputPixelType;
    typedef int16_t IntPixelType;
    typedef unsigned char ScalarPixelType;

    typedef itk::Image <ScalarPixelType, 2> ScalarImageType;
    typedef itk::Image <InputPixelType, 2>  InputImageType;
    typedef itk::Image <OutputPixelType, 2>  OutputImageType;
    typedef itk::Image<IntPixelType,2> DCMImageType;

    typedef itk::ImageFileReader <InputImageType>  ReaderType;
    //typedef itk::ImageSeriesReader <InputImageType> ReaderType;
    typedef itk::ImageFileReader <DCMImageType> DCMReaderType;
    //typedef itk::ImageSeriesReader <DCMImageType> DCMReaderType;

    typedef itk::ImageRegionIterator<InputImageType> InputIteratorType;
    typedef itk::ImageRegionIterator<DCMImageType> DCMIteratorType;

    typedef itk::ImageFileWriter <DCMImageType> DCMWriterType;
    //typedef itk::ImageSeriesWriter <DCMImageType> 
    typedef itk::ImageFileWriter <ScalarImageType> WriterType;
    
    typedef itk::RecursiveGaussianImageFilter <InputImageType, OutputImageType>  FilterType;
    
    ReaderType::Pointer reader = ReaderType::New();
    //argv[1] = "C:\\Users\\lyf\\Documents\\ITK_workspace\\CannyEdge_bin\\Debug\\1.2.840.113619.2.55.3.2797609686.440.1500426707.400.130.dcm";
    reader -> SetFileName(inputFile);
    typedef itk::GDCMImageIO           ImageIOType;
    ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
    //reader->SetImageIO( gdcmImageIO );
    

    try{
        reader->Update();
    }
    catch(itk::ExceptionObject &e){
        std::cerr<<"In reader ";
        std::cerr<<e.what()<<std::endl;
        return EXIT_FAILURE;
    }

    
    FilterType::Pointer filterX = FilterType::New();
    FilterType::Pointer filterY = FilterType::New();
    filterX -> SetDirection(0);
    filterY -> SetDirection(1);
    filterX -> SetOrder(FilterType::ZeroOrder);
    filterY -> SetOrder(FilterType::ZeroOrder);
    filterX -> SetNormalizeAcrossScale(false);
    filterY -> SetNormalizeAcrossScale(false);
    filterX -> SetInput(reader -> GetOutput());
    filterY -> SetInput(filterX -> GetOutput());
    //argv[3] = "0.1";
    //const double sigma = atof(argv[3]);
    const double sigma = 0.1;
    filterX -> SetSigma(sigma);
    filterY -> SetSigma(sigma);
    double gaussianTime;
    try{
        filterY -> Update();
        gaussianTime = clock();
        std::cout<<"Gaussian duration: "<<(double)(gaussianTime - startTime) / CLOCKS_PER_SEC<<'\n';
    }
    catch(itk::ExceptionObject &e){
        std::cerr<<"In Gaussian filterY ";
        std::cerr<<e.what()<<std::endl;
        return EXIT_FAILURE;
    }


    typedef itk::RescaleIntensityImageFilter <OutputImageType, ScalarImageType> RescaleFilterType;
    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler -> SetOutputMinimum(0);
    rescaler -> SetOutputMaximum(255);
    rescaler -> SetInput(filterY -> GetOutput());
    try{
        rescaler -> Update();
    }
    catch(itk::ExceptionObject &e){
        std::cerr<<"In rescaler ";
        std::cerr<<e.what()<<std::endl;
        return EXIT_FAILURE;
    }
    
    typedef itk::ScalarImageKmeansImageFilter<ScalarImageType> KMeansFilterType;
    KMeansFilterType::Pointer kmeansfilter = KMeansFilterType::New();
    kmeansfilter -> SetInput(rescaler -> GetOutput());

    //argv[4] = "1";
    //argv[5] = "3";
    //argv[6] = "14";
    //argv[7] = "90";
    //argv[8] = "130";
    //const int numberOfInitialClasses = atoi(argv[5]);
    const int numberOfInitialClasses = 3;
    //const int useNonContiguousLabels = atoi(argv[4]);
    const int useNonContiguousLabels = 1;
    const int argoffset = 6;

    kmeansfilter -> SetUseNonContiguousLabels(useNonContiguousLabels);

    //for( unsigned k=0; k < numberOfInitialClasses; k++ ){
    //    const double userProvidedInitialMean = atof(parameters[k+argoffset] );
    //    kmeansfilter -> AddClassWithInitialMean(userProvidedInitialMean);
    //}
	kmeansfilter -> AddClassWithInitialMean(14);
	kmeansfilter -> AddClassWithInitialMean(90);
	kmeansfilter -> AddClassWithInitialMean(130);
    double kmeansTime;
    try{
        kmeansfilter->Update();
        kmeansTime = clock();
        std::cout<<"Kmeans duration: "<<(double)(kmeansTime - gaussianTime)/CLOCKS_PER_SEC<<'\n';
    }
    catch(itk::ExceptionObject &e){
        std::cerr<<"In kmeans filter ";
        std::cerr<<e.what()<<std::endl;
        return EXIT_FAILURE;
    }
 
    typedef double RealPixelType;
    typedef itk::Image<RealPixelType, 2> RealImageType;
    typedef itk::CastImageFilter<ScalarImageType, RealImageType> CastToRealFilterType;
    //typedef itk::CastImageFilter<RealImageType, ScalarImageType> CastToCharFilterType;

    CastToRealFilterType::Pointer torealfilter = CastToRealFilterType::New();
    torealfilter -> SetInput(kmeansfilter -> GetOutput());
    
    float variance = 0.0;
    float upperThreshold = 0.0;
    float lowerThreshold = 0.0;
    //float maximumError = 0.1;

    typedef itk::CannyEdgeDetectionImageFilter<RealImageType, RealImageType> CannyFilterType;
    CannyFilterType::Pointer cannyfilter = CannyFilterType::New();
    cannyfilter -> SetInput(torealfilter -> GetOutput());
    cannyfilter -> SetVariance(variance);
    cannyfilter -> SetUpperThreshold(upperThreshold);
    cannyfilter -> SetLowerThreshold(lowerThreshold);
    double cannyedgeTime;
    try{
        cannyfilter -> Update();
        cannyedgeTime = clock();
        std::cout<<"Canny edge duration: "<<(double)(cannyedgeTime - kmeansTime)/CLOCKS_PER_SEC<<'\n';
    }
    catch(itk::ExceptionObject &e){
        std::cerr<<"In cannyfilter ";
        std::cerr<<e.what()<<std::endl;
        return EXIT_FAILURE;
    }
    
    //CastToCharFilterType::Pointer tocharfilter = CastToCharFilterType::New();
    //tocharfilter -> SetInput(cannyfilter -> GetOutput());
    rescaler -> SetInput(cannyfilter -> GetOutput());
    rescaler -> Update();

    std::vector<int> cstRun, cenRun, rRun, tag; 
    std::vector<std::pair<int,int>> equivalence;
    std::vector<std::pair<int,int>> sortedArea;
    recordRuns(rescaler -> GetOutput(),cstRun,cenRun,rRun);
    firstPass(cstRun, cenRun, rRun, tag, equivalence);
    secondPass(tag,equivalence);
    calculateArea(cstRun,cenRun,tag,sortedArea);

    //segment roi
    int thirdLargestTag = sortedArea[2].first;
    int fourthLargestTag = sortedArea[3].first;
    std::cout<<sliceNumber<<":\n"<<"thirdLargest:"<<thirdLargestTag<<' '<<sortedArea[2].second<<'\n'
             <<"fourthLargest:"<<fourthLargestTag<<' '<<sortedArea[3].second<<'\n';
    std::vector<int> ROITag;
    //for(int i=0;i<sortedArea.size();i++){
    //    if(sortedArea[i].second<1800&&sortedArea[i].second>100){
    //        ROITag.push_back(sortedArea[i].first);
    //    }
    //    else if(sortedArea[i].second<=100)
    //        break;
    //}
    inline bool IsROITag(int tag, const std::vector<int>&ROITag);


    DCMReaderType::Pointer reader1 = DCMReaderType::New();
    reader1 -> SetFileName(inputFile);
   // itk::GDCMImageIO::Pointer gdcmimageio = itk::GDCMImageIO::New();
    reader1 ->SetImageIO(gdcmImageIO);
    try{
        reader1 -> Update();
    }
    catch(itk::ExceptionObject &e){
        std::cerr<<"In reader1 ";
        std::cerr<<e.what()<<std::endl;
        return EXIT_FAILURE;
    }
    typedef itk::MetaDataDictionary   DictionaryType;
    DictionaryType & dictionary = reader1->GetOutput()->GetMetaDataDictionary();
    //Slice Location
    //failed
    //or maybe i should use 3d data as input and slice to 2d to process, translate index into physical position
    //https://itk.org/Wiki/ITK/Examples/DICOM/ResampleDICOM
    //std::ostringstream value;
    //value.str("");
    //double sliceLocation = -321-0.625*sliceNumber;//order and spacing need to be modified,just a key test
    //value<<sliceLocation;
    //itk::EncapsulateMetaData<std::string>( dictionary, "0020|1041", value.str() );

    
    
    DCMImageType *input = reader1 -> GetOutput();//do test on a copy
    DCMIteratorType iter(input, input -> GetRequestedRegion());
    iter.GoToBegin();

    int index = 0;
    for(int i=1;i<=512;i++){
        for(int j=1;j<=512;j++,voliter++,iter++){
            if(//IsROITag(tag[index],ROITag)&&
				(tag[index]==thirdLargestTag||tag[index]==fourthLargestTag)&&
                rRun[index] == i && j>=cstRun[index] && j<= cenRun[index]){
                if(iter.Get()>-2000&&iter.Get()<2000){
                    voliter.Set(iter.Get());
					
				}
                else{
                    voliter.Set(-2000);
					std::cout<<"Error: Out of range!"<<'\n';
				}
            }
            else{
                voliter.Set(-2000);
            }
			if(i==250&&j==250&&sliceNumber==10)
				std::cout<<"testPixel is set to:"<<voliter.Get()<<'\n';
            if(rRun[index] == i && cenRun[index] == j){
                index +=1;
            }
        }
    }
    //myresult->Update();
    //DCMWriterType::Pointer writer1 = DCMWriterType::New();
    //writer1 -> SetFileName("C:\\Users\\lyf\\Documents\\ITK_workspace\\CannyEdge_bin\\Debug\\segment.3rd4thfull.my.dcm");
    //char *temp = strcat(outputFolder,'\\');
    //char *outputFile = strcat(temp,(char*)(inputFile));

	//writer1 -> SetFileName(std::string(outputFolder).append("\\").append(inputFile,inputFile.length()-7,7));
    //writer1 -> SetInput(myresult);
    //writer1 -> UseInputMetaDataDictionaryOff ();
    //writer1 -> SetImageIO(gdcmImageIO);
    //try{
    //    writer1 -> Update();
    //}
    //catch(itk::ExceptionObject &ex)
    //{
    //    std::cerr<<ex.what()<<std::endl;
	//	return EXIT_FAILURE;
    //}
    finishTime = clock();
    double duration = (double)(finishTime - cannyedgeTime) / CLOCKS_PER_SEC;
    std::cout<<"Connected region duration:"<<duration<<std::endl;
	return EXIT_SUCCESS;
}

//parameters:input_folder isovalue
//input: a folder containing CT images
//output: vtk reconstruction result
int main(int argc, char *argv[])
{
	argv[1] = "C:\\SRT\\testCT\\LungTest";
	argv[2] = "-900";
    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
    namesGenerator -> SetInputDirectory(argv[1]);
	try{
        namesGenerator -> Update();
	}
	catch(itk::ExceptionObject &ex){
		std::cout<<ex.what()<<std::endl;
		return EXIT_FAILURE;
	}
    std::vector<std::string> inputFileNames = namesGenerator->GetInputFileNames();

    
    VolumeImageType::Pointer reconstruct = VolumeImageType::New();
    VolumeImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    VolumeImageType::SizeType size;
    size[0] = 512;
    size[1] = 512;
    size[2] = inputFileNames.size();
    VolumeImageType::RegionType region;
    region.SetIndex(start);
    region.SetSize(size);
    reconstruct->SetRegions(region);
    reconstruct->Allocate();
    VolumeIteratorType voliter(reconstruct,reconstruct->GetRequestedRegion());
    voliter.GoToBegin();
    for(int i=0;i<size[2];i++){
        SegmentSingleDCM(i,inputFileNames[i],voliter);
		reconstruct->Update();
		std::cout<<"Slice "<<i<<"succeed!"<<'\n';
    }
    reconstruct->Update(); 

    //vtk marching cube
    typedef itk::ImageToVTKImageFilter< VolumeImageType > itktovtkFilterType;
    itktovtkFilterType::Pointer itktovtkFilter = itktovtkFilterType::New();
    itktovtkFilter->SetInput(reconstruct);
    try{
        itktovtkFilter->Update();
    }
    catch(itk::ExceptionObject &e){
        std::cerr<<"In itk to vtk image data "<<e.what();
        return EXIT_FAILURE;
    }


    vtkSmartPointer<vtkImageData> volume = vtkSmartPointer<vtkImageData>::New();
    double isoValue;

    isoValue = atof(argv[2]);
    volume->DeepCopy(itktovtkFilter->GetOutput());

	vtkSmartPointer<vtkMarchingCubes> surface = vtkSmartPointer<vtkMarchingCubes>::New();
    #if VTK_MAJOR_VERSION <= 5
    surface->SetInput(volume);
    #else
    surface->SetInputData(volume);
    #endif
    surface->ComputeNormalsOn();
    surface->SetValue(0, isoValue);
 
    vtkSmartPointer<vtkRenderer> renderer =  vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(.1, .2, .3);
 
    vtkSmartPointer<vtkRenderWindow> renderWindow =  vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> interactor =  vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);
 
    vtkSmartPointer<vtkPolyDataMapper> mapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(surface->GetOutputPort());
    mapper->ScalarVisibilityOff();
 
    vtkSmartPointer<vtkActor> actor =  vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
 
    renderer->AddActor(actor);
 
    renderWindow->Render();
    interactor->Start();
  
    return EXIT_SUCCESS;
}