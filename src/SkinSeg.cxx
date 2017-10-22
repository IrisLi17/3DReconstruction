#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesWriter.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkGrayscaleFillholeImageFilter.h"
#include "itkExceptionObject.h"
#include "itkImageToVTKImageFilter.h"
#include <fstream>
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
//VTK_MODULE_INIT(vtkRenderingOpenGL2); 
//VTK_MODULE_INIT(vtkInteractionStyle);
typedef int16_t DCMPixelType;
typedef unsigned char CharPixelType;
typedef itk::Image<DCMPixelType, 2> DCMImageType;
typedef itk::Image<CharPixelType, 2> CharImageType;
typedef itk::ImageFileReader <DCMImageType>  DCMReaderType;
typedef itk::ImageFileWriter <CharImageType> CharWriterType;
typedef itk::ImageFileWriter <DCMImageType> DCMWriterType;
typedef itk::ImageRegionIterator<DCMImageType> DCMIteratorType;
typedef itk::ImageRegionIterator<CharImageType> CharIteratorType;

typedef int16_t VolumePixelType;
typedef itk::Image<VolumePixelType, 3> VolumeImageType;
typedef itk::ImageRegionIterator<VolumeImageType> VolumeIteratorType;

int main(int argc, char *argv[]){
    ofstream diary("diary.txt",std::ios::out);
    void SingleSlice(int sliceNumber,std::string inputFile, VolumeIteratorType& voliter,ofstream& diary);
    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
    argv[1] = "C:\\SRT\\testCT\\LungTest";
	argv[2] = "0";
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
	clock_t t1=clock(),t2;
    for(int i=0;i<size[2];i++){
        SingleSlice(i,inputFileNames[i],voliter,diary);
		t2 = clock();
		diary<<inputFileNames[i]<<'\t'<<(t2-t1)/CLOCKS_PER_SEC<<'\n';
		t1 = t2;
        reconstruct->Update();
        std::cout<<"Slice "<<i<<" succeed!"<<'\n';
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
void SingleSlice(int sliceNumber,std::string inputFile, VolumeIteratorType& voliter,ofstream& diary){

    DCMReaderType::Pointer reader = DCMReaderType::New();
    CharWriterType::Pointer writer = CharWriterType::New();
	DCMWriterType::Pointer dcmwriter = DCMWriterType::New();
    
    //char *inputfile = inputFile;
    reader->SetFileName(inputFile);
    try{
        reader->Update();
    }
    catch(itk::ExceptionObject &e){
        diary<<e.what()<<'\n';
        std::cerr<<"Error"<<'\n';
    }
    //DCMIteratorType testiter(reader->GetOutput(),reader->GetOutput()->GetRequestedRegion());
    //testiter.GoToBegin();
    int count = 0;
	//int min = 0,max = 0;
    //for(count=0;count<512*512;count++){
    //    outfile<<testiter.Get()<<'\n';
	//	min = (testiter.Get()<min?testiter.Get():min);
	//	max = (testiter.Get()>max?testiter.Get():max);
    //    testiter++;
    //}
	//std::cout<<min<<' '<<max;
    //BinaryThreshold
    //not good
	itk::BinaryThresholdImageFilter<DCMImageType,CharImageType>::Pointer binarythresholdfilter
		= itk::BinaryThresholdImageFilter<DCMImageType,CharImageType>::New();
    binarythresholdfilter->SetInput(reader->GetOutput());
    binarythresholdfilter->SetLowerThreshold(-250);
    binarythresholdfilter->SetUpperThreshold(2000);
    binarythresholdfilter->SetInsideValue(255);
    binarythresholdfilter->SetOutsideValue(0);
    try{
        binarythresholdfilter->Update();
    }
    catch(itk::ExceptionObject &e){
        diary<<e.what()<<'\n';
        std::cerr<<"Error"<<'\n';
    }

    //writer->SetInput(binarythresholdfilter->GetOutput());
    //writer->SetFileName("temp.png");
    //writer->Update();

    //fillhole
    typedef itk::GrayscaleFillholeImageFilter<CharImageType, CharImageType> FillHoleFilterType;
    FillHoleFilterType::Pointer fillholefilter = FillHoleFilterType::New();
    fillholefilter -> SetInput(binarythresholdfilter -> GetOutput());
    try{
        fillholefilter -> Update();
    }
    catch(itk::ExceptionObject &e){
        diary<<e.what()<<'\n';
        std::cerr<<"Error"<<'\n';
    }

    //writer->SetInput(fillholefilter->GetOutput());
    //writer->SetFileName("temp_fillhole.png");
    //writer->Update();

    //expand minus erode
    typedef itk::BinaryBallStructuringElement<CharPixelType, 2> StructuringElementType;  
    StructuringElementType structuringElement;  
    structuringElement.SetRadius(3);  
    structuringElement.CreateStructuringElement();  
  
    typedef itk::BinaryMorphologicalOpeningImageFilter <CharImageType, CharImageType, StructuringElementType> BinaryMorphologicalOpeningImageFilterType;  
    BinaryMorphologicalOpeningImageFilterType::Pointer mathopenfilter  
            = BinaryMorphologicalOpeningImageFilterType::New();  
    mathopenfilter->SetInput(fillholefilter->GetOutput());  
    mathopenfilter->SetKernel(structuringElement);
    try{
        mathopenfilter->Update();
    }  
    catch(itk::ExceptionObject &e){
        diary<<e.what()<<'\n';
        std::cerr<<"Error"<<'\n';
    }

    typedef itk::BinaryErodeImageFilter <CharImageType, CharImageType, StructuringElementType>  
            BinaryErodeImageFilterType;  
    BinaryErodeImageFilterType::Pointer matherodefilter = BinaryErodeImageFilterType::New();  
    matherodefilter->SetInput(fillholefilter->GetOutput());  
    matherodefilter->SetKernel(structuringElement);
    try{
        matherodefilter->Update();
    }
    catch(itk::ExceptionObject &e){
        diary<<e.what()<<'\n';
        std::cerr<<"Error"<<'\n';
    }

    typedef itk::SubtractImageFilter<CharImageType, CharImageType> SubtractImageFilterType;
    SubtractImageFilterType::Pointer mathsubtractfilter = SubtractImageFilterType::New();
    mathsubtractfilter->SetInput1(mathopenfilter->GetOutput());
    mathsubtractfilter->SetInput2(matherodefilter->GetOutput());
    try{
        mathsubtractfilter->Update();
    }
    catch(itk::ExceptionObject &e){
        diary<<e.what()<<'\n';
        std::cerr<<"Error"<<'\n';
    }

    //output
    writer->SetInput(mathsubtractfilter->GetOutput());
    writer->SetFileName("minus3.png");
    try{
        writer->Update();
    }
    catch(itk::ExceptionObject &e){
        diary<<e.what()<<'\n';
        std::cerr<<"Error"<<'\n';
    }

    DCMImageType *input = reader->GetOutput();
    DCMIteratorType inputiter(input,input->GetRequestedRegion());
    inputiter.GoToBegin();

    CharImageType *bin = mathsubtractfilter->GetOutput();
    CharIteratorType biniter(bin,bin->GetRequestedRegion());
    biniter.GoToBegin();

    //DCMImageType::Pointer output = DCMImageType::New();
    //DCMImageType::IndexType start;
    //start[0] = 0;
    //start[1] = 0;
    //DCMImageType::SizeType size;
    //size[0] = 512;
    //size[1] = 512;
    //DCMImageType::RegionType region;
    //region.SetIndex(start);
    //region.SetSize(size);
    //output->SetRegions(region);
    //output->Allocate();
    //DCMIteratorType outputiter(output,output->GetRequestedRegion());
    //outputiter.GoToBegin();

	count=0;
    while(count<512*512){
        if(biniter.Get()==255){
            voliter.Set(inputiter.Get());
        }
        else{
            voliter.Set(-3024);
        }
		count++;
        biniter++;
        inputiter++;
        voliter++;
    }

    //dcmwriter->SetInput(output);
    //dcmwriter->SetFileName("SkinSeg.dcm");
    //try{
    //    dcmwriter->Update();
    //}
    //catch(itk::ExceptionObject &e){
    //    std::cerr<<e.what()<<'\n';
    //    return EXIT_FAILURE;
    //}

}

