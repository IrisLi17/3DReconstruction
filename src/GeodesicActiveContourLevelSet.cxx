/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkScaleTransform.h"
#include "itkResampleImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include <io.h>

const unsigned int                Dimension = 2;
typedef float                                    InputPixelType;
typedef itk::Image< InputPixelType, Dimension >  InputImageType;
typedef unsigned char                            OutputPixelType;
typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

typedef  itk::ImageFileReader< InputImageType >  ReaderType;
typedef  itk::ImageFileWriter< OutputImageType > WriterType;

typedef  itk::FastMarchingImageFilter< InputImageType, InputImageType > FastMarchingFilterType;
typedef FastMarchingFilterType::NodeContainer  NodeContainer;
typedef FastMarchingFilterType::NodeType       NodeType;

void getFiles( std::string path, std::vector<std::string>& files )  
{  
    //文件句柄  
    long   hFile   =   0;  
    //文件信息  
    struct _finddata_t fileinfo;  
    std::string p;  
    if((hFile = _findfirst(p.assign(path).append("\\*").c_str(),&fileinfo)) !=  -1)  
    {  
        do  
        {  
            //如果是目录,迭代之  
            //如果不是,加入列表  
            if((fileinfo.attrib &  _A_SUBDIR))  
            {  
                if(strcmp(fileinfo.name,".") != 0  &&  strcmp(fileinfo.name,"..") != 0)  
                    getFiles( p.assign(path).append("\\").append(fileinfo.name), files );  
            }  
            else  
            {  
                files.push_back(p.assign(path).append("\\").append(fileinfo.name) );  
            }  
        }while(_findnext(hFile, &fileinfo)  == 0);  
        _findclose(hFile);  
    }  
}

void updateSeeds(OutputImageType *contour,NodeContainer::Pointer seeds){
	itk::ImageRegionIterator<OutputImageType> it(contour,contour->GetRequestedRegion());
	it.GoToBegin();
	while(!it.IsAtEnd()){
		std::cout<<it.Value()<<'\n';
		it ++;
	}
}

bool cmp(std::string fn1,std::string fn2){
	std::string t1=fn1.substr(0,fn1.find_last_of("."));
	t1 = t1.substr(t1.find_last_of(".")+1);
	t1.reserve();
	int n1 = 0;
	for(int i=0;i<t1.length();i++){
		n1 *= 10;
		n1 += int(t1[i]-'0');
	}
	std::string t2=fn2.substr(0,fn2.find_last_of("."));
	t2 = t2.substr(t2.find_last_of(".")+1);
	t2.reserve();
	int n2 = 0;
	for(int i=0;i<t2.length();i++){
		n2 *= 10;
		n2 += int(t2[i]-'0');
	}
	return n1<n2;
}

int main( int argc, char* argv[] )
{
/*
  if( argc != 11 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " <InputFileName>  <OutputFileName>";
    std::cerr << " <seedX> <seedY> <InitialDistance>";
    std::cerr << " <Sigma> <SigmoidAlpha> <SigmoidBeta>";
    std::cerr << " <PropagationScaling> <NumberOfIterations>"  << std::endl;
    return EXIT_FAILURE;
    }

  /*const char * inputFileName =      argv[1];
  const char * outputFileName =     argv[2];
  const int seedPosX =              atoi( argv[3] );
  const int seedPosY =              atoi( argv[4] );

  const double initialDistance =    atof( argv[5] );
  const double sigma =              atof( argv[6] );
  const double alpha =              atof( argv[7] );
  const double beta  =              atof( argv[8] );
  const double propagationScaling = atof( argv[9] );
  const double numberOfIterations = atoi( argv[10] );
  const double seedValue =          - initialDistance;*/
  
  //const char * outputFileName =     "output_test.png";
  int seedPosX =              atoi( "360" );
  int seedPosY =              atoi( "310" );
  int initialSliceIndex =	  atoi( "184" );

  const double initialDistance =    atof( "14" );
  const double sigma =              atof( "1.0" );
  const double alpha =              atof( "-0.1" );
  const double beta  =              atof( "50.0" );
  const double propagationScaling = atof( "2.0" );
  const double numberOfIterations = atoi( "1000" );
  double seedValue =          - initialDistance;
  
  std::string path = "C:\\SRT\\testCT\\LungTest\\shenmingmou";
  std::vector<std::string> inputFiles;
  getFiles(path,inputFiles);
  std::sort(inputFiles.begin(),inputFiles.end(),cmp);
  
  FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
  

  InputImageType::IndexType  seedPosition1, seedPosition2;
  seedPosition1[0] = seedPosX;
  seedPosition1[1] = seedPosY;
  //seedPosition2[0] = 387;
  //seedPosition2[1] = 327;

  NodeContainer::Pointer seeds = NodeContainer::New();
  NodeType node1,node2;
  node1.SetValue( seedValue );
  node1.SetIndex( seedPosition1 );
  //node2.SetValue( 1.4 );
  //node2.SetIndex( seedPosition2);

    seeds->Initialize();
    seeds->InsertElement( 0, node1 );
  //seeds->SetElement(0, node1);
  //seeds->InsertElement( 1, node2 );

  
    fastMarching->SetSpeedConstant( 1.0 );
  
    typedef  itk::GeodesicActiveContourLevelSetImageFilter< InputImageType, InputImageType >  GeodesicActiveContourFilterType;
    GeodesicActiveContourFilterType::Pointer geodesicActiveContour = GeodesicActiveContourFilterType::New();
    geodesicActiveContour->SetPropagationScaling( propagationScaling );
    geodesicActiveContour->SetCurvatureScaling( 1.0 );
    geodesicActiveContour->SetAdvectionScaling( 1.0 );
    geodesicActiveContour->SetMaximumRMSError( 0.01 );
    geodesicActiveContour->SetNumberOfIterations( numberOfIterations );

	bool reverse = 1;
    for(int i=0;i<inputFiles.size();i++){
        ReaderType::Pointer reader = ReaderType::New();		
		reader->SetFileName( inputFiles[(2*reverse-1)*i+initialSliceIndex-1] );
        reader->Update();
	   //InputImageType *out = InputImageType::New();
	  //itk::ImageSource<InputImageType>::OutputImageType *feature = FeatureImage(inputFiles[i],reader->GetOutput(),out,sigma,alpha,beta);
	    typedef  itk::CurvatureAnisotropicDiffusionImageFilter< InputImageType, InputImageType > SmoothingFilterType;
	    SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();
        smoothing->SetTimeStep( 0.05 );
        smoothing->SetNumberOfIterations( 5 );
        smoothing->SetConductanceParameter( 9.0 );
	    smoothing->SetInput( reader->GetOutput() );

      typedef  itk::GradientMagnitudeRecursiveGaussianImageFilter< InputImageType, InputImageType > GradientFilterType;
      GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
      gradientMagnitude->SetSigma( sigma );
      gradientMagnitude->SetInput( smoothing->GetOutput() );

      typedef  itk::SigmoidImageFilter< InputImageType, InputImageType > SigmoidFilterType;
      SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
      sigmoid->SetOutputMinimum( 0.0 );
      sigmoid->SetOutputMaximum( 1.0 );
      sigmoid->SetAlpha( alpha );
      sigmoid->SetBeta( beta );
      sigmoid->SetInput( gradientMagnitude->GetOutput() );
	/*
      typedef itk::RescaleIntensityImageFilter<InputImageType, OutputImageType> CastFilter;
      CastFilter::Pointer caster = CastFilter::New();
      caster->SetInput(sigmoid->GetOutput());
      WriterType::Pointer writer = WriterType::New();
      //std::string _t = filename.substr(0,filename.find_last_of("."));
      //std::string index = _t.substr(_t.find_last_of(".")+1);
      writer->SetFileName("feature.png");
      writer->SetInput(caster->GetOutput());

      try{
	      writer->Update();
      }
      catch(itk::ExceptionObject &e){
	      std::cerr<<e.what()<<'\n';
	      return NULL;
      }
	  */
      try{
	      sigmoid->Update();
      }
      catch(itk::ExceptionObject &e){
	      std::cerr<<e.what();
	      return NULL;
      }

	  fastMarching->SetTrialPoints( seeds );
	  fastMarching->SetOutputSize(sigmoid->GetOutput()->GetLargestPossibleRegion().GetSize());
	  fastMarching->SetOutputOrigin(sigmoid->GetOutput()->GetOrigin());
	  fastMarching->SetOutputSpacing(sigmoid->GetOutput()->GetSpacing());
	  /*
	  typedef itk::RescaleIntensityImageFilter<InputImageType, OutputImageType> CastFilter;
      //CastFilter::Pointer caster = CastFilter::New();
	  caster->SetInput(fastMarching->GetOutput());
	  //WriterType::Pointer writer = WriterType::New();
	  writer->SetInput(caster->GetOutput());
	  writer->SetFileName("fastmarching.png");
	  writer->Update();
	  */
	  geodesicActiveContour->SetInput( fastMarching->GetOutput() );
	  geodesicActiveContour->SetFeatureImage( sigmoid->GetOutput() );
	  InputImageType *contour = geodesicActiveContour->GetOutput();
	  
	  typedef itk::BinaryThresholdImageFilter< InputImageType, OutputImageType > ThresholdingFilterType;
	  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
      thresholder->SetLowerThreshold( -1000.0 );
      thresholder->SetUpperThreshold( 0.0 );
      thresholder->SetOutsideValue( itk::NumericTraits< OutputPixelType >::min() );
      thresholder->SetInsideValue( itk::NumericTraits< OutputPixelType >::max() );
      thresholder->SetInput( contour );
	  try{
		  thresholder->Update();
	  }
	  catch(itk::ExceptionObject &e){
		  std::cerr<<e.what()<<'\n';
		  return EXIT_FAILURE;
	  }
	  
	  
	  WriterType::Pointer writer = WriterType::New();
	  writer->SetInput(thresholder->GetOutput());	  
	  writer->SetFileName("shenmingmou\\output"+std::to_string((2*reverse-1)*i+initialSliceIndex)+".png");
	  try{
          writer->Update();
	  }
	  catch(itk::ExceptionObject &e){
		  std::cerr<<e.what()<<'\n';
		  return EXIT_FAILURE;
	  }
	  
	  typedef itk::BinaryBallStructuringElement<OutputPixelType,2> BallStructureType;
	  BallStructureType kernel;
	  kernel.SetRadius(1);
	  kernel.CreateStructuringElement();

	  typedef itk::BinaryErodeImageFilter<OutputImageType,OutputImageType,BallStructureType> ErodeFilterType;
	  ErodeFilterType::Pointer erode = ErodeFilterType::New();
	  erode->SetInput(thresholder->GetOutput());
	  erode->SetKernel(kernel);

	  typedef itk::BinaryMorphologicalOpeningImageFilter<OutputImageType,OutputImageType,BallStructureType> OpenFilterType;
	  OpenFilterType::Pointer opener = OpenFilterType::New();
	  opener->SetInput(thresholder->GetOutput());
	  opener->SetKernel(kernel);

	  typedef itk::SubtractImageFilter<OutputImageType,OutputImageType> SubstractImageType;
	  SubstractImageType::Pointer substracter = SubstractImageType::New();
	  substracter->SetInput1(opener->GetOutput());
	  substracter->SetInput2(erode->GetOutput());
	  try{
		  substracter->Update();
	  }
	  catch(itk::ExceptionObject &e){
		  std::cerr<<e.what()<<'\n';
		  return EXIT_FAILURE;
	  }
	  
	  //updateSeeds(thresholder->GetOutput(),seeds);
	  /*itk::ImageRegionIterator<OutputImageType> it(thresholder->GetOutput(),thresholder->GetOutput()->GetRequestedRegion());
	  it.GoToBegin();
	  while(!it.IsAtEnd()){
		  std::cout<<it.Get()<<'\n';
		  it ++ ;
	  }*/
	  std::cout<<(2*reverse-1)*i+initialSliceIndex<<':'<<node1.GetIndex()[0]<<'\t'<<node1.GetIndex()[1]<<'\t'<<node1.GetValue()<<'\n';
	  //update seed positon to the middle
	  itk::ImageRegionIterator<OutputImageType> it(thresholder->GetOutput(),thresholder->GetOutput()->GetLargestPossibleRegion());
	  it.GoToBegin();
	  int left = 512, right = 0, up = 512, down = 0;
	  while(!it.IsAtEnd()){
		  if(it.Get() == 255){
			  OutputImageType::IndexType temp = it.GetIndex();
			  left = temp[0]<left?temp[0]:left;
			  right = temp[0]>right?temp[0]:right;
			  up = temp[1]<up?temp[1]:up;
			  down = temp[1]>down?temp[1]:down;
		  }
		  it ++;
	  }
	  if(reverse == 0 && (left>=right-2 || up>=down-2)) return EXIT_SUCCESS;
	  if(left>=right-2 || up>=down-2){
		  reverse = !reverse;
		  node1.SetIndex(seedPosition1);
		  node1.SetValue(-initialDistance);
		  seeds->SetElement(0,node1);
		  i=1;
		  continue;
	  }
	  seedPosX = (left + right)/2;
	  seedPosY = (up + down)/2;
	  OutputImageType::IndexType middlepoint;
	  middlepoint[0] = seedPosX;
	  middlepoint[1] = seedPosY;

	  OutputImageType::IndexType outpos = middlepoint;

	  if(thresholder->GetOutput()->GetPixel(outpos) == 255){
		  bool isIntersect = false;
		  int manhaton=1;
		  while(!isIntersect){
			  for(int i=0;i<manhaton;i++){
				  outpos[0] = middlepoint[0]+i;
				  outpos[1] = middlepoint[1]-manhaton+i;
				  if(thresholder->GetOutput()->GetPixel(outpos) == 0){
					  isIntersect = true;
					  break;
				  }
				  outpos[0] = middlepoint[0]+manhaton-i;
				  outpos[1] = middlepoint[1]+i;
				  if(thresholder->GetOutput()->GetPixel(outpos) == 0){
					  isIntersect = true;
					  break;
				  }
				  outpos[0] = middlepoint[0]-i;
				  outpos[1] = middlepoint[1]+manhaton-i;
				  if(thresholder->GetOutput()->GetPixel(outpos) == 0){
					  isIntersect = true;
					  break;
				  }
				  outpos[0] = middlepoint[0]-manhaton+i;
				  outpos[1] = middlepoint[1]-i;
				  if(thresholder->GetOutput()->GetPixel(outpos) == 0){
					  isIntersect = true;
					  break;
				  }
				  manhaton ++;
			  }
		  }
		  seedValue = -pow((outpos[0]-middlepoint[0])*(outpos[0]-middlepoint[0])+(outpos[1]-middlepoint[1])*(outpos[1]-middlepoint[1]),0.5)+1;
	  }
	  else if(thresholder->GetOutput()->GetPixel(outpos) == 0){
		  bool isIntersect = false;
		  int manhaton=1;
		  while(!isIntersect){
			  for(int i=0;i<manhaton;i++){
				  outpos[0] = middlepoint[0]+i;
				  outpos[1] = middlepoint[1]-manhaton+i;
				  if(thresholder->GetOutput()->GetPixel(outpos) == 255){
					  isIntersect = true;
					  break;
				  }
				  outpos[0] = middlepoint[0]+manhaton-i;
				  outpos[1] = middlepoint[1]+i;
				  if(thresholder->GetOutput()->GetPixel(outpos) == 255){
					  isIntersect = true;
					  break;
				  }
				  outpos[0] = middlepoint[0]-i;
				  outpos[1] = middlepoint[1]+manhaton-i;
				  if(thresholder->GetOutput()->GetPixel(outpos) == 255){
					  isIntersect = true;
					  break;
				  }
				  outpos[0] = middlepoint[0]-manhaton+i;
				  outpos[1] = middlepoint[1]-i;
				  if(thresholder->GetOutput()->GetPixel(outpos) == 255){
					  isIntersect = true;
					  break;
				  }
				  manhaton ++;
			  }
		  }
		  seedValue = pow((outpos[0]-middlepoint[0])*(outpos[0]-middlepoint[0])+(outpos[1]-middlepoint[1])*(outpos[1]-middlepoint[1]),0.5)-1;
	  }
	  node1.SetIndex(middlepoint);
	  node1.SetValue(seedValue);
	  seeds->SetElement(0,node1);
  }
  


  /*
  typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType > CastFilterType;

  CastFilterType::Pointer caster1 = CastFilterType::New();
  CastFilterType::Pointer caster2 = CastFilterType::New();
  CastFilterType::Pointer caster3 = CastFilterType::New();
  CastFilterType::Pointer caster4 = CastFilterType::New();

  WriterType::Pointer writer1 = WriterType::New();
  WriterType::Pointer writer2 = WriterType::New();
  WriterType::Pointer writer3 = WriterType::New();
  WriterType::Pointer writer4 = WriterType::New();

  caster1->SetInput( smoothing->GetOutput() );
  writer1->SetInput( caster1->GetOutput() );
  writer1->SetFileName("GeodesicActiveContourImageFilterOutput1_test.png");
  caster1->SetOutputMinimum( itk::NumericTraits< OutputPixelType >::min() );
  caster1->SetOutputMaximum( itk::NumericTraits< OutputPixelType >::max() );
  writer1->Update();

  caster2->SetInput( gradientMagnitude->GetOutput() );
  writer2->SetInput( caster2->GetOutput() );
  writer2->SetFileName("GeodesicActiveContourImageFilterOutput2_test.png");
  caster2->SetOutputMinimum( itk::NumericTraits< OutputPixelType >::min() );
  caster2->SetOutputMaximum( itk::NumericTraits< OutputPixelType >::max() );
  writer2->Update();

  caster3->SetInput( sigmoid->GetOutput() );
  writer3->SetInput( caster3->GetOutput() );
  writer3->SetFileName("GeodesicActiveContourImageFilterOutput3_test.png");
  caster3->SetOutputMinimum( itk::NumericTraits< OutputPixelType >::min() );
  caster3->SetOutputMaximum( itk::NumericTraits< OutputPixelType >::max() );
  writer3->Update();

  caster4->SetInput( fastMarching->GetOutput() );
  writer4->SetInput( caster4->GetOutput() );
  writer4->SetFileName("GeodesicActiveContourImageFilterOutput4_test.png");
  caster4->SetOutputMinimum( itk::NumericTraits< OutputPixelType >::min() );
  caster4->SetOutputMaximum( itk::NumericTraits< OutputPixelType >::max() );

  fastMarching->SetOutputSize(
	  reader->GetOutput()->GetLargestPossibleRegion().GetSize() );

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName );
  writer->SetInput( thresholder->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << std::endl;
  std::cout << "Max. no. iterations: " << geodesicActiveContour->GetNumberOfIterations() << std::endl;
  std::cout << "Max. RMS error: " << geodesicActiveContour->GetMaximumRMSError() << std::endl;
  std::cout << std::endl;
  std::cout << "No. elpased iterations: " << geodesicActiveContour->GetElapsedIterations() << std::endl;
  std::cout << "RMS change: " << geodesicActiveContour->GetRMSChange() << std::endl;

  try
    {
    writer4->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
*/
  return EXIT_SUCCESS;
}
