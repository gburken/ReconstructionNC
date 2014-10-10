#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkEuler3DTransform.h>
#include <itkResampleImageFilter.h>
#include <itkImageFileWriter.h>

static double pi = 3.14159265359;

int main(int argc, char** argv)
{
  // some additional stuff
  static double XYPlaneElements[16] = {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1 };


  static double XZPlaneElements[16] = {
    1, 0, 0, 0,
    0, 0, 1, 0,
    0,-1, 0, 0,
    0, 0, 0, 1 };

  static double YZPlaneElements[16] = {
    0, 0,-1, 0,
    1, 0, 0, 0,
    0,-1, 0, 0,
    0, 0, 0, 1 };

  // fundamental typedefs
  typedef unsigned char  PixelType;
  typedef itk::Image<PixelType, 2>  Char2dImageType;
  typedef itk::Image<PixelType, 3>  Char3dImageType;
      

  // read image
  std::string fileName = "D:/Data/ExperimentalData/EPExperiments/EP1_24-Feb-2014/consecutiveInputDataNC/IM_0007/600Images/2DinputImage_0584.png";
  typedef itk::ImageFileReader<Char3dImageType> CharReaderType;
  CharReaderType::Pointer reader = CharReaderType::New();
  reader->SetFileName( fileName );
  try
  {
    reader->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  Char3dImageType::PointType inOrigin = reader->GetOutput()->GetOrigin();
  Char3dImageType::SpacingType inSpacing = reader->GetOutput()->GetSpacing();
  Char3dImageType::SizeType inSize = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  inOrigin[0] = (double)(inSize[0]/-2.0)*inSpacing[0];
  reader->GetOutput()->SetOrigin( inOrigin );
  reader->GetOutput()->Print( std::cout );
    
  // create transform
  typedef itk::Euler3DTransform<double> TransformTyp;
  TransformTyp::Pointer transform = TransformTyp::New();
  TransformTyp::ParametersType transformParameters = transform->GetParameters();
  transformParameters[0] = 90.0/180.0*pi; transformParameters[1] = 45.0/180.0*pi; transformParameters[2] = 45.0/180.0*pi;
  std::cout << transformParameters << std::endl;
  
  // create the output size, spacing etc
  Char3dImageType::IndexType outIndex;
  outIndex.Fill( 0 );
  Char3dImageType::SizeType outSize;
  outSize.Fill( 600 );
  Char3dImageType::SpacingType outSpacing;
  outSpacing.Fill( 1.0 );
  Char3dImageType::PointType outOrigin;
  outOrigin.Fill( -300.0 );

  // create resample filter
  typedef itk::ResampleImageFilter<Char3dImageType, Char3dImageType, double> ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetInput( reader->GetOutput() );
  resampler->SetTransform( transform );
  resampler->SetOutputStartIndex( outIndex );
  resampler->SetSize( outSize );
  resampler->SetOutputSpacing( outSpacing );
  resampler->SetOutputOrigin( outOrigin );
  try
  {
    resampler->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  // write 3Dimage to file
  typedef itk::ImageFileWriter<Char3dImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( resampler->GetOutput() );
  writer->SetFileName( "outImage.mhd" );
  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  std::cin.ignore();
  return 0;
}