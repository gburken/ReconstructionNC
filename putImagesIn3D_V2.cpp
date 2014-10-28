#include <itkNumericSeriesFileNames.h>
#include <itkImage.h>
#include <itkImageSeriesReader.h>
#include <itkChangeInformationImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.hxx>
#include <itkMersenneTwisterRandomVariateGenerator.h>

#include <itkEuler3DTransform.h>
#include <itkResampleImageFilter.h>
#include <itkImageFileWriter.h>

// read plus file
#include <iostream>
#include <fstream>
#include <string>

static double pi = 3.14159265359;

int main(int argc, char** argv)
{
 
  // fundamental typedefs
  typedef unsigned char  PixelType;
  typedef itk::Image<PixelType, 2>  Char2dImageType;  
  typedef itk::Image<PixelType, 3>  Char3dImageType;
    
  // test reading Plus images
  std::string plusFileName = "C:/Users/ahaak/PlusApp-2.1.2.3381-Win32/data/TrackedImageSequence_20141014_163846.mha";
  typedef itk::ImageFileReader<Char3dImageType> PlusReaderType;
  PlusReaderType::Pointer plusReader = PlusReaderType::New();
  plusReader->SetFileName( plusFileName );
  try
  {
    plusReader->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << plusReader->GetUseStreaming();
  Char3dImageType::Pointer test = plusReader->GetOutput();

  // more testing on those Plus images
  std::ifstream plusFile( plusFileName );
  std::string line;
  if ( plusFile.is_open() )
  {
    for ( int itr = 0; itr < 15; ++itr )
    {
      std::getline( plusFile, line );
      std::cout << line << std::endl;
    }
    while(1)
    {
      std::getline( plusFile, line );
      std::size_t pos = 0;
      std::size_t len = 21;
      int idx1 = line.compare( pos, len, "ProbeToTrackerTransform" );
      std::cout << line << std::endl;
    }
    
  }

  // read image
  std::string baseFileName = "D:/Data/ExperimentalData/EPExperiments/EP1_24-Feb-2014/consecutiveInputDataNC/IM_0007/600Images/2DinputImage_";
  typedef itk::NumericSeriesFileNames FileNameType;
  FileNameType::Pointer fileNames = FileNameType::New();
  fileNames->SetSeriesFormat( baseFileName + "%04d.png" );
  fileNames->SetStartIndex( 1 );
  fileNames->SetEndIndex( 10 );
    
  
  typedef itk::ImageSeriesReader<Char3dImageType> CharReaderType;
  CharReaderType::Pointer reader = CharReaderType::New();
  reader->SetFileNames( fileNames->GetFileNames() );
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
  Char3dImageType::Pointer inImage = reader->GetOutput();

  // create transform
  typedef itk::Euler3DTransform<double> TransformTyp;
  TransformTyp::Pointer transform = TransformTyp::New();
  TransformTyp::ParametersType transformParameters = transform->GetParameters();
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize();
  
  
  transformParameters[0] = generator->GetUniformVariate(0.0, 180.0)/180.0*pi; 
  transformParameters[1] = generator->GetUniformVariate(0.0, 180.0)/180.0*pi; 
  transformParameters[2] = generator->GetUniformVariate(0.0, 180.0)/180.0*pi;
  transform->SetParameters( transformParameters );
 
  // create the output size, spacing etc.
  Char3dImageType::IndexType outIndex;
  outIndex.Fill( 0 );
  Char3dImageType::SizeType outSize;
  outSize.Fill( 600 );
  Char3dImageType::RegionType outRegion;
  outRegion.SetIndex( outIndex );
  outRegion.SetSize( outSize );
  Char3dImageType::SpacingType outSpacing;
  outSpacing.Fill( 1.0 );
  Char3dImageType::PointType outOrigin;
  outOrigin.Fill( -300.0 );

  Char3dImageType::Pointer cmap = Char3dImageType::New();
  cmap->SetRegions( outRegion );
  cmap->SetSpacing( outSpacing );
  cmap->SetOrigin( outOrigin );
  cmap->Allocate();
  cmap->FillBuffer( 0 );

  Char3dImageType::Pointer imap = Char3dImageType::New();
  imap->SetRegions( outRegion );
  imap->SetSpacing( outSpacing );
  imap->SetOrigin( outOrigin );
  imap->Allocate();
  imap->FillBuffer( 0 );

  // region iterator, points, and regions
  typedef itk::ImageRegionConstIteratorWithIndex<Char3dImageType> inIteratorType;
  typedef itk::ImageRegionIteratorWithIndex<Char3dImageType> outIteratorType;

  Char3dImageType::RegionType inSliceRegion = reader->GetOutput()->GetLargestPossibleRegion();
  Char3dImageType::SizeType inSliceSize = inSliceRegion.GetSize();
  inSliceSize[2] = 1;
  inSliceRegion.SetSize( inSliceSize );
  Char3dImageType::IndexType inSliceIndex = inSliceRegion.GetIndex();
  
  Char3dImageType::PointType inPt;
  Char3dImageType::PointType outPt;

  // iterate through the input image stag
  for( int itr = 0; itr < inSize[2]; itr++)
  {
    // create random transform
    transformParameters[0] = generator->GetUniformVariate(0.0, 180.0)/180.0*pi; 
    transformParameters[1] = generator->GetUniformVariate(0.0, 180.0)/180.0*pi; 
    transformParameters[2] = generator->GetUniformVariate(0.0, 180.0)/180.0*pi;
    transform->SetParameters( transformParameters );

    // iterate through all pixel in input slice 
    inSliceIndex[2] = itr;
    inSliceRegion.SetIndex( inSliceIndex );
    inIteratorType inItr( reader->GetOutput(), inSliceRegion );
    inItr.GoToBegin();

    while( !inItr.IsAtEnd() )
    {
      inImage->TransformIndexToPhysicalPoint( inItr.GetIndex(), inPt );
      outPt = transform->TransformPoint( inPt );
      bool isInside = cmap->TransformPhysicalPointToIndex(outPt, outIndex);
      if ( isInside )
      {
        PixelType cValue = cmap->GetPixel( outIndex ); 
        if ( cValue > 0 )
        {
          PixelType iValue = imap->GetPixel( outIndex );
          PixelType inValue = inItr.Get();
          imap->SetPixel( outIndex, (PixelType)(iValue + inValue)/2.0 );
          cmap->SetPixel( outIndex, 1 );
        } 
        else
        {
          imap->SetPixel( outIndex, inItr.Get() );
          cmap->SetPixel( outIndex, 1 );
        }
      }
      ++inItr;
    }
  }
  

  // write 3Dimage to file
  typedef itk::ImageFileWriter<Char3dImageType> WriterType;

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( imap );
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