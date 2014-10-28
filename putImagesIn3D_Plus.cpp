#include <itkNumericSeriesFileNames.h>
#include <itkImage.h>
#include <itkImageSeriesReader.h>
#include <itkChangeInformationImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

#include "vtkTrackedFrameList.h"
#include "TrackedFrame.h"
#include <itkEuler3DTransform.h>
#include <itkAffineTransform.h>
#include "vtkMetaImageSequenceIO.h"
#include <itkImageFileWriter.h>
//#include <itkGPUDiscreteGaussianImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkDivideImageFilter.h>

// read plus file
#include <iostream>
#include <fstream>
#include <string>

static double pi = 3.14159265359;
itk::Euler3DTransform<double>::MatrixType vtkMatrix4x4ToItkRotMatrix( vtkSmartPointer<vtkMatrix4x4> transformMatrix );
itk::Euler3DTransform<double>::OffsetType vtkMatrix4x4ToItkOffset( vtkSmartPointer<vtkMatrix4x4> transformMatrix );
std::vector<double> computeExtent( vtkSmartPointer<vtkTrackedFrameList> frames );

int main(int argc, char** argv)
{
 
  // fundamental typedefs
  typedef unsigned char  PixelType;
  typedef itk::Image<PixelType, 2>  Char2dImageType;  
  typedef itk::Image<PixelType, 3>  Char3dImageType;
    
  // reading Plus images
  std::string plusFileName = "C:/Users/ahaak/PlusApp-2.1.2.3381-Win32/data/initalShortSweep3.mha";
  vtkSmartPointer<vtkTrackedFrameList> frames = vtkSmartPointer<vtkTrackedFrameList>::New();
  frames->ReadFromSequenceMetafile( plusFileName.c_str() );
  int numFrames = frames->GetNumberOfTrackedFrames();
  TrackedFrame* frame = frames->GetTrackedFrame( 0 );
  std::vector<double> extendVolume = computeExtent( frames );
   
  // create the output size, spacing etc.
  typedef itk::Image<float, 3> Float3dImageType;
  Float3dImageType::IndexType outIndex;
  outIndex.Fill( 0 );
  Float3dImageType::SizeType outSize;
  for (int i = 0; i < 3; ++i)
  {
    outSize[i] = ceil(extendVolume[2*i+1] - extendVolume[2*i]);
  }
  
  Float3dImageType::RegionType outRegion;
  outRegion.SetIndex( outIndex );
  outRegion.SetSize( outSize );
  Float3dImageType::SpacingType outSpacing;
  outSpacing.Fill( 1.0 );
  Float3dImageType::PointType outOrigin;
  outOrigin[0] = 100;//extendVolume[0]*1;
  outOrigin[1] = extendVolume[2]*1;
  outOrigin[2] = extendVolume[4]*1;
  
  /*outOrigin[0] = -300;
  outOrigin[1] = -180;
  outOrigin[2] = -450;*/

  /*curved2.mha
  outOrigin[0] = -300;
  outOrigin[1] = -80;
  outOrigin[2] = -450;*/
  /*curved1.mha 
  outOrigin[0] = -100;
  outOrigin[1] = -80;
  outOrigin[2] = -650;*/
  /*outOrigin[0] = -262.712;
  outOrigin[1] = -80.4578;
  outOrigin[2] = -525.568;*/
    
  Float3dImageType::Pointer cmap = Float3dImageType::New();
  cmap->SetRegions( outRegion );
  cmap->SetSpacing( outSpacing );
  cmap->SetOrigin( outOrigin );
  cmap->Allocate();
  cmap->FillBuffer( 0 );

  Float3dImageType::Pointer imap = Float3dImageType::New();
  imap->SetRegions( outRegion );
  imap->SetSpacing( outSpacing );
  imap->SetOrigin( outOrigin );
  imap->Allocate();
  imap->FillBuffer( 0 );

  // create transform
  std::vector<PlusTransformName> transformNames;
  vtkSmartPointer<vtkMatrix4x4> transformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  typedef itk::AffineTransform<double> TransformTyp;
  //typedef itk::Euler3DTransform<double> TransformTyp;
  TransformTyp::Pointer transform = TransformTyp::New();
  transform->SetIdentity();

  // region iterator, points, and regions
  typedef itk::ImageRegionConstIteratorWithIndex<Char3dImageType> inConstIteratorType;
  typedef itk::ImageRegionIteratorWithIndex<Char3dImageType> inIteratorType;
  typedef itk::ImageRegionConstIteratorWithIndex<Char2dImageType> inIterator2dType;
  typedef itk::ImageRegionIteratorWithIndex<Float3dImageType> outIteratorType;

  Char2dImageType::RegionType region2dFrame = frames->GetTrackedFrame( 0 )->GetImageData()->GetDisplayableImage()->GetLargestPossibleRegion();
  Char3dImageType::RegionType inSliceRegion; 
  Char3dImageType::SizeType inSliceSize = inSliceRegion.GetSize();
  inSliceSize[0] = region2dFrame.GetSize()[0];
  inSliceSize[1] = region2dFrame.GetSize()[1]; 
  inSliceSize[2] = 1;
  inSliceRegion.SetSize( inSliceSize );
  Char3dImageType::IndexType inSliceIndex = inSliceRegion.GetIndex();
    
  Char3dImageType::PointType inPt;
  Char2dImageType::PointType inFramePt;
  Float3dImageType::PointType outPt;

  // create mask image
  Char3dImageType::Pointer maskFrame = Char3dImageType::New();
  maskFrame->SetRegions( inSliceRegion );
  Char3dImageType::PointType maskOrigin;
  maskOrigin.Fill( 0.0 );
  maskOrigin[0] = -(double)inSliceSize[0]/2.0;
  maskFrame->SetOrigin( maskOrigin );
  maskFrame->Allocate();
  maskFrame->FillBuffer( 0 );

  double openingAngleTan = atan( 75.0 / 180 * pi );
  double radius = 528;
  double eps = std::numeric_limits<double>::epsilon();

  inIteratorType maskItr( maskFrame, inSliceRegion );
  maskItr.GoToBegin();
  Char3dImageType::PointType maskPt;

  while( !maskItr.IsAtEnd() )
  {
    maskFrame->TransformIndexToPhysicalPoint( maskItr.GetIndex(), maskPt );
    double currentTan = abs( maskPt[1] / (maskPt[0] + eps) );
    double currentRadius = sqrt( maskPt[0]*maskPt[0] + maskPt[1]*maskPt[1] );
    if ( currentTan >= openingAngleTan && currentRadius <= radius )
    {
      maskItr.Set( 1 );
    }
    ++maskItr;
  }

  // write mask for testing
  typedef itk::ImageFileWriter<Char3dImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( maskFrame );
  writer->SetFileName( "maskFrame.mhd" );
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


  // iterate through the input image stag
  for( int itr = 0; itr < frames->GetNumberOfTrackedFrames(); itr++)
  {
    std::cout << "Frame: " << itr << std::endl;
    // get frame
    TrackedFrame* frame = frames->GetTrackedFrame( itr );
    Char2dImageType::Pointer itkFrame2D = frame->GetImageData()->GetDisplayableImage();

    // get transform
    frame->GetCustomFrameTransformNameList( transformNames );
    frame->GetCustomFrameTransform( transformNames[0], transformMatrix);
    double time = frame->GetTimestamp();
    transform->SetIdentity();
    transform->SetMatrix( vtkMatrix4x4ToItkRotMatrix(transformMatrix) );
    transform->SetOffset( vtkMatrix4x4ToItkOffset(transformMatrix) );

    // iterate through all pixel in input slice 
    inIterator2dType inItr( itkFrame2D, itkFrame2D->GetLargestPossibleRegion() );
    inConstIteratorType maskItr( maskFrame, maskFrame->GetLargestPossibleRegion() );
    inItr.GoToBegin();
    maskItr.GoToBegin();

    while( !inItr.IsAtEnd() )
    {
      if ( maskItr.Get() > 0 )
      {
	      maskFrame->TransformIndexToPhysicalPoint( maskItr.GetIndex(), inPt );
        outPt = transform->TransformPoint( inPt );
	      bool isInside = cmap->TransformPhysicalPointToIndex(outPt, outIndex);
	      if ( isInside )
	      {
	        Float3dImageType::PixelType cValue = cmap->GetPixel( outIndex ); 
	        if ( cValue > 0.0 )
	        {
	          Float3dImageType::PixelType iValue = imap->GetPixel( outIndex );
	          Float3dImageType::PixelType inValue = inItr.Get();
	          imap->SetPixel( outIndex, (Float3dImageType::PixelType)(iValue + inValue)/2.0 );
	          cmap->SetPixel( outIndex, 1.0 );
	        } 
	        else
	        {
	          imap->SetPixel( outIndex, (Float3dImageType::PixelType) inItr.Get() );
	          cmap->SetPixel( outIndex, 1.0 );
	        }
	      }
      }
      ++maskItr;
      ++inItr;
    }
  }
  
  // do the blurring
  typedef itk::DiscreteGaussianImageFilter<Float3dImageType, Float3dImageType> GaussianFilterType;
  GaussianFilterType::Pointer iBlurrer = GaussianFilterType::New();
  GaussianFilterType::Pointer cBlurrer = GaussianFilterType::New();
  double variance = 20;
  iBlurrer->SetVariance( variance );
  cBlurrer->SetVariance( variance );
  iBlurrer->SetInput( imap );//iCaster->GetOutput()
  cBlurrer->SetInput( cmap );//cCaster->GetOutput()
  iBlurrer->SetNumberOfThreads( 16 );
  cBlurrer->SetNumberOfThreads( 16 );

  try
  {
    iBlurrer->Update();
    cBlurrer->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  // save data
  typedef itk::ImageFileWriter<Float3dImageType> FloatWriterType;
  FloatWriterType::Pointer floatWriter = FloatWriterType::New();
  floatWriter->SetInput( iBlurrer->GetOutput() );
  floatWriter->SetFileName( "iBlurrer.mhd" );
  try
  {
    floatWriter->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  // save data
  floatWriter->SetInput( cBlurrer->GetOutput() );
  floatWriter->SetFileName( "cBlurrer.mhd" );
  try
  {
    floatWriter->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  
  // divide the two maps
  Float3dImageType::Pointer ncMap = Float3dImageType::New();
  ncMap->SetRegions( imap->GetLargestPossibleRegion() );
  ncMap->SetSpacing( imap->GetSpacing() );
  ncMap->SetOrigin( imap->GetOrigin() );
  ncMap->Allocate();
  ncMap->FillBuffer( 0.0 );

  outIteratorType iMapItr( imap, imap->GetLargestPossibleRegion() );
  outIteratorType cMapItr( cmap, cmap->GetLargestPossibleRegion() );
  outIteratorType ncMapItr( ncMap, ncMap->GetLargestPossibleRegion() );

  iMapItr.GoToBegin();
  cMapItr.GoToBegin();
  ncMapItr.GoToBegin();

  while( !iMapItr.IsAtEnd() )
  {
    if ( cMapItr.Get() > std::numeric_limits<float>::epsilon() )
    {
      ncMapItr.Set( iMapItr.Get() / cMapItr.Get() );
    }
    ++iMapItr;
    ++cMapItr;
    ++ncMapItr;
  }

    
  // save data
  floatWriter->SetInput( ncMap );
  floatWriter->SetFileName( "ncImage.mhd" );
  try
  {
    floatWriter->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Done!" << std::endl;
  std::cin.ignore();
  return 0;
}

itk::Euler3DTransform<double>::MatrixType vtkMatrix4x4ToItkRotMatrix( vtkSmartPointer<vtkMatrix4x4> transformMatrix)
{
  itk::Euler3DTransform<double>::MatrixType rotationMatrix;
  for ( int row = 0; row < 3; row++ )
  {
    for ( int col = 0; col < 3; col++ )
    {
      rotationMatrix[row][col] = transformMatrix->GetElement( row, col ); 
    }
  }
  return rotationMatrix;
}

itk::Euler3DTransform<double>::OffsetType vtkMatrix4x4ToItkOffset( vtkSmartPointer<vtkMatrix4x4> transformMatrix)
{
  itk::Euler3DTransform<double>::OffsetType translation;
  translation.Fill( 0.0 );
  for ( int row = 0; row < 3; row++ )
  {
    translation[row] =  transformMatrix->GetElement( row, 3 ); 
  }
  return translation;
}

std::vector<double> computeExtent( vtkSmartPointer<vtkTrackedFrameList> frames )
{
  std::vector<double> returnValues;
  std::vector<double> minExtend(3, 1000), maxExtend(3, -1000);
  
  // define extend of 2D plane (Values below are for Ultransonix with 15 cm depth)
  typedef itk::Image<unsigned char, 3>::PointType PointType;
  PointType pt1, pt2, pt3, pt4, outPt;
  pt1[0] = 0; pt1[1] = 415; pt1[2] = 0;
  pt2[0] = 324; pt2[1] = 528; pt2[2] = 0;
  pt3[0] = 646; pt3[1] = 415; pt3[2] = 0;
  pt4[0] = 324; pt4[1] = 0; pt4[2] = 0;
  std::vector<PointType> pts;
  pts.push_back(pt1); pts.push_back(pt2); pts.push_back(pt3); pts.push_back(pt4); 
  
  // transform
  std::vector<PlusTransformName> transformNames;
  vtkSmartPointer<vtkMatrix4x4> transformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  typedef itk::AffineTransform<double> TransformTyp;
  TransformTyp::Pointer transform = TransformTyp::New();

  for (int itr = 0; itr < frames->GetNumberOfTrackedFrames(); ++itr )
  {
    // get current transform
    frames->GetTrackedFrame( itr )->GetCustomFrameTransformNameList( transformNames );
    frames->GetTrackedFrame( itr )->GetCustomFrameTransform( transformNames[0], transformMatrix);

    transform->SetMatrix( vtkMatrix4x4ToItkRotMatrix(transformMatrix) );
    transform->SetOffset( vtkMatrix4x4ToItkOffset(transformMatrix) );

    // loop through the 4 land mark points and check limits
    for ( int ptItr = 0; ptItr < 4; ++ptItr )
    {
      outPt = transform->TransformPoint( pts[ptItr] );
      for ( int compItr = 0; compItr < 3; ++compItr )
      {
        minExtend[compItr] = std::min( minExtend[compItr], outPt[compItr]);
        maxExtend[compItr] = std::max( maxExtend[compItr], outPt[compItr]);
      }
    }
  }
  returnValues.push_back( minExtend[0] ); returnValues.push_back( maxExtend[0] );
  returnValues.push_back( minExtend[1] ); returnValues.push_back( maxExtend[1] );
  returnValues.push_back( minExtend[2] ); returnValues.push_back( maxExtend[2] );
  return returnValues;
}