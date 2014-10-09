#include "vtkSmartPointer.h"
#include "vtkPNGReader.h"
#include "vtkLookupTable.h" // not needed?
#include <vtkImageMapper3D.h>
#include <vtkImageMapToColors.h>
#include "vtkImageActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleImage.h"

#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkImageViewer.h"





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

  vtkSmartPointer<vtkMatrix4x4> resliceAxes = vtkSmartPointer<vtkMatrix4x4>::New();
  resliceAxes->DeepCopy(XZPlaneElements);
  resliceAxes->SetElement(0, 3, 0.0 );resliceAxes->SetElement(1, 3, 0.0 );resliceAxes->SetElement(2, 3, 0.0 );
  
  
  vtkSmartPointer<vtkTransform> resliceAxesTransform = vtkSmartPointer<vtkTransform>::New();
  resliceAxesTransform->SetMatrix( resliceAxes );
  resliceAxesTransform->Update();
  

  resliceAxes->DeepCopy(YZPlaneElements);
  resliceAxes->SetElement(0, 3, 0.0 );resliceAxes->SetElement(1, 3, 0.0 );resliceAxes->SetElement(2, 3, 0.0 );

  vtkSmartPointer<vtkTransform> resliceAxesTransform2 = vtkSmartPointer<vtkTransform>::New();
  resliceAxesTransform2->SetMatrix( resliceAxes );
  resliceAxesTransform2->Update();

  // read image
  std::string fileName = "D:/Data/ExperimentalData/EPExperiments/EP1_24-Feb-2014/consecutiveInputDataNC/IM_0007/600Images/2DinputImage_0584.png";
  vtkSmartPointer<vtkPNGReader> reader = vtkPNGReader::New();
  reader->SetFileName( fileName.c_str() );
  reader->Update();
  
  // Create a greyscale lookup table
  vtkSmartPointer<vtkLookupTable> table = vtkSmartPointer<vtkLookupTable>::New();
  table->SetRange(0, 255); // image intensity range
  table->SetValueRange(0.0, 1.0); // from black to white
  table->SetSaturationRange(0.0, 0.0); // no color saturation
  table->SetRampToLinear();
  table->Build();

  // mapper
  /*vtkSmartPointer<vtkImageMapper3D> mapper = vtkSmartPointer<vtkImageMapper3D>::New();
  mapper->SetInputConnection( reader->GetOutputPort() ); 
 */
  // mapper 2
  vtkSmartPointer<vtkImageMapToColors> mapper2 = vtkSmartPointer<vtkImageMapToColors>::New();
  mapper2->SetLookupTable( table );
  mapper2->SetInputConnection( reader->GetOutputPort() );
  // actor
  vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();
  actor->SetUserTransform( resliceAxesTransform );
  actor->GetMapper()->SetInputConnection( mapper2->GetOutputPort() );
  
  vtkSmartPointer<vtkImageActor> actor2 = vtkSmartPointer<vtkImageActor>::New();
  actor2->SetUserTransform( resliceAxesTransform2 );
  actor2->GetMapper()->SetInputConnection( mapper2->GetOutputPort() );

  // renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(actor);
  renderer->AddActor(actor2);
  renderer->SetViewport(0.5,0.0,1,0.5);

  //render window
  vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
  window->AddRenderer(renderer);

  // Set up the interaction
  vtkSmartPointer<vtkInteractorStyleImage> imageStyle = vtkSmartPointer<vtkInteractorStyleImage>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> interactor =  vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetInteractorStyle( imageStyle );
  window->SetInteractor( interactor );
  window->Render();
  window->Start();

  interactor->Start();

  std::cin.ignore();
  return 0;
}