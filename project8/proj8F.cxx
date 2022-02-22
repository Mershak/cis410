
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkContourFilter.h>
#include <vtkCutter.h>
#include <vtkPlane.h>

int main(int argc, char *argv[])
{
   vtkDataSetReader *reader = vtkDataSetReader::New();
   reader->SetFileName("noise.vtk");
   reader->Update();

   vtkPlane *plane = vtkPlane::New();
   plane->SetNormal(0,0,1);

   vtkCutter *cutter = vtkCutter::New();
   cutter->SetInputConnection(reader->GetOutputPort());
   cutter->SetCutFunction(plane);
   cutter->Update();


   vtkContourFilter *cf = vtkContourFilter::New();
   cf->SetInputConnection(reader->GetOutputPort());
   cf->SetNumberOfContours(2);
   cf->SetValue(0, 2.4);
   cf->SetValue(1, 4);

   vtkLookupTable *lut = vtkLookupTable::New();
   for (int i = 0; i < 256; i++){
      lut->SetTableValue(i, (255*i / 255), 0, (255*(255-i)/255));
   }
   lut->Build();

   vtkDataSetMapper *mapper = vtkDataSetMapper::New();
   mapper->SetInputConnection(cutter->GetOutputPort());
   mapper->SetResolveCoincidentTopologyToPolygonOffset();
   mapper->SetLookupTable(lut);
   mapper->SetScalarRange(1, 6);


   vtkDataSetMapper *mapper2 = vtkDataSetMapper::New();
   mapper2->SetInputConnection(cf->GetOutputPort());
   mapper2->SetLookupTable(lut);
   mapper2->SetScalarRange(1, 6);

   vtkActor *actor = vtkActor::New();
   actor->SetMapper(mapper);

   vtkActor *actor2 = vtkActor::New();
   actor2->SetMapper(mapper2);

   vtkRenderer *ren = vtkRenderer::New();
   ren->AddActor(actor);
   ren->SetViewport(0,0,0.5,1);

   vtkRenderer *ren2 = vtkRenderer::New();
   ren2->AddActor(actor2);
   ren2->SetViewport(0.5,0,1,1);


   vtkRenderWindow *renwin = vtkRenderWindow::New();
   renwin->AddRenderer(ren);
   renwin->AddRenderer(ren2);
   renwin->SetSize(768,768);

   vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
   iren->SetRenderWindow(renwin);
   renwin->Render();
   iren->Start();
}


