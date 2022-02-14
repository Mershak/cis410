/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include "LUT.h"

// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    // idx[0] = pointId%dims[0];
    // idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}


class TriangleList
{
   public:
                   TriangleList() { maxTriangles = 1000000; triangleIdx = 0; pts = new float[9*maxTriangles]; };
     virtual      ~TriangleList() { delete [] pts; };

     void          AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxTriangles;
     int           triangleIdx;
};

void
TriangleList::AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3)
{
    pts[9*triangleIdx+0] = X1;
    pts[9*triangleIdx+1] = Y1;
    pts[9*triangleIdx+2] = Z1;
    pts[9*triangleIdx+3] = X2;
    pts[9*triangleIdx+4] = Y2;
    pts[9*triangleIdx+5] = Z2;
    pts[9*triangleIdx+6] = X3;
    pts[9*triangleIdx+7] = Y3;
    pts[9*triangleIdx+8] = Z3;
    triangleIdx++;
}

vtkPolyData *
TriangleList::MakePolyData(void)
{
    int ntriangles = triangleIdx;
    int numPoints = 3*(ntriangles);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *tris = vtkCellArray::New();
    tris->EstimateSize(numPoints,4);
    for (int i = 0 ; i < ntriangles ; i++)
    {
        double pt[3];
        pt[0] = pts[9*i];
        pt[1] = pts[9*i+1];
        pt[2] = pts[9*i+2];
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[9*i+3];
        pt[1] = pts[9*i+4];
        pt[2] = pts[9*i+5];
        vtk_pts->SetPoint(ptIdx+1, pt);
        pt[0] = pts[9*i+6];
        pt[1] = pts[9*i+7];
        pt[2] = pts[9*i+8];
        vtk_pts->SetPoint(ptIdx+2, pt);
        vtkIdType ids[3] = { ptIdx, ptIdx+1, ptIdx+2 };
        tris->InsertNextCell(3, ids);
        ptIdx += 3;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetPolys(tris);
    tris->Delete();
    vtk_pts->Delete();

    return pd;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj7.vtk");
    rdr->Update();
    if (rdr->GetOutput() == NULL || rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Could not find input file." << endl;
        exit(EXIT_FAILURE);
    }

    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    int dims[3];
    rgrid->GetDimensions(dims);
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    // These were useful to me
    int edgeToVertex[12][2] =
        {
            {  0,  1 },
            {  2,  1 },
            {  2,  3 },
            {  0,  3 },
            {  4,  5 },
            {  5,  6 },
            {  6,  7 },
            {  4,  7 },
            {  0,  4 },
            {  1,  5 },
            {  3,  7 },
            {  2,  6 }
         };
    // This follows the convention in Lecture 11 slides (and project 6)
    // X is left-to-right, Y is up-and-down, Z is front-and-back.
    int offsetsI[8] = { 0, 1, 1, 0, 0, 1, 1, 0 };
    int offsetsJ[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };
    int offsetsK[8] = { 0, 0, 1, 1, 0, 0, 1, 1 };

    TriangleList tl;
    int ncells = rgrid->GetNumberOfCells();
    cerr << "Number of cells to isosurface is " << ncells << endl;
    float isoval = 3.2;
    for (int i = 0 ; i < ncells ; i++)
    {
      int idx[3];
      GetLogicalCellIndex(idx,i,dims);
      int logicalv0[3] = {idx[0], idx[1], idx[2]};
      int logicalv1[3] = {idx[0]+1, idx[1], idx[2]};
      int logicalv2[3] = {idx[0]+1, idx[1], idx[2]+1};
      int logicalv3[3] = {idx[0], idx[1], idx[2]+1};
      int logicalv4[3] = {idx[0], idx[1]+1, idx[2]};
      int logicalv5[3] = {idx[0]+1, idx[1]+1, idx[2]};
      int logicalv6[3] = {idx[0]+1, idx[1]+1, idx[2]+1};
      int logicalv7[3] = {idx[0], idx[1]+1, idx[2]+1};

      int indexv0 = GetPointIndex(logicalv0, dims);
      int indexv1 = GetPointIndex(logicalv1, dims);
      int indexv2 = GetPointIndex(logicalv2, dims);
      int indexv3 = GetPointIndex(logicalv3, dims);
      int indexv4 = GetPointIndex(logicalv4, dims);
      int indexv5 = GetPointIndex(logicalv5, dims);
      int indexv6 = GetPointIndex(logicalv6, dims);
      int indexv7 = GetPointIndex(logicalv7, dims);

      float f0 = F[indexv0];
      float f1 = F[indexv1];
      float f2 = F[indexv2];
      float f3 = F[indexv3];
      float f4 = F[indexv4];
      float f5 = F[indexv5];
      float f6 = F[indexv6];
      float f7 = F[indexv7];

      int caseId = ((f7>isoval) * 128 + 
                    (f6>isoval) * 64 + 
                    (f5>isoval) * 32 + 
                    (f4>isoval) * 16 + 
                    (f3>isoval) * 8 + 
                    (f2>isoval) * 4 + 
                    (f1>isoval) * 2 + 
                    (f0>isoval) * 1 + 
                    )
      int lup[16] = lookupTable[caseId];
      int k = 0;
      while(edges[k] != -1){
        int edge1 = lup[k];
        int edge2 = lup[k+1];
        int edge3 = lup[k+2];
        int edges[3] = {edge1, edge2, edge3};

        k+=3;
        float pts[3][3];

        for (int b =0; b < 3; b++){
          if (edge[b] == 0){
            pts[b][0] = X[logicalv0[0]] + ((isoval-f0)/(f1-f0)) * (X[logicalv1[0]] - logicalv0[0]);
            pts[b][1] = Y[logicalv0[1]];
            pts[b][2] = Z[logicalv0[2]];
          }
          else if (edge[b] == 1){
            pts[b][0] = X[logicalv1[0]];
            pts[b][1] = Y[logicalv1[1]];
            pts[b][2] = Z[logicalv1[2]] + ((isoval - f1)/(f2-f1)) * (Z[logicalv2[2]] - Z[logicalv1[2]]);
          }
          else if(edge[b] == 2){
            pts[b][0] = X[logicalv3[0]] + ((isoval-f3)/(f2-f3)) * (X[logicalv2[0]] - logicalv3[0]);
            pts[b][1] = Y[logicalv3[1]];
            pts[b][2] = Z[logicalv3[2]];
          }
          else if (edge[b] == 3){
            pts[b][0] = X[logicalv0[0]];
            pts[b][1] = Y[logicalv0[1]];
            pts[b][2] = Z[logicalv0[2]] + ((isoval - f0)/(f3-f0)) * (Z[logicalv3[2]] - Z[logicalv0[2]]);
          }
          else if (edge[b] == 4){
            pts[b][0] = X[logicalv4[0]] + ((isoval-f4)/(f5-f4)) * (X[logicalv5[0]] - logicalv4[0]);
            pts[b][1] = Y[logicalv4[1]];
            pts[b][2] = Z[logicalv4[2]];
          }
          else if (edge[b] == 5){
            pts[b][0] = X[logicalv5[0]];
            pts[b][1] = Y[logicalv5[1]];
            pts[b][2] = Z[logicalv5[2]] + ((isoval - f5)/(f6-f5)) * (Z[logicalv6[2]] - Z[logicalv5[2]]);
          }
          else if (edge[b] == 6){
            pts[b][0] = X[logicalv7[0]] + ((isoval-f7)/(f6-f7)) * (X[logicalv6[0]] - logicalv7[0]);
            pts[b][1] = Y[logicalv7[1]];
            pts[b][2] = Z[logicalv7[2]];
          }
          else if (edge[b] == 7){
            pts[b][0] = X[logicalv4[0]];
            pts[b][1] = Y[logicalv4[1]];
            pts[b][2] = Z[logicalv4[2]] + ((isoval - f4)/(f7-f4)) * (Z[logicalv7[2]] - Z[logicalv4[2]]);
          }
          else if (edge[b] == 8){
            pts[b][0] = X[logicalv0[0]];
            pts[b][1] = Y[logicalv0[1]] + ((isoval-f0)/(f4-f0)) * (Y[logicalv4[1]] - Y[logicalv0[1]]);
            pts[b][2] = Z[logicalv0[2]];
          }
          else if (edge[b] == 9){
            pts[b][0] = X[logicalv1[0]];
            pts[b][1] = Y[logicalv1[1]] + ((isoval-f1)/(f5-f1)) * (Y[logicalv5[1]] - Y[logicalv1[1]]);
            pts[b][2] = Z[logicalv1[2]];
          }
          else if (edge[b] == 10){
            pts[b][0] = X[logicalv3[0]];
            pts[b][1] = Y[logicalv3[1]] + ((isoval-f3)/(f7-f3)) * (Y[logicalv7[1]] - Y[logicalv3[1]]);
            pts[b][2] = Z[logicalv3[2]];
          }
        }
      }
         // YOUR CODE GOES HERE
         // My advice:
         //   (1) collect all of the info about a cell (8 vertex locations, 8 field values) and put them in arrays you can use
         //          -- and then add print statements and debug
         //   (2) figure out which case in the lookup table this cell goes to 
         //          -- and then add print statements and debug
         //   (3) determine the triangles within the cell using the lookup table and info about the cell
         //          -- hints:
         //              check to make sure T is between 0 and 1.  call abort() if not and debug.
         //              you will need to calculate the exact position of each triangle vertex along an edge.
         //              I put my three vertices in "float Xt[3], Yt[3], Zt[3];"
         //              And then called:             tl.AddTriangle(Xt[0], Yt[0], Zt[0], Xt[1], Yt[1], Zt[1], Xt[2], Yt[2], Zt[2]);
         // 
         // My print statements for cell 4771:
/*
	Cell Log Idx = 18, 48, 1
	Cell's 8 points are
		Pt[0] = 18, 48, 1
		 ptIdx = 4918, field = 3.29052, loc = (-2.65306, 9.59184, -9.59184)
		Pt[1] = 19, 48, 1
		 ptIdx = 4919, field = 3.18607, loc = (-2.2449, 9.59184, -9.59184)
		Pt[2] = 19, 48, 2
		 ptIdx = 7419, field = 3.16074, loc = (-2.2449, 9.59184, -9.18367)
		Pt[3] = 18, 48, 2
		 ptIdx = 7418, field = 3.25718, loc = (-2.65306, 9.59184, -9.18367)
		Pt[4] = 18, 49, 1
		 ptIdx = 4968, field = 3.10792, loc = (-2.65306, 10, -9.59184)
		Pt[5] = 19, 49, 1
		 ptIdx = 4969, field = 3.02942, loc = (-2.2449, 10, -9.59184)
		Pt[6] = 19, 49, 2
		 ptIdx = 7469, field = 3.02574, loc = (-2.2449, 10, -9.18367)
		Pt[7] = 18, 49, 2
		 ptIdx = 7468, field = 3.10095, loc = (-2.65306, 10, -9.18367)
	Triangle case index is 9
		 Working on triangle with vertices along edges 0, 10, 8
			 Edge 0 has vertices at 0 and 1, and the new vertex should be placed at t=0.866634
			 Interpolated to get point -2.29933, 9.59184, -9.59184
			 Edge 10 has vertices at 3 and 7, and the new vertex should be placed at t=0.365998
			 Interpolated to get point -2.65306, 9.74123, -9.18367
			 Edge 8 has vertices at 0 and 4, and the new vertex should be placed at t=0.495728
			 Interpolated to get point -2.65306, 9.79418, -9.59184
		 Working on triangle with vertices along edges 10, 2, 0
			 Edge 10 has vertices at 3 and 7, and the new vertex should be placed at t=0.365998
			 Interpolated to get point -2.65306, 9.74123, -9.18367
			 Edge 2 has vertices at 2 and 3, and the new vertex should be placed at t=0.407094
			 Interpolated to get point -2.41106, 9.59184, -9.18367
			 Edge 0 has vertices at 0 and 1, and the new vertex should be placed at t=0.866634
			 Interpolated to get point -2.29933, 9.59184, -9.59184
 */
    }

    vtkPolyData *pd = tl.MakePolyData();

/*
    //This can be useful for debugging
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("proj6_out.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkCleanPolyData *cpd = vtkCleanPolyData::New();
    cpd->SetInputData(pd);
    cpd->SetAbsoluteTolerance(0);
    cpd->PointMergingOn();
    cpd->Update();
    vtkPolyDataNormals *pdn = vtkPolyDataNormals::New();
    pdn->SetInputData(cpd->GetOutput());
    //pdn->SetInputData(pd);
    pdn->Update();

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pdn->GetOutput());
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0., 0., 0.);
    ren1->GetActiveCamera()->SetPosition(0,0,-62);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(1, 100);
    //ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
