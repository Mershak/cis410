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

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


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
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
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
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always zero if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always zero if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)*idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always zero if 2D)
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
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always zero if 2D)
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
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}


// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//           The first number is the size of the array in argument X, 
//           the second the size of Y.
//     X: an array (size is specified by dims).  
//        This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//        This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  
//        Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float
EvaluateFieldAtLocation(const float *pt, const int *dims, 
                        const float *X, const float *Y, const float *F)
{
    int xP = -10;
    int yP = -10;
    for (int i = 0; i < dims[0]; i++){
        if (X[i] > pt[0]){
            xP = i - 1;
            break;
        }
    }

    for (int i = 0; i < dims[1]; i++){
        if (Y[i] > pt[1]){
            yP = i - 1;
            break;
        }
    }

    if (xP == -10 || yP == -10){
        return 0;
    }
    int tLP[] = {xP,(yP+1)};
    int tRP[] = {(xP+1),(yP+1)};
    int bLP[] = {xP,yP};
    int bRP[] = {(xP+1),yP};
    float topLeft = F[GetPointIndex(tLP,dims)];
    float topRight = F[GetPointIndex(tRP,dims)];
    float bottomLeft = F[GetPointIndex(bLP,dims)];
    float bottomRight = F[GetPointIndex(bRP,dims)];

    float topInterpolate = topLeft + ((pt[0] - X[xP]) / (X[xP+1] - X[xP])) *  (topRight - topLeft);
    float bottomInterpolate = bottomLeft + ((pt[0] - X[xP]) / (X[xP+1] - X[xP])) *  (bottomRight - bottomLeft);
    // float pointInterpolate = (float) malloc(sizeof(float));
    float pointInterpolate = bottomInterpolate + ((pt[1] - Y[yP]) / (Y[yP+1] - Y[yP])) * (topInterpolate - bottomInterpolate);

    return pointInterpolate; // IMPLEMENT ME!
}


// ****************************************************************************
//  Function: AreaForCell
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     cellId: a cellIndex (I.e., between 0 and GetNumberOfCells(dims))
//  
//  Returns:  the area of the cell.  Each cell is a rectangle, and the
//            area of a rectangle is width*height.
//            If an invalid cell is specified, then the function should return 0.
//
// ****************************************************************************

float
AreaForCell(const float *X, const float *Y, const int *dims, int cellId)
{
    // IMPLEMENT ME!
    if (cellId > GetNumberOfCells(dims) || cellId < 0){
        return 0;
    }
    int idx[2];
    GetLogicalCellIndex(idx,cellId,dims);
    int posX = idx[0];
    int posY = idx[1];
    
    float difY = Y[posY+1] - Y[posY];
    float difX = X[posX+1] - X[posX];

    // float result = (float) malloc(sizeof(float));
    float result = difX * difY;
    return result;
}

// ****************************************************************************
//  Function: CountNumberOfCellsWithinThreshold
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//     min: the minimum range for a threshold
//     max: the maximum range for a threshold
//
//  Returns:  the number of cells within a threshold range.
//
//  Example: assume the desired threshold range is (L, H).  Consider cell C,
//    which has 4 vertices, and assume these 4 vertices have values F1, F2,
//    F3, and F4.
//    Then we define C to be within the threshold range if:
//       L <= F1 <= H  *AND*
//       L <= F2 <= H  *AND*
//       L <= F3 <= H  *AND*
//       L <= F4 <= H 
//
//   Your goal is to count the number of cells that are within the threshold
//   range. 
//
// ****************************************************************************

int
CountNumberOfCellsWithinThreshold(const float *X, const float *Y, const int *dims,
                             const float *F, float L, float H)
{   
    int result = 0;
    for (int i = 0; i <GetNumberOfCells(dims); i++){
        int idx[2];
        GetLogicalCellIndex(idx,i,dims);

        int tLP[] = {idx[0],idx[1]+1};
        int tRP[] = {idx[0]+1,idx[1]+1};
        int bLP[] = {idx[0],idx[1]};
        int bRP[] = {idx[0]+1,idx[1]};

        float fTopLeft = F[GetPointIndex(tLP,dims)];
        float fTopRight = F[GetPointIndex(tRP,dims)];
        float fBottomLeft = F[GetPointIndex(bLP,dims)];
        float fBottomRight = F[GetPointIndex(bRP,dims)];
        if (fTopLeft <=H && fTopLeft >= L 
        && fTopRight <=H && fTopRight >= L 
        && fBottomLeft <=H && fBottomLeft >= L 
        &&fBottomRight <=H && fBottomRight >= L){
            result += 1;
        }
    }
    return result;
    // IMPLEMENT ME!
}

int main()
{
    int  i;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj2_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    float ranges[10] = { -1, 1, -1, 0, 0, 1, -0.75, -0.25, 0.25, 0.75 };
    for (int i = 0 ; i < 5 ; i++)
    {
        int numCells = CountNumberOfCellsWithinThreshold(X, Y, dims, F, ranges[2*i], ranges[2*i+1]);
        cerr << "The number of cells between " << ranges[2*i] << " and " << ranges[2*i+1] << " is " << numCells << endl;
    }

    const int ncells = 5;
    int cellIds[ncells] = { 0, 50, 678, 1000, 1200 };
    for (i = 0 ; i < ncells ; i++)
    {
        float area = AreaForCell(X, Y, dims, cellIds[i]);
        cerr << "The area for cell " << cellIds[i] << " is " << area << endl;
    }

    const int npts = 10;
    float pt[npts][3] = 
         {
            {1.01119, 0.122062, 0},
            {0.862376, 1.33829, 0},
            {0.155026, 0.126123, 0},
            {0.69736, 0.0653565, 0},
            {0.2, 0.274117, 0},
            {0.893699, 1.04111, 0},
            {0.608791, -0.0533753, 0},
            {1.00543, 0.158024, 0},
            {0.384158, -0.0768977, 0},
            {0.666757, 0.60259, 0},
         };

    

    for (i = 0 ; i < npts ; i++)
    {
        float f = EvaluateFieldAtLocation(pt[i], dims, X, Y, F);
        cerr << "Evaluated field at (" << pt[i][0] <<"," << pt[i][1] << ") as "
             << f << endl;
    }
    
   
    cerr << "Infinite loop here, else Windows people may have the terminal "
         << "disappear before they see the output."
         << " Remove these lines if they annoy you." << endl;
    cerr << "(press Ctrl-C to exit program)" << endl;
    while (1) ; 
}




