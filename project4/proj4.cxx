#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSphereSource.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

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
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
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
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}


// ****************************************************************************
//  Function: EvaluateVectorFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//        The first value in the field is the x-component for the first point.
//        The second value in the field is the y-component for the first point.
//
//     rv (output): the interpolated field value. (0,0) if the location is out of bounds.
//
// ****************************************************************************

void EvaluateVectorFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F, float *rv)
{
    // IMPLEMENT ME!
    

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

    int tLP[] = {xP,(yP+1)};
    int tRP[] = {(xP+1),(yP+1)};
    int bLP[] = {xP,yP};
    int bRP[] = {(xP+1),yP};
    float topLeftX = F[2 * (GetPointIndex(tLP,dims))];
    float topLeftY = F[2 * (GetPointIndex(tLP,dims)) + 1];
    float topRightX = F[2 * (GetPointIndex(tRP,dims))];
    float topRightY = F[2 * (GetPointIndex(tRP,dims)) + 1];
    float bottomLeftX = F[2 * (GetPointIndex(bLP,dims))];
    float bottomLeftY = F[2 * (GetPointIndex(bLP,dims)) + 1];
    float bottomRightX = F[2 * (GetPointIndex(bRP,dims))];
    float bottomRightY = F[2 * (GetPointIndex(bRP,dims)) + 1];

    float topInterpolateX = topLeftX + ((pt[0] - X[xP]) / (X[xP+1] - X[xP])) *  (topRightX - topLeftX);
    float bottomInterpolateX = bottomLeftX + ((pt[0] - X[xP]) / (X[xP+1] - X[xP])) *  (bottomRightX - bottomLeftX);

    float pointInterpolateX = bottomInterpolateX + ((pt[1] - Y[yP]) / (Y[yP+1] - Y[yP])) * (topInterpolateX - bottomInterpolateX);

    float topInterpolateY = topLeftY + ((pt[0] - X[xP]) / (X[xP+1] - X[xP])) *  (topRightY - topLeftY);
    float bottomInterpolateY = bottomLeftY + ((pt[0] - X[xP]) / (X[xP+1] - X[xP])) *  (bottomRightY - bottomLeftY);

    float pointInterpolateY = bottomInterpolateY + ((pt[1] - Y[yP]) / (Y[yP+1] - Y[yP])) * (topInterpolateY - bottomInterpolateY);

    rv[0] = pointInterpolateX; // setting the x-component of the velocity
    rv[1] = pointInterpolateY; // setting the y-component of the velocity
}

// ****************************************************************************
//  Function: AdvectWithEulerStep
//
//  Arguments:
//     pt: the seed location (two-dimensions)
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//     h: The size of the Euler step
//     nsteps: The number of Euler steps to take
//     output_locations (output): An array of size 2*(nsteps+1).  It's first entry
//        should be the seed location.  The second entry should be the result
//        of the first advection step.  The final entry should be the result
//        of the final advection step.
//
// ****************************************************************************

void
AdvectWithEulerStep(const float *pt, const int *dims, const float *X, 
                    const float *Y, const float *F, 
                    float h, int nsteps, float *output_locations)
{
    // IMPLEMENT ME!
    output_locations[0] = pt[0];
    output_locations[1] = pt[1];

    float lastPos[2] = { pt[0], pt[1]};

    for (int i = 1; i <= nsteps; i++){
        float velocity[2];
        EvaluateVectorFieldAtLocation(lastPos, dims, X, Y, F, velocity);

        float nextPos[2];
        nextPos[0] = lastPos[0] + h * velocity[0];
        nextPos[1] = lastPos[1] + h * velocity[1];
        output_locations[2 * i] = nextPos[0];
        output_locations[2 * i  + 1] = nextPos[1];
        lastPos[0] = nextPos[0];
        lastPos[1] = nextPos[1];
    }
}

void
AdvectWithRK4Step(const float *seed, const int *dims, const float*X, 
                const float *Y, const float*F, float h, int nsteps, float *RK4_output_locations)
{

    RK4_output_locations[0] = seed[0];
    RK4_output_locations[1] = seed[1];
    float lastPos[2] = {seed[0], seed[1]};

    for (int i = 1; i <= nsteps; i++){
        float vel1[2];
        float k1Pos[2] = {lastPos[0], lastPos[1]};
        EvaluateVectorFieldAtLocation(k1Pos, dims, X, Y, F, vel1);

        float vel2[2];
        float k2Pos[2] = {(lastPos[0] + h/2*vel1[0]), (lastPos[1] + h/2*vel1[1])};
        EvaluateVectorFieldAtLocation(k2Pos, dims, X, Y, F, vel2);

        float vel3[2];
        float k3Pos[2] = {(lastPos[0] + h/2*vel2[0]), (lastPos[1] + h/2*vel2[1])};
        EvaluateVectorFieldAtLocation(k3Pos, dims, X, Y, F, vel3);

        float vel4[2];
        float k4Pos[2] = {(lastPos[0] + h*vel3[0]), (lastPos[1] + h*vel3[1])};
        EvaluateVectorFieldAtLocation(k4Pos, dims, X, Y, F, vel4);


        float k1X = vel1[0];
        float k2X = vel2[0];
        float k3X = vel3[0];
        float k4X = vel4[0];
        float k1Y = vel1[1];
        float k2Y = vel2[1];
        float k3Y = vel3[1];
        float k4Y = vel4[1];

        float newPos[2];
        newPos[0] = lastPos[0] + h * (k1X + 2*k2X + 2* k3X + k4X)/6;
        newPos[1] = lastPos[1] + h * (k1Y + 2*k2Y + 2* k3Y + k4Y)/6;
        RK4_output_locations[2*i] = newPos[0];
        RK4_output_locations[2*i+1] = newPos[1];
        lastPos[0] = newPos[0];
        lastPos[1] = newPos[1];
    }

}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// VTK files are only 3D, so the vector data is all of the form (X,Y,0).
// Remove the 0's since it is counter-intuitive for students who are 
// thinking of this as 2D data.
float *
Convert3DVectorDataTo2DVectorData(const int *dims, const float *F)
{
    float *rv = new float[dims[0]*dims[1]*2];
    int index3D = 0;
    int index2D = 0;
    for (int i = 0 ; i < dims[0] ; i++)
       for (int j = 0 ; j < dims[1] ; j++)
       {
           rv[index2D]   = F[index3D];
           rv[index2D+1] = F[index3D+1];
           index2D += 2;
           index3D += 3;
       }

    return rv;
}

void PrintSteps(const char *solver, int nsteps, float *locations)
{
   cerr << "Printing output for solver " << solver << endl;
   for (int j = 0 ; j < nsteps+1 ; j++)
   {
       cerr << j << ": (" << locations[2*j] << ", " << locations[2*j+1] << ")" << endl;
   }
}

int main()
{
    int  i, j;

    // HANK'S CODE TO SET THINGS UP -- DON'T MODIFY THIS
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj4_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    if (dims[0] <= 0 || dims[1] <= 0)
    {
        cerr << "Was not able to successfully open file \"proj4_data.vtk\"" << endl;
        exit(EXIT_FAILURE);
    }
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F_3D = (float *) rgrid->GetPointData()->GetVectors()->GetVoidPointer(0);
    float *F = Convert3DVectorDataTo2DVectorData(dims, F_3D);
    
    float seed[2] = { 1, -5 };
    
    // SANITY CHECK TO MAKE SURE VECTOR FIELD EVALUATION IS WORKING
    float vec[2];
    EvaluateVectorFieldAtLocation(seed, dims, X, Y, F, vec);
    cerr << "Velocity at (" << seed[0] <<", " << seed[1] << ") is (" << vec[0] << ", " << vec[1] << ")" << endl;

    float h = 0.01;
    const int nsteps = 100;
    float *output_locations = new float[2*(nsteps+1)];
    AdvectWithEulerStep(seed, dims, X, Y, F, h, nsteps, output_locations);
    PrintSteps("Euler", nsteps, output_locations);



    float *RK4_output_locations = new float[2*(nsteps+1)];
    AdvectWithRK4Step(seed, dims, X, Y, F, h, nsteps, RK4_output_locations);
    PrintSteps("RK4", nsteps, RK4_output_locations);
    
}
