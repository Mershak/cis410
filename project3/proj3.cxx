#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
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
//  Function: EvaluateFieldAtLocation
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
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F)
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

    return pointInterpolate;
}


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
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

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,0.5) 
//        F=1: (1.0,1.0,1.0) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
//  Note: color map is in float (0 to 1) and output is in unsigned char
//        (0 to 255), so the floating point values you calculate for 
//        the three color channels need to be multiplied by 255 when
//        being assigned to RGB.
//
// ****************************************************************************

void ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    RGB[0] = F*255;
    RGB[1] = F*255;
    RGB[2] = 0.5*255 + F * 0.5*255;
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color 
//     using a difference colormap.
//
//     The difference color map has:
//        F=0: (0,0,0.5) 0,0,127
//        F=0.5: (1.0,1.0,1.0) 255, 255 ,255
//        F=1: (0.5, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
//  Note: color map is in float (0 to 1) and output is in unsigned char
//        (0 to 255), so the floating point values you calculate for 
//        the three color channels need to be multiplied by 255 when
//        being assigned to RGB.
//
// ****************************************************************************

void ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    if (F <= 0.5){
        RGB[0] = 255.0*(2.0*F);
        RGB[1] = 255.0*(2.0*F);
        RGB[2] = 255.0/2.0 + 255.0*F;
    }
    else{
        RGB[0] = 255.0/2 + (1.0-F)*225.0;
        RGB[1] = 0.0 + (1.0-F)*2.0*255.0;
        RGB[2] = 0.0 + (1.0-F)*2.0*255.0;
        // RGB[0] = 255.0 - (F-0.5)*255.0;
        // RGB[1] = 255.0 - (F-0.5)*2.0*255.0;
        // RGB[2] = 255.0 - (F-0.5)*2.0*225.0;
    }
}

// ****************************************************************************
//  Function: ApplyHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0 <= F <= 1) to a color using 
//     an HSV rainbow colormap.
//
//     The rainbow colormap uses a saturation = 1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees.
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1  
//       RGB (output):  the location to store the color
//      
//  Note: as with the first two functions, make sure to multiple by 255 
//        when converting floats to unsigned chars.
//
// ****************************************************************************

void ApplyHSVColorMap(float F, unsigned char *RGB)  
{   
    float r;
    float g;
    float b;

    float saturation = 1.0;
    float value = 1.0;
    float v = 1.0;
    float hue = (360.0 * F);
    if (hue == 360.0){
        hue = 0.0;
    }
    hue /= 60.f;
    // sector 0 to 5  
    int i = floor( hue );
    float f = hue - i;
    // factorial part of h
    float p = value * ( 1 - saturation);
    float q = value * ( 1 - saturation * f );
    float t = value * ( 1 - saturation * ( 1 - f ) );      
switch (i)
		{
		case 0:
			r = v;
			g = t;
			b = p;
			break;

		case 1:
			r = q;
			g = v;
			b = p;
			break;

		case 2:
			r = p;
			g = v;
			b = t;
			break;

		case 3:
			r = p;
			g = q;
			b = v;
			break;

		case 4:
			r = t;
			g = p;
			b = v;
			break;

		default:
			r = v;
			g = p;
			b = q;
			break;
		}
    RGB[0] = r*255;
    RGB[1] = g*255;
    RGB[2] = b*255;
}


int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3 ; i++)
       for (j = 0 ; j < 3*nx*ny ; j++)
            buffer[i][j] = 0;

    for (i = 0 ; i < nx ; i++)
        for (j = 0 ; j < ny ; j++)
        {
            // ITERATE OVER PIXELS
            float pt[2];
            pt[0] = -9.0 + i*(18.0/500.0);
            pt[1] = -9.0 + j*(18.0/500.0);

            float f = EvaluateFieldAtLocation(pt,dims,X,Y,F);
            float normalizedF = (f-1.2)/ 3.8;
            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
