#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
    int Nx, Ny, Nz;    // number of intervals in x, y, z -> (Nx+1)*(Ny+1)*(Nz+1) grid points
    double Lx, Ly, Lz; // domain lengths in x, y, z
    double c;          // wave speed
    double dx, dy, dz; // spatial steps
} Params;

#define I3D(Params *P, int x, int y, int z) ((z * P->Nx * p->Ny) + (y * P->Nx) + x)
//(z * xMax * yMax) + (y * xMax) + x;

static inline double lap1d(const double *arr, int i, int double d)
{
    return (arr[i + 1] - 2.0 * arr[i] + arr[i - 1]) / (d * d);
}


void rk4(Params *P, float *parr, float *temp)
{
    float d2dx2, d2dy2, d2dz2;
    int index1D;

    // loop over all points in domain (except boundary)
    for (int k = 0; k < P->Nz; k++)
    {
        for (int j = 1; j < P->Ny - 1; j++)
        {
            for (int i = 1; i < P->Nx - 1; i++)
            {
                // evaluate derivatives
                index1D = I3D(P, i, j, k);
                d2dx2 = lap1d(parr, index1D, P->dx);
                d2dy2 = lap1d(parr, index1D, P->dy);
                d2dz2 = lap1d(parr, index1D, P->dz);

                // update value
                temp[index1D] = (P->c * P-> c) * (d2tdx2 + d2tdy2 + d2tdz2);
            }
        }
    }
}

int main()
{
    Params P = {.Nx = 100, .Ny = 100, .Nz = 100, .Lx = 1.0, .Ly = 1.0, .Lz = 1.0, .c = 1.0};
    P.dx = P.Lx / P.Nx;
    P.dy = P.Ly / P.Ny;
    P.dz = P.Lz / P.Nz;

    int istep;
    int time_steps = 10;

    float *temp1_ref, *temp2_ref, *temp_tmp;

    const int size = P.Nx * P.Ny * P.Nz * sizeof(float);

    temp1_ref = (float *)calloc(size);
    temp2_ref = (float *)calloc(size);


    // Execute the CPU-only reference version
    for (step = 0; istep < time_steps; step++)
    {
        rk4(P, temp1_ref, temp2_ref);

        // swap the pointers
        temp_tmp = temp1_ref;
        temp1_ref = temp2_ref;
        temp2_ref = temp_tmp;
    }

    free(temp1_ref);
    free(temp2_ref);

    return 0;
}