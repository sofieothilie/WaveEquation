#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
    int Nx, Ny, Nz;
    double Lx, Ly, Lz;
    double c;
    double dx, dy, dz;
} Params;

static inline int get_index(const Params *P, int x, int y, int z)
{
    return (z * P->Ny * P->Nx) + (y * P->Nx) + x;
}

// Laplacian at (i,j,k)
double laplacian(const Params *P, float *arr, int i, int j, int k)
{
    int idx = get_index(P, i, j, k);
    double lap = 0.0;
    // Second-order central differences for 3D Laplacian
    // d2u/dx2
    lap += (arr[get_index(P, i + 1, j, k)] - 2 * arr[idx] + arr[get_index(P, i - 1, j, k)]) / (P->dx * P->dx);
    // + d2u/dy2
    lap += (arr[get_index(P, i, j + 1, k)] - 2 * arr[idx] + arr[get_index(P, i, j - 1, k)]) / (P->dy * P->dy);
    // + d2u/dz2
    lap += (arr[get_index(P, i, j, k + 1)] - 2 * arr[idx] + arr[get_index(P, i, j, k - 1)]) / (P->dz * P->dz);
    return lap;
}

// Initial condition: Spherical outgoing wave
void set_initial_conditions(const Params *P, float *u, float *v)
{
    double c = P->c;
    double eps = 1e-6; // avoid dividing by 0 in centre.
    for (int k = 0; k < P->Nz; k++)
    {
        double z = k * P->dz - 0.5 * P->Lz;
        for (int j = 0; j < P->Ny; j++)
        {
            double y = j * P->dy - 0.5 * P->Ly;
            for (int i = 0; i < P->Nx; i++)
            {
                double x = i * P->dx - 0.5 * P->Lx;
                int idx = get_index(P, i, j, k);
                double r = sqrt(x * x + y * y + z * z);
                u[idx] = sin(10.0 * r) / (r + eps);
                v[idx] = -10.0 * c * cos(10.0 * r) / (r + eps); // du/dt at t=0
            }
        }
    }
}

// Write 3D data to file
void write_3d_data(const Params *P, float *arr, const char *filename)
{
    FILE *f = fopen(filename, "w");
    if (!f)
    {
        perror("Failed to open file for writing");
        return;
    }
    for (int k = 0; k < P->Nz; k++)
    {
        double z = k * P->dz - 0.5 * P->Lz;
        for (int j = 0; j < P->Ny; j++)
        {
            double y = j * P->dy - 0.5 * P->Ly;
            for (int i = 0; i < P->Nx; i++)
            {
                double x = i * P->dx - 0.5 * P->Lx;
                int idx = get_index(P, i, j, k);
                fprintf(f, "%f %f %f %f\n", x, y, z, arr[idx]);
            }
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

// RK4 step for the wave equation
void rk4_wave(const Params *P, float *u, float *v, float *u_tmp, float *v_tmp, double dt)
{
    int N = P->Nx * P->Ny * P->Nz;
    float *ku1 = (float *)calloc(N, sizeof(float));
    float *kv1 = (float *)calloc(N, sizeof(float));
    float *ku2 = (float *)calloc(N, sizeof(float));
    float *kv2 = (float *)calloc(N, sizeof(float));
    float *ku3 = (float *)calloc(N, sizeof(float));
    float *kv3 = (float *)calloc(N, sizeof(float));
    float *ku4 = (float *)calloc(N, sizeof(float));
    float *kv4 = (float *)calloc(N, sizeof(float));

    // k1
    for (int k = 1; k < P->Nz - 1; k++)
    {
        for (int j = 1; j < P->Ny - 1; j++)
        {
            for (int i = 1; i < P->Nx - 1; i++)
            {
                int idx = get_index(P, i, j, k);
                ku1[idx] = v[idx];
                kv1[idx] = P->c * P->c * laplacian(P, u, i, j, k);
            }
        }
    }
    // k2
    for (int idx = 0; idx < N; idx++)
    {
        u_tmp[idx] = u[idx] + 0.5 * dt * ku1[idx];
        v_tmp[idx] = v[idx] + 0.5 * dt * kv1[idx];
    }
    for (int k = 1; k < P->Nz - 1; k++)
    {
        for (int j = 1; j < P->Ny - 1; j++)
        {
            for (int i = 1; i < P->Nx - 1; i++)
            {
                int idx = get_index(P, i, j, k);
                ku2[idx] = v_tmp[idx];
                kv2[idx] = P->c * P->c * laplacian(P, u_tmp, i, j, k);
            }
        }
    }
    // k3
    for (int idx = 0; idx < N; idx++)
    {
        u_tmp[idx] = u[idx] + 0.5 * dt * ku2[idx];
        v_tmp[idx] = v[idx] + 0.5 * dt * kv2[idx];
    }
    for (int k = 1; k < P->Nz - 1; k++)
    {
        for (int j = 1; j < P->Ny - 1; j++)
        {
            for (int i = 1; i < P->Nx - 1; i++)
            {
                int idx = get_index(P, i, j, k);
                ku3[idx] = v_tmp[idx];
                kv3[idx] = P->c * P->c * laplacian(P, u_tmp, i, j, k);
            }
        }
    }
    // k4
    for (int idx = 0; idx < N; idx++)
    {
        u_tmp[idx] = u[idx] + dt * ku3[idx];
        v_tmp[idx] = v[idx] + dt * kv3[idx];
    }
    for (int k = 1; k < P->Nz - 1; k++)
    {
        for (int j = 1; j < P->Ny - 1; j++)
        {
            for (int i = 1; i < P->Nx - 1; i++)
            {
                int idx = get_index(P, i, j, k);
                ku4[idx] = v_tmp[idx];
                kv4[idx] = P->c * P->c * laplacian(P, u_tmp, i, j, k);
            }
        }
    }
    // Update u and v
    for (int idx = 0; idx < N; idx++)
    {
        u[idx] += (dt / 6.0) * (ku1[idx] + 2 * ku2[idx] + 2 * ku3[idx] + ku4[idx]);
        v[idx] += (dt / 6.0) * (kv1[idx] + 2 * kv2[idx] + 2 * kv3[idx] + kv4[idx]);
    }

    free(ku1);
    free(kv1);
    free(ku2);
    free(kv2);
    free(ku3);
    free(kv3);
    free(ku4);
    free(kv4);
}

int main()
{
    Params P = {
        .Nx = 30,
        .Ny = 30,
        .Nz = 30,
        .Lx = 2.0,
        .Ly = 2.0,
        .Lz = 2.0,
        .c = 1.0};
    P.dx = P.Lx / (P.Nx - 1);
    P.dy = P.Ly / (P.Ny - 1);
    P.dz = P.Lz / (P.Nz - 1);

    int steps = 20;
    double dt = 1.0 / steps;

    // stability check
    double cmax = sqrt(2);
    double h = sqrt(1/(P.dx*P.dx) + 1/(P.dy*P.dy) + 1/(P.dz*P.dz));


    if ((P.c * dt) / h > cmax)
    {
        printf("Unstable simulation. Exiting. \n");
        return 0;
    }

    int N = P.Nx * P.Ny * P.Nz;
    float *u = (float *)calloc(N, sizeof(float));
    float *v = (float *)calloc(N, sizeof(float));
    float *u_tmp = (float *)calloc(N, sizeof(float));
    float *v_tmp = (float *)calloc(N, sizeof(float));

    set_initial_conditions(&P, u, v);

    for (int s = 0; s < steps; s++)
    {
        char fname[64];
        snprintf(fname, sizeof(fname), "timestep/wave3d_%03d.txt", s);
        write_3d_data(&P, u, fname);
        rk4_wave(&P, u, v, u_tmp, v_tmp, dt);
    }

    char fname[64];
    snprintf(fname, sizeof(fname), "timestep/wave3d_%03d.txt", steps);
    write_3d_data(&P, u, fname);

    free(u);
    free(v);
    free(u_tmp);
    free(v_tmp);
    return 0;
}