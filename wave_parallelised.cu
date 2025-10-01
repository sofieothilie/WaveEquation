#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dimx 8
#define dimy 8
#define dimz 4

typedef struct
{
    int Nx, Ny, Nz;
    double Lx, Ly, Lz;
    double c;
    double dx, dy, dz;
} Params;

static inline __host__ __device__ int get_index(const Params *P, int x, int y, int z)
{
    // rewrite to fit with threadid and stuff.
    return (z * P->Ny * P->Nx) + (y * P->Nx) + x;
}

// Laplacian at (i,j,k)
__device__ double laplacian(const Params *P, float *arr, int i, int j, int k)
{
    if (i <= 0 || i >= P->Nx - 1 ||
        j <= 0 || j >= P->Ny - 1 ||
        k <= 0 || k >= P->Nz - 1)
        return 0.0;
    
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
    // vectorise
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

static inline __device__ int get_coords(const Params *P, int *x_p, int *y_p, int *z_p)
{
    // Compute the offset in each dimension
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    if (x >= P->Nx || y >= P->Ny || z >= P->Nz)
        return 0;

    *x_p = x;
    *y_p = y;
    *z_p = z;
    return 1;
}

__global__ void vec_add_uv(const Params *P, float *ku, float *kv, float *u, float *v)
{
    int x, y, z;
    printf("hello\n");

    if (get_coords(P, &x, &y, &z) == 0)
    {
        printf("out of bounds");
        return; // Thread out of bounds
    }
    int idx = get_index(P, x, y, z);
    printf("%i", idx);
    ku[idx] = v[idx];
    printf("inside gpu: ku[%i]: %f\n", idx, ku[idx]);
    kv[idx] = P->c * P->c * laplacian(P, u, x, y, z);
    printf("inside gpu: kv[%i]: %f\n", idx, kv[idx]);
}

__global__ void vec_add_rk(const Params *P, float *u, float *v, float *u_tmp, float *v_tmp, double dt, float *ku, float *kv)
{
    int x, y, z;
    get_coords(P, &x, &y, &z);
    int idx = get_index(P, x, y, z);
    u_tmp[idx] = u[idx] + dt * ku[idx];
    v_tmp[idx] = v[idx] + dt * kv[idx];
}

__global__ void vec_add_frk(const Params *P, float *u, float *v, double dt, float *ku1, float *kv1, float *ku2, float *kv2, float *ku3, float *kv3, float *ku4, float *kv4)
{
    int x, y, z;
    get_coords(P, &x, &y, &z);
    int idx = get_index(P, x, y, z);
    u[idx] += (dt / 6.0) * (ku1[idx] + 2 * ku2[idx] + 2 * ku3[idx] + ku4[idx]);
    v[idx] += (dt / 6.0) * (kv1[idx] + 2 * kv2[idx] + 2 * kv3[idx] + kv4[idx]);
}

// RK4 step for the wave equation
void rk4_wave_parallelised(const Params *P, float *u, float *v, float *u_tmp, float *v_tmp, double dt)
{
    dim3 block(dimx, dimy, dimz);
    dim3 grid((P->Nx + dimx - 1) / dimx, (P->Ny + dimy - 1) / dimy, (P->Nz + dimz - 1) / dimz);
    int N = P->Nx * P->Ny * P->Nz;
    int size = N * sizeof(float);
    float *ku1;
    cudaMalloc(&ku1, size);
    float *kv1;
    cudaMalloc(&kv1, size);
    float *ku2;
    cudaMalloc(&ku2, size);
    float *kv2;
    cudaMalloc(&kv2, size);
    float *ku3;
    cudaMalloc(&ku3, size);
    float *kv3;
    cudaMalloc(&kv3, size);
    float *ku4;
    cudaMalloc(&ku4, size);
    float *kv4;
    cudaMalloc(&kv4, size);

    // k1
    vec_add_uv<<<grid, block>>>(P, ku1, kv1, u, v);

    // k2
    vec_add_rk<<<grid, block>>>(P, u, v, u_tmp, v_tmp, 0.5 * dt, ku1, kv1);
    vec_add_uv<<<grid, block>>>(P, ku2, kv2, u_tmp, v_tmp);

    // k3
    vec_add_rk<<<grid, block>>>(P, u, v, u_tmp, v_tmp, 0.5 * dt, ku2, kv2);
    vec_add_uv<<<grid, block>>>(P, ku3, kv3, u_tmp, v_tmp);

    // k4
    vec_add_rk<<<grid, block>>>(P, u, v, u_tmp, v_tmp, dt, ku3, kv3);
    vec_add_uv<<<grid, block>>>(P, ku4, kv4, u_tmp, v_tmp);

    // Update u and v
    vec_add_frk<<<grid, block>>>(P, u, v, dt, ku1, kv1, ku2, kv2, ku3, kv3, ku4, kv4);

    cudaFree(ku1);
    cudaFree(kv1);
    cudaFree(ku2);
    cudaFree(kv2);
    cudaFree(ku3);
    cudaFree(kv3);
    cudaFree(ku4);
    cudaFree(kv4);
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
    double dt = 1.0 / steps; // must be changed

    // stability check
    double cmax = sqrt(2);
    double h = sqrt(1 / (P.dx * P.dx) + 1 / (P.dy * P.dy) + 1 / (P.dz * P.dz));

    if ((P.c * dt) / h > cmax)
    {
        printf("Unstable simulation. Exiting. \n");
        return 0;
    }

    Params *d_P;
    cudaMalloc(&d_P, sizeof(Params));
    cudaMemcpy(d_P, &P, sizeof(Params), cudaMemcpyHostToDevice);

    int N = P.Nx * P.Ny * P.Nz;
    int size = N * sizeof(float);
    float *h_u = (float *)calloc(N, sizeof(float));
    float *h_v = (float *)calloc(N, sizeof(float));

    set_initial_conditions(&P, h_u, h_v);

    float *d_u;
    cudaMalloc(&d_u, size);
    float *d_v;
    cudaMalloc(&d_v, size);

    cudaMemcpy(d_u, h_u, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, h_v, size, cudaMemcpyHostToDevice);

    float *d_u_tmp;
    cudaMalloc(&d_u_tmp, size);
    float *d_v_tmp;
    cudaMalloc(&d_v_tmp, size);

    system("mkdir -p timestep");

    for (int s = 0; s < steps; s++)
    {
        // save results. Done before in order to also save init
        char fname[64];
        snprintf(fname, sizeof(fname), "timestep/wave3d_%03d.txt", s);
        write_3d_data(&P, h_u, fname);

        // note, last timestep isn't written to file lol.
        rk4_wave_parallelised(d_P, d_u, d_v, d_u_tmp, d_v_tmp, dt);
        cudaDeviceSynchronize();

        cudaMemcpy(h_u, d_u, size, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_v, d_v, size, cudaMemcpyDeviceToHost);
    }

    // save final timestep
    char fname[64];
    snprintf(fname, sizeof(fname), "timestep/wave3d_%03d.txt", steps);
    write_3d_data(&P, h_u, fname);

    // free host memory
    free(h_u);
    free(h_v);

    // free device memory
    cudaFree(d_u);
    cudaFree(d_v);
    cudaFree(d_u_tmp);
    cudaFree(d_v_tmp);

    return 0;
}