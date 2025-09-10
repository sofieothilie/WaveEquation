// rk4_wave1d.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int    N;      // number of intervals -> N+1 grid points
    double L;      // domain length
    double c;      // wave speed
    double dx;     // spatial step
} Params;

// ---------- Problem setup (edit these to your liking) ----------
static inline double f0(double x, double L) { // initial displacement u(x,0)
    return sin(M_PI * x / L);                 // try a half-sine
}

static inline double g0(double x) {           // initial velocity v(x,0) = u_t(x,0)
    return 0.0;
}

// Apply Dirichlet BCs: fixed ends u(0,t)=u(L,t)=0
// TODO 
static inline void apply_bc(double *u, double *v, const Params *P) {
    (void)v;               // velocity BC not strictly needed for Dirichlet u
    u[0] = 0.0;
    u[P->N] = 0.0;
}

// Discrete Laplacian at i using centred differences
// TODO Make it 3D
static inline double lap1d(const double *u, int i, const Params *P) {
    return (u[i+1] - 2.0*u[i] + u[i-1]) / (P->dx * P->dx);
}

// Right-hand side of first-order system:
// u_t = v
// v_t = c^2 * u_xx
static void rhs(const double *u, const double *v,
                double *du, double *dv, const Params *P)
{
    // interior
    for (int i = 1; i < P->N; ++i) {
        du[i] = v[i];
        dv[i] = (P->c * P->c) * lap1d(u, i, P);
    }
    // boundaries: choose something consistent with apply_bc
    du[0]   = 0.0; dv[0]   = 0.0;
    du[P->N]= 0.0; dv[P->N]= 0.0;
}

// One RK4 step: (u,v) -> (u,v) + dt * Φ
static void rk4_step(double *u, double *v, double dt,
                     double *k1u, double *k1v,
                     double *k2u, double *k2v,
                     double *k3u, double *k3v,
                     double *k4u, double *k4v,
                     double *utmp, double *vtmp,
                     const Params *P)
{
    // k1
    rhs(u, v, k1u, k1v, P);

    // k2: evaluate at (u + dt/2 k1u, v + dt/2 k1v)
    for (int i = 0; i <= P->N; ++i) {
        utmp[i] = u[i] + 0.5 * dt * k1u[i];
        vtmp[i] = v[i] + 0.5 * dt * k1v[i];
    }
    apply_bc(utmp, vtmp, P);
    rhs(utmp, vtmp, k2u, k2v, P);

    // k3
    for (int i = 0; i <= P->N; ++i) {
        utmp[i] = u[i] + 0.5 * dt * k2u[i];
        vtmp[i] = v[i] + 0.5 * dt * k2v[i];
    }
    apply_bc(utmp, vtmp, P);
    rhs(utmp, vtmp, k3u, k3v, P);

    // k4
    for (int i = 0; i <= P->N; ++i) {
        utmp[i] = u[i] + dt * k3u[i];
        vtmp[i] = v[i] + dt * k3v[i];
    }
    apply_bc(utmp, vtmp, P);
    rhs(utmp, vtmp, k4u, k4v, P);

    // combine
    for (int i = 0; i <= P->N; ++i) {
        u[i] += (dt/6.0) * (k1u[i] + 2.0*k2u[i] + 2.0*k3u[i] + k4u[i]);
        v[i] += (dt/6.0) * (k1v[i] + 2.0*k2v[i] + 2.0*k3v[i] + k4v[i]);
    }
    apply_bc(u, v, P);
}

int main(void) {
    Params P = { .N = 200, .L = 1.0, .c = 1.0 };
    P.dx = P.L / P.N;

    const double safety = 0.9;
    const double dt = safety * P.dx / P.c;  // CFL-ish for 1-D
    const double Tfinal = 2.0;
    const int steps = (int)ceil(Tfinal / dt);

    // allocate arrays
    size_t n = (size_t)(P.N + 1);
    double *u = calloc(n, sizeof *u);
    double *v = calloc(n, sizeof *v);

    // RK work arrays (reuse to avoid realloc every step)
    double *k1u = calloc(n, sizeof *k1u), *k1v = calloc(n, sizeof *k1v);
    double *k2u = calloc(n, sizeof *k2u), *k2v = calloc(n, sizeof *k2v);
    double *k3u = calloc(n, sizeof *k3u), *k3v = calloc(n, sizeof *k3v);
    double *k4u = calloc(n, sizeof *k4u), *k4v = calloc(n, sizeof *k4v);
    double *utmp = calloc(n, sizeof *utmp), *vtmp = calloc(n, sizeof *vtmp);

    if (!u||!v||!k1u||!k1v||!k2u||!k2v||!k3u||!k3v||!k4u||!k4v||!utmp||!vtmp) {
        perror("calloc"); return 1;
    }

    // initial conditions
    for (int i = 0; i <= P.N; ++i) {
        double x = i * P.dx;
        u[i] = f0(x, P.L);
        v[i] = g0(x);
    }
    apply_bc(u, v, &P);

    // time loop
    for (int nstep = 0; nstep < steps; ++nstep) {
        rk4_step(u, v, dt, k1u,k1v, k2u,k2v, k3u,k3v, k4u,k4v, utmp,vtmp, &P);

        // (Optional) sample output every k steps
        if (nstep % 50 == 0) {
            double t = (nstep+1)*dt;
            printf("# t = %.6f\n", t);
            for (int i = 0; i <= P.N; ++i) {
                double x = i * P.dx;
                printf("%.6f %.6f\n", x, u[i]);
            }
            printf("\n");
        }
    }

    free(u); free(v);
    free(k1u); free(k1v); free(k2u); free(k2v);
    free(k3u); free(k3v); free(k4u); free(k4v);
    free(utmp); free(vtmp);
    return 0;
}
