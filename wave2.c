// rk4_cpml_1d.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int    N;        // number of intervals -> N+1 points
    int    PML;      // PML thickness (cells) on EACH side
    double L;        // domain length
    double c;        // wave speed
    double dx;       // grid spacing
} Params;

// ---- Initial conditions (edit) ----
static inline double p0(double x, double L) {
    // e.g. small Gaussian pulse in pressure
    double xc = 0.3 * L, w = 0.03 * L;
    double s = (x - xc)/w;
    return exp(-s*s);
}
static inline double v0(double x, double L) {
    return 0.0;
}

// ---- Convenience ----
static inline int clampi(int a, int lo, int hi){ return a<lo?lo:(a>hi?hi:a); }
static inline int idx(int i){ return i; } // 1-D

// ---- Build PML profiles σ(x), a(x), b(x) (host) ----
static void build_pml(const Params *P, double sigma_max,
                      double *sigma, double *a, double *b)
{
    const int Np1 = P->N + 1;
    for (int i = 0; i < Np1; ++i) {
        sigma[i] = 0.0;
        a[i] = 0.0;
        b[i] = 0.0;
    }
    // cubic ramp on each side
    for (int s = 0; s < P->PML; ++s) {
        double r = (double)(s + 1) / (double)P->PML; // (0,1]
        double ramp = r*r*r;

        // left side cell index
        int iL = s;
        // right side cell index
        int iR = (P->N) - s;

        sigma[iL] = fmax(sigma[iL], sigma_max * ramp);
        sigma[iR] = fmax(sigma[iR], sigma_max * ramp);
    }

    // Simple CPML ODE coefficients:
    // You can refine these; as a starter choose a=b=σ (works reasonably).
    for (int i = 0; i < Np1; ++i) {
        a[i] = sigma[i];  // TODO: try frequency-shifted CPML if you like
        b[i] = sigma[i];
    }
}

// ---- Spatial ops (centred differences with simple edge handling) ----
static inline double ddx_center(const double *f, int i, const Params *P) {
    const int N = P->N;
    const double dx = P->dx;
    int im1 = clampi(i-1, 0, N);
    int ip1 = clampi(i+1, 0, N);
    return (f[idx(ip1)] - f[idx(im1)]) / (2.0*dx);
}

// ---- RHS of ODE system (p, v, psi_p, psi_v) ----
typedef struct {
    double *p, *v;
    double *psi_p, *psi_v;
} State;

static void rhs(const State U, State dU, const Params *P,
                const double *sigma, const double *a, const double *b)
{
    const int Np1 = P->N + 1;
    const double kappa = P->c * P->c;

    // Spatial derivatives (store temps or recompute inline)
    for (int i = 0; i < Np1; ++i) {
        double dpdx = ddx_center(U.p, i, P);
        double dvdx = ddx_center(U.v, i, P);

        // Memory variable ODEs
        dU.psi_p[i] = a[i] * dpdx - b[i] * U.psi_p[i];
        dU.psi_v[i] = a[i] * dvdx - b[i] * U.psi_v[i];

        // Primary fields
        dU.p[i] = -kappa * dvdx - sigma[i]*U.p[i] + kappa * U.psi_v[i];
        dU.v[i] = -dpdx        - sigma[i]*U.v[i] +        U.psi_p[i];
    }

    // Basic boundary conditions (start simple):
    // zero-gradient for v at ends; zero-gradient for p.
    // (You can replace by higher-quality one-sided stencils if you wish.)
    dU.p[0]      = dU.p[1];
    dU.p[P->N]   = dU.p[P->N-1];
    dU.v[0]      = dU.v[1];
    dU.v[P->N]   = dU.v[P->N-1];
    dU.psi_p[0]  = dU.psi_p[1];
    dU.psi_p[P->N] = dU.psi_p[P->N-1];
    dU.psi_v[0]  = dU.psi_v[1];
    dU.psi_v[P->N] = dU.psi_v[P->N-1];
}

// ---- One RK4 step on the whole state ----
static void rk4_step(State U, double dt, const Params *P,
                     const double *sigma, const double *a, const double *b,
                     // work buffers (pre-allocated)
                     State k1, State k2, State k3, State k4,
                     State W1, State W2, State W3)
{
    const int Np1 = P->N + 1;

    // k1 = f(U)
    rhs(U, k1, P, sigma, a, b);

    // W1 = U + dt/2 * k1
    for (int i=0;i<Np1;++i){
        W1.p[i]     = U.p[i]     + 0.5*dt*k1.p[i];
        W1.v[i]     = U.v[i]     + 0.5*dt*k1.v[i];
        W1.psi_p[i] = U.psi_p[i] + 0.5*dt*k1.psi_p[i];
        W1.psi_v[i] = U.psi_v[i] + 0.5*dt*k1.psi_v[i];
    }

    // k2 = f(W1)
    rhs(W1, k2, P, sigma, a, b);

    // W2 = U + dt/2 * k2
    for (int i=0;i<Np1;++i){
        W2.p[i]     = U.p[i]     + 0.5*dt*k2.p[i];
        W2.v[i]     = U.v[i]     + 0.5*dt*k2.v[i];
        W2.psi_p[i] = U.psi_p[i] + 0.5*dt*k2.psi_p[i];
        W2.psi_v[i] = U.psi_v[i] + 0.5*dt*k2.psi_v[i];
    }

    // k3 = f(W2)
    rhs(W2, k3, P, sigma, a, b);

    // W3 = U + dt * k3
    for (int i=0;i<Np1;++i){
        W3.p[i]     = U.p[i]     + dt*k3.p[i];
        W3.v[i]     = U.v[i]     + dt*k3.v[i];
        W3.psi_p[i] = U.psi_p[i] + dt*k3.psi_p[i];
        W3.psi_v[i] = U.psi_v[i] + dt*k3.psi_v[i];
    }

    // k4 = f(W3)
    rhs(W3, k4, P, sigma, a, b);

    // Combine
    for (int i=0;i<Np1;++i){
        U.p[i]     += (dt/6.0)*(k1.p[i] + 2.0*k2.p[i] + 2.0*k3.p[i] + k4.p[i]);
        U.v[i]     += (dt/6.0)*(k1.v[i] + 2.0*k2.v[i] + 2.0*k3.v[i] + k4.v[i]);
        U.psi_p[i] += (dt/6.0)*(k1.psi_p[i] + 2.0*k2.psi_p[i] + 2.0*k3.psi_p[i] + k4.psi_p[i]);
        U.psi_v[i] += (dt/6.0)*(k1.psi_v[i] + 2.0*k2.psi_v[i] + 2.0*k3.psi_v[i] + k4.psi_v[i]);
    }
}

// ---- Main (wire it together) ----
int main(void){
    Params P = { .N = 200, .PML = 20, .L = 1.0, .c = 1.0 };
    P.dx = P.L / P.N;

    // CFL-ish time step (RK4 still needs CFL safety)
    const double safety = 0.7; // you can push towards 0.9; start safe
    const double dt = safety * P.dx / P.c;
    const double Tfinal = 1.0;
    const int steps = (int)ceil(Tfinal / dt);

    const int Np1 = P.N + 1;

    // Allocate fields
    double *p     = calloc(Np1, sizeof *p);
    double *v     = calloc(Np1, sizeof *v);
    double *psi_p = calloc(Np1, sizeof *psi_p);
    double *psi_v = calloc(Np1, sizeof *psi_v);

    // Work arrays for RK4
    double *k1p=calloc(Np1,sizeof* k1p), *k1v=calloc(Np1,sizeof* k1v);
    double *k1psip=calloc(Np1,sizeof* k1psip), *k1psiv=calloc(Np1,sizeof* k1psiv);
    double *k2p=calloc(Np1,sizeof* k2p), *k2v=calloc(Np1,sizeof* k2v);
    double *k2psip=calloc(Np1,sizeof* k2psip), *k2psiv=calloc(Np1,sizeof* k2psiv);
    double *k3p=calloc(Np1,sizeof* k3p), *k3v=calloc(Np1,sizeof* k3v);
    double *k3psip=calloc(Np1,sizeof* k3psip), *k3psiv=calloc(Np1,sizeof* k3psiv);
    double *k4p=calloc(Np1,sizeof* k4p), *k4v=calloc(Np1,sizeof* k4v);
    double *k4psip=calloc(Np1,sizeof* k4psip), *k4psiv=calloc(Np1,sizeof* k4psiv);

    double *w1p=calloc(Np1,sizeof* w1p), *w1v=calloc(Np1,sizeof* w1v);
    double *w1psip=calloc(Np1,sizeof* w1psip), *w1psiv=calloc(Np1,sizeof* w1psiv);
    double *w2p=calloc(Np1,sizeof* w2p), *w2v=calloc(Np1,sizeof* w2v);
    double *w2psip=calloc(Np1,sizeof* w2psip), *w2psiv=calloc(Np1,sizeof* w2psiv);
    double *w3p=calloc(Np1,sizeof* w3p), *w3v=calloc(Np1,sizeof* w3v);
    double *w3psip=calloc(Np1,sizeof* w3psip), *w3psiv=calloc(Np1,sizeof* w3psiv);

    double *sigma = calloc(Np1, sizeof *sigma);
    double *a = calloc(Np1, sizeof *a);
    double *b = calloc(Np1, sizeof *b);

    if (!p||!v||!psi_p||!psi_v||!k1p||!k1v||!k1psip||!k1psiv||
        !k2p||!k2v||!k2psip||!k2psiv||!k3p||!k3v||!k3psip||!k3psiv||
        !k4p||!k4v||!k4psip||!k4psiv||!w1p||!w1v||!w1psip||!w1psiv||
        !w2p||!w2v||!w2psip||!w2psiv||!w3p||!w3v||!w3psip||!w3psiv||
        !sigma||!a||!b) {
        perror("calloc");
        return 1;
    }

    // Tie states
    State U   = { p, v, psi_p, psi_v };
    State k1  = { k1p, k1v, k1psip, k1psiv };
    State k2  = { k2p, k2v, k2psip, k2psiv };
    State k3  = { k3p, k3v, k3psip, k3psiv };
    State k4  = { k4p, k4v, k4psip, k4psiv };
    State W1  = { w1p, w1v, w1psip, w1psiv };
    State W2  = { w2p, w2v, w2psip, w2psiv };
    State W3  = { w3p, w3v, w3psip, w3psiv };

    // Build PML profiles
    const double sigma_max = 60.0 * (P.c / P.L); // TODO: tune; rule-of-thumb scale
    build_pml(&P, sigma_max, sigma, a, b);

    // Initial conditions
    for (int i = 0; i < Np1; ++i) {
        double x = i * P.dx;
        p[i] = p0(x, P.L);
        v[i] = v0(x, P.L);
        psi_p[i] = 0.0;
        psi_v[i] = 0.0;
    }

    // Time loop
    for (int n = 0; n < steps; ++n) {
        // (Optional) add source on p or v here (e.g. Ricker at centre)
        // TODO: source injection before/after rk4_step (be consistent)

        rk4_step(U, dt, &P, sigma, a, b, k1,k2,k3,k4, W1,W2,W3);

        // Simple sampling
        if (n % 50 == 0) {
            double t = (n+1)*dt;
            printf("# t = %.6f\n", t);
            for (int i=0;i<Np1;++i) {
                double x = i * P.dx;
                printf("%.6f %.6f\n", x, p[i]);
            }
            printf("\n");
        }
    }

    // cleanup
    free(p); free(v); free(psi_p); free(psi_v);
    free(k1p); free(k1v); free(k1psip); free(k1psiv);
    free(k2p); free(k2v); free(k2psip); free(k2psiv);
    free(k3p); free(k3v); free(k3psip); free(k3psiv);
    free(k4p); free(k4v); free(k4psip); free(k4psiv);
    free(w1p); free(w1v); free(w1psip); free(w1psiv);
    free(w2p); free(w2v); free(w2psip); free(w2psiv);
    free(w3p); free(w3v); free(w3psip); free(w3psiv);
    free(sigma); free(a); free(b);
    return 0;
}
