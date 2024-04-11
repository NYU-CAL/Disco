#include "../paul.h"
#include "../geometry.h"

static double gamma = 0.0;
static double shockMach = 10.0;
static double shockAngleDeg = 30.0;
static double rho0 = 1.0;
static double cs0 = 1.0;
static double l0 = 1.0/6.0;
static double *t_ptr = NULL;


void setICparams(struct domain *theDomain)
{
    gamma = theDomain->theParList.Adiabatic_Index;
    t_ptr = &(theDomain->t);
}

void initial(double *prim, const double *xc)
{
    double xyz[3];
    get_xyz(xc, xyz);

    double t = *t_ptr;

    double cp = cos(-shockAngleDeg * M_PI/180.0);
    double sp = sin(-shockAngleDeg * M_PI/180.0);

    double x =  cp * xyz[0] + sp * xyz[2];

    double x0 = l0 * cp;

    double shockV = shockMach * cs0;

    double rho, vx, P, f;

    if(x < x0 + shockV * t)
    {
        double M2 = shockMach * shockMach;

        rho = (gamma+1)*M2 / ((gamma-1)*M2 + 2) * rho0;
        vx = 2*(M2-1.0) / ((gamma+1) * shockMach) * cs0;
        P = (2*gamma*M2 - (gamma-1)) / (gamma+1) * cs0*cs0*rho0/gamma;
        f = 1;
    }
    else
    {
        rho = rho0;
        vx = 0.0;
        P = rho0*cs0*cs0/gamma;
        f = 0;
    }

    double vxyz[3] = {cp*vx, 0.0, sp*vx};
    double v[3];

    get_vec_from_xyz(xc, vxyz, v);
    get_vec_contravariant(xc, v, v);
    
    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = v[0];
    prim[UPP] = v[1];
    prim[UZZ] = v[2];
    
    if(NUM_N > 0)
    {
        prim[NUM_C] = f;
    }
}
