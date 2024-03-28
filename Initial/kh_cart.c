#include "../paul.h"
#include "../geometry.h"

static int mode = 0;
static double rhoIn = 0.0;
static double dv = 0.0;
static double P0 = 0.0;
static double phi_over_pi = 0.0;
static double n = 0.0;
static double amp = 0.0;
static double width = 0.0;

void setICparams(struct domain *theDomain)
{
    mode = theDomain->theParList.initPar0;

    rhoIn = theDomain->theParList.initPar1;
    dv = theDomain->theParList.initPar2;
    P0 = theDomain->theParList.initPar3;
    phi_over_pi = theDomain->theParList.initPar4;
    n = theDomain->theParList.initPar5;
    amp = theDomain->theParList.initPar6;
    width = theDomain->theParList.initPar7;
}

void initial(double *prim, const double *xc)
{
    double xyz[3];
    get_xyz(xc, xyz);

    double phi = phi_over_pi * M_PI;
    double cp = cos(phi);
    double sp = sin(phi);
    
    double x =  cp * xyz[0] + sp * xyz[1];
    double y = -sp * xyz[0] + cp * xyz[1];
    
    double rhoOut = 1.0;

    double vOut = 0.5*dv;
    double vIn = -0.5*dv;

    double wIn, wOut;

    if(mode == 1)
    {
        wIn = 0.0;
        wOut = 0.0;
    }
    else
    {
        while(y < 0.0)
            y += 1.0;
        while(y > 1.0)
            y -= 1.0;

        if(y < 0.25)
        {
            double w = exp((y-0.25)/width);
            wOut = 1.0 - 0.5*w;
            wIn = 0.5*w;
        }
        else if(y < 0.5)
        {
            double w = exp(-(y-0.25)/width);
            wOut = 0.5*w;
            wIn = 1-0.5*w;
        }
        else if(y < 0.75)
        {
            double w = exp((y-0.75)/width);
            wOut = 0.5*w;
            wIn = 1-0.5*w;
        }
        else
        {
            double w = exp(-(y-0.75)/width);
            wOut = 1-0.5*w;
            wIn = 0.5*w;
        }
    }

    double rho = wOut * rhoOut + wIn * rhoIn;
    double vx = wOut * vOut + wIn * vIn;
    double vy = amp * sin(2*M_PI*n*x);
    double f = wIn < 0.5 ? 0.0 : 1.0;

    double vxyz[3] = {cp*vx-sp*vy, sp*vx+cp*vy, 0.0};
    double v[3];

    get_vec_from_xyz(xc, vxyz, v);

    
    prim[RHO] = rho;
    prim[PPP] = P0;
    prim[URR] = v[0];
    prim[UPP] = v[1];
    prim[UZZ] = v[2];
    
    if(NUM_N > 0)
    {
        prim[NUM_C] = f;
    }
}
