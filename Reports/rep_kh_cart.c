#include "../paul.h"
#include "../geometry.h"

static double phi_over_pi = 0.0;

void setReportParams(struct domain *theDomain)
{
    phi_over_pi = theDomain->theParList.initPar4;
}

int num_shared_reports()
{
    return 0;
}

int num_distributed_aux_reports()
{
    return 0;
}

int num_distributed_integral_reports()
{
    return 19;
}

void get_shared_reports(double *Q, struct domain *theDomain)
{
    UNUSED(Q);
    UNUSED(theDomain);
}

void get_distributed_aux_reports(double *Q, struct domain *theDomain)
{
    UNUSED(Q);
    UNUSED(theDomain);
}

void get_distributed_integral_reports(const double *xc, const double *prim,
                                      double *Q, struct domain *theDomain)
{
    double xyz[3];
    get_xyz(xc, xyz);

    double phi = phi_over_pi * M_PI;
    double cp = cos(phi);
    double sp = sin(phi);
    
    double x =  cp * xyz[0] + sp * xyz[1];
    double y = -sp * xyz[0] + cp * xyz[1];

    double v[3] = {prim[URR], prim[UPP], prim[UZZ]};
    get_vec_covariant(xc, v, v);
    double vxyz[3];
    get_vec_xyz(xc, v, vxyz);

    double vx =  cp * vxyz[0] + sp * vxyz[1];
    double vy = -sp * vxyz[0] + cp * vxyz[1];
    
    while(y < 0.0)
        y += 1.0;
    while(y > 1.0)
        y -= 1.0;

    double w;

    if(y < 0.5)
        w = exp(-4*M_PI*fabs(y-0.25));
    else
        w = exp(-4*M_PI*fabs(y-0.75));

    Q[0] = w;

    Q[1] = w*vx;
    Q[2] = w*vy;

    Q[3] = w*vx*cos(2*M_PI*x);
    Q[4] = w*vx*sin(2*M_PI*x);
    Q[5] = w*vy*cos(2*M_PI*x);
    Q[6] = w*vy*sin(2*M_PI*x);

    Q[7] = w*vx*cos(4*M_PI*x);
    Q[8] = w*vx*sin(4*M_PI*x);
    Q[9] = w*vy*cos(4*M_PI*x);
    Q[10] = w*vy*sin(4*M_PI*x);

    Q[11] = w*vx*cos(6*M_PI*x);
    Q[12] = w*vx*sin(6*M_PI*x);
    Q[13] = w*vy*cos(6*M_PI*x);
    Q[14] = w*vy*sin(6*M_PI*x);

    Q[15] = w*vx*cos(8*M_PI*x);
    Q[16] = w*vx*sin(8*M_PI*x);
    Q[17] = w*vy*cos(8*M_PI*x);
    Q[18] = w*vy*sin(8*M_PI*x);
}
