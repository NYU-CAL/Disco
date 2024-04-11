
#include "../paul.h"
#include "../boundary.h"

void boundary_trans( struct domain * theDomain , int dim )
{
    if(dim == 1)
    {
        boundary_fixed_rinn(theDomain);
        boundary_zerograd_rout(theDomain, 0);
    }
    else if(dim == 2)
    {
        boundary_reflect_zbot(theDomain);
        boundary_fixed_zbot_rrange(theDomain, -1.0, 1.0/6.0);
        boundary_fixed_ztop(theDomain);
    }
}

