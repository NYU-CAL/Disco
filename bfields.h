#ifndef DISCO_BFIELDS_H
#define DISCO_BFIELDS_H

#include "paul.h"
 
void setBfieldsParams(struct domain *theDomain);
void set_B_fields(struct domain *theDomain);
void B_faces_to_cells(struct domain *theDomain, int type);
void update_B_fluxes(struct domain *theDomain, double dt);
void add_E_phi(double *phiL, double *phiR, double *phiD, double *phiU,
               double Edldt);

void avg_Efields(struct domain *theDomain);
void avg_Efields_Extra_Duffell16( struct domain * theDomain );
void avg_Efields_Extra_GardinerStone05( struct domain * theDomain );

void subtract_advective_B_fluxes(struct domain *theDomain);
void check_flipped(struct domain *theDomain, int dim);
void get_phi_pointer(struct face *f, double **P, double **RK_P, int dim);
void flip_fluxes(struct domain *theDomain, int dim);

int phi_switch(double dphi, double Pmax, int mode);
int get_which4(double phi, double phiR, double phiU, double phiUR,
               int *LR_alt, int *UD_alt, int mode, double Pmax);

void make_edge_adjust(struct domain *theDomain, double dt);


#endif
