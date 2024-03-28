
#include "paul.h"
#include "geometry.h"
#include "hydro.h"
#include "bfields.h"

void initial( double * , double * ); 
void setup_faces( struct domain * , int );
int get_num_rzFaces( int , int , int );
void calc_prim(struct domain *);

static int CT_Solver = 0;
 
void setBfieldsParams(struct domain *theDomain)
{
    CT_Solver = theDomain->theParList.CT_Solver;
}

void set_B_fields( struct domain * theDomain ){

   int i,j,k;
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

#if NUM_FACES > 0
    for( k=0 ; k<Nz ; ++k ){
      double z = get_centroid(z_kph[k], z_kph[k-1], 2);
      for( j=0 ; j<Nr ; ++j ){
         double r = get_centroid(r_jph[j], r_jph[j-1], 1);
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = { r_jph[j]   , phip , z_kph[k]   };
            double xm[3] = { r_jph[j-1] , phim , z_kph[k-1] };
            double x[3] = { r , phip , z};
            double prim[NUM_Q];
            initial( prim , x ); 
            double dA = get_dA( xp , xm , 0 );
            double Phi = 0.0;
            if( NUM_Q > BPP ) Phi = prim[BPP]*dA;
                c->Phi[0] = Phi;
         }    
      }    
    }
#endif

   int NRZ1 = theDomain->N_ftracks_r;
   setup_faces( theDomain , 1 ); 

#if NUM_FACES >= 3
   int n;
   for( n=0 ; n<theDomain->fIndex_r[NRZ1] ; ++n ){
      struct face * f = theDomain->theFaces_1 + n;
      double prim[NUM_Q];
      initial( prim , f->cm );
      double Phi = 0.0;
      if( NUM_Q > BRR ) Phi = prim[BRR]*f->dA;
      if( f->LRtype == 0 ){
         f->L->Phi[2] = Phi;
      }else{
         f->R->Phi[1] = Phi;
      }
   }
#endif

#if NUM_FACES == 5
   if(theDomain->Nz > 1 ){
      int NRZ2 = theDomain->N_ftracks_z;
      setup_faces( theDomain , 2 );
      for( n=0 ; n<theDomain->fIndex_z[NRZ2] ; ++n ){
         struct face * f = theDomain->theFaces_2 + n; 
         double prim[NUM_Q];
         initial( prim , f->cm );
         double Phi = 0.0; 
         if( NUM_Q > BZZ ) Phi = prim[BZZ]*f->dA;

         if( f->LRtype == 0 ){ 
            f->L->Phi[4] = Phi; 
         }else{
            f->R->Phi[3] = Phi; 
         }    
      }
   }
#endif

   B_faces_to_cells( theDomain , 0 );
   B_faces_to_cells( theDomain , 1 );
 
   free( theDomain->theFaces_1 );
   if( theDomain->theFaces_2 ) free( theDomain->theFaces_2 );

}

void B_faces_to_cells( struct domain * theDomain , int type ){

   if( NUM_Q > BZZ ){
      struct cell ** theCells = theDomain->theCells;
      struct face * theFaces_1 = theDomain->theFaces_1;
      struct face * theFaces_2 = theDomain->theFaces_2;
   
      int Nr = theDomain->Nr;
      int Nz = theDomain->Nz;
      //int NgRa = theDomain->NgRa;
      //int NgRb = theDomain->NgRb;
      //int NgZa = theDomain->NgZa;
      //int NgZb = theDomain->NgZb;
      int * Np = theDomain->Np;
      double * r_jph = theDomain->r_jph;
      double * z_kph = theDomain->z_kph;

      int Nf = theDomain->fIndex_r[theDomain->N_ftracks_r];  
 
      int i,j,k;

      for( j=0 ; j<Nr ; ++j ){
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               c->tempDoub = 0.0;
               if( type==0 ){
                  c->prim[BRR] = 0.0;
                  c->prim[BPP] = 0.0;
                  if( NUM_FACES==5 ) c->prim[BZZ] = 0.0;
               }else{
                  c->cons[BRR] = 0.0;
                  c->cons[BPP] = 0.0;
                  if( NUM_FACES==5 ) c->cons[BZZ] = 0.0;
               }
            }
         }
      }   

      int n;
      for( n=0 ; n<Nf ; ++n ){
         struct face * f = theFaces_1 + n;
         struct cell * cL = f->L;
         struct cell * cR = f->R;
         double Phi;
         if( f->LRtype==0 ){
            Phi = cL->Phi[2];
         }else{
            Phi = cR->Phi[1];
         }
         cL->tempDoub += f->dA;
         cR->tempDoub += f->dA;

         if( type==0 ){
            cL->prim[BRR] += Phi;
            cR->prim[BRR] += Phi;
         }else{
            cL->cons[BRR] += Phi;
            cR->cons[BRR] += Phi;
         }
      }

      for( j=0 ; j<Nr ; ++j ){
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               int im = i-1;
               if( im==-1 ) im = Np[jk]-1;

               struct cell * c = &(theCells[jk][i]);
               struct cell * cm = &(theCells[jk][im]);

               double xp[3] = { r_jph[j]   , c->piph  , z_kph[k]   };
               double xm[3] = { r_jph[j-1] , cm->piph , z_kph[k-1] };
               double x[3];
               get_centroid_arr(xp, xm, x);
               double dA = get_dA( xp , xm , 0 );
               double dV = get_dV( xp , xm );

               if( type==0 ){
                  c->prim[BRR] /= c->tempDoub;
                  c->prim[BPP] = .5*(c->Phi[0]+cm->Phi[0])/dA;
               }else{
                  double rfac = bfield_scale_factor(x[0], 0);
                  double hr = get_scale_factor(x, 1);
                  double hp = get_scale_factor(x, 0);
                  c->cons[BRR] *= dV*rfac/(c->tempDoub * hr) ;
                  c->cons[BPP] = .5*(c->Phi[0]+cm->Phi[0])*dV/(dA * hp);
               }
            }
         }
      }
      /*
      if(NgRa == 0)
      {
          j=0;
          for(k=0; k<Nz; k++)
          {
              int jk = j + Nr*k;
              for(i=0; i<Np[jk]; i++)
              {
                  if(type == 0)
                      theCells[jk][i].prim[BRR] *= 0.5;
                  else
                      theCells[jk][i].cons[BRR] *= 0.5;
              }
          }
      }

      if(NgRb == 0)
      {
          j=Nr-1;
          for(k=0; k<Nz; k++)
          {
              int jk = j + Nr*k;
              for(i=0; i<Np[jk]; i++)
              {
                  if(type == 0)
                      theCells[jk][i].prim[BRR] *= 0.5;
                  else
                      theCells[jk][i].cons[BRR] *= 0.5;
              }
          }
      }
      */
      

      if( NUM_FACES == 5 && Nz>1 ){
         for( j=0 ; j<Nr ; ++j ){
            for( k=0 ; k<Nz ; ++k ){
               int jk = j+Nr*k;
               for( i=0 ; i<Np[jk] ; ++i ){
                  theCells[jk][i].tempDoub = 0.0;
               }
            }
         }
         int Nfz = theDomain->fIndex_z[theDomain->N_ftracks_z];  
         for( n=0 ; n<Nfz ; ++n ){
            struct face * f = theFaces_2 + n;
            struct cell * cL = f->L;
            struct cell * cR = f->R;
            double Phi;
            if( f->LRtype==0 ){
               Phi = cL->Phi[4];
            }else{
               Phi = cR->Phi[3];
            }
            cL->tempDoub += f->dA;
            cR->tempDoub += f->dA;

            if( type==0 ){
               cL->prim[BZZ] += Phi;
               cR->prim[BZZ] += Phi;
            }else{
               cL->cons[BZZ] += Phi;
               cR->cons[BZZ] += Phi;
            }
         }

         for( j=0 ; j<Nr ; ++j ){
            for( k=0 ; k<Nz ; ++k ){
               int jk = j+Nr*k;
               for( i=0 ; i<Np[jk] ; ++i ){
                  int im = i-1; 
                  if( im==-1 ) im = Np[jk]-1;

                  struct cell * c = &(theCells[jk][i]);
                  struct cell * cm = &(theCells[jk][im]);

                  double xp[3] = { r_jph[j]   , c->piph  , z_kph[k]   };   
                  double xm[3] = { r_jph[j-1] , cm->piph , z_kph[k-1] };
                  double x[3];
                  get_centroid_arr(xp, xm, x);
                  double dV = get_dV( xp , xm );

                  if( type==0 ){
                     c->prim[BZZ] /= c->tempDoub;
                  }else{
                     double zfac = bfield_scale_factor(x[2], 2);
                     double hz = get_scale_factor(x, 2);
                     c->cons[BZZ] *= dV*zfac/(c->tempDoub * hz);
                  }    
               }    
            }    
         }  
        /*
        if(NgZa == 0)
        {
            k=0;
            for(j=0; j<Nr; j++)
            {
                int jk = j + Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    if(type == 0)
                        theCells[jk][i].prim[BZZ] *= 0.5;
                    else
                        theCells[jk][i].cons[BZZ] *= 0.5;
                }
            }
        }

        if(NgZb == 0)
        {
            k=Nz-1;
            for(j=0; j<Nr; j++)
            {
                int jk = j + Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    if(type == 0)
                        theCells[jk][i].prim[BZZ] *= 0.5;
                    else
                        theCells[jk][i].cons[BZZ] *= 0.5;
                }
            }
        }
        */
      }
   }
}

void update_B_fluxes( struct domain * theDomain , double dt ){

   struct face * theFaces = theDomain->theFaces_1;
   int * Nf = theDomain->fIndex_r;
   int Njk = theDomain->N_ftracks_r;
   int n;
   int jk;
   int n0 = 0;
   if(NUM_FACES >= 3 && NUM_EDGES >= 2)
   {
       for( jk=0 ; jk<Njk ; ++jk ){
          for( n=0 ; n<Nf[jk]-n0 ; ++n ){
             struct face * f  = theFaces + n0 + n;
             int np = n+1;
             if( np==Nf[jk]-n0 ) np=0; 
             struct face * fp = theFaces + n0 + np;

             double E;
             double dl = f->dl;
             if( f->LRtype == 0 ){
                E = f->L->E[1];
                f->L->Phi[2] -= E*dl*dt;
                f->L->Phi[0] += E*dl*dt;
             }else{
                E = f->R->E[0];
                f->R->Phi[1] -= E*dl*dt;
                f->R->Phi[0] -= E*dl*dt;
             }
             if( fp->LRtype == 0 ){
                fp->L->Phi[2] += E*dl*dt;
             }else{
                fp->R->Phi[1] += E*dl*dt;
             }
          }
          n0 = Nf[jk];
       }
   }

   if( NUM_FACES == 5 && NUM_EDGES == 8 ){
      theFaces = theDomain->theFaces_2;
      Nf = theDomain->fIndex_z;
      Njk = theDomain->N_ftracks_z;
      n0 = 0; 
      for( jk=0 ; jk<Njk ; ++jk ){
         for( n=0 ; n<Nf[jk]-n0 ; ++n ){
            struct face * f  = theFaces + n0 + n; 
            int np = n+1; 
            if( np==Nf[jk]-n0 ) np=0; 
            struct face * fp = theFaces + n0 + np;

            double E;
            double dl = f->dl;
            if( f->LRtype == 0 ){ 
               E = f->L->E[5];
               f->L->Phi[4] += E*dl*dt;
               f->L->Phi[0] -= E*dl*dt;
            }else{
               E = f->R->E[4];
               f->R->Phi[3] += E*dl*dt;
               f->R->Phi[0] += E*dl*dt;
            }    
            if( fp->LRtype == 0 ){ 
               fp->L->Phi[4] -= E*dl*dt;
            }else{
               fp->R->Phi[3] -= E*dl*dt;
            }    
         }    
         n0 = Nf[jk];
      }

      if( NUM_AZ_EDGES == 4 && theDomain->Nz>1 ) 
          make_edge_adjust( theDomain , dt );
   } 
}

void add_E_phi(double *phiL, double *phiR, double *phiD, double *phiU,
               double Edldt )
{
    *phiL -= Edldt;
    *phiR += Edldt;
    *phiU -= Edldt;
    *phiD += Edldt;
}

void avg_Efields( struct domain * theDomain ){

    if(NUM_EDGES < 4)
        return;

    int i,j,k;
    struct cell ** theCells = theDomain->theCells;
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;

    for( j=0 ; j<Nr ; ++j )
    {
        for( k=0 ; k<Nz ; ++k )
        {
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i )
            {
                struct cell * c  = theCells[jk]+i;
                int ip = (i+1)%Np[jk];
                struct cell * cp = theCells[jk]+ip;

                double El_avg = .5*( c->E[0] + cp->E[2] );
                double Er_avg = .5*( c->E[1] + cp->E[3] );

                 c->E[0] = El_avg;
                 c->E[1] = Er_avg;
                cp->E[2] = El_avg;
                cp->E[3] = Er_avg;

                double Bl_avg = .5*( c->B[0] + cp->B[2] );
                double Br_avg = .5*( c->B[1] + cp->B[3] );

                 c->B[0] = Bl_avg;
                 c->B[1] = Br_avg;
                cp->B[2] = Bl_avg;
                cp->B[3] = Br_avg;

                if( NUM_EDGES == 8 )
                {
                    El_avg = .5*( c->E[4] + cp->E[6] );
                    Er_avg = .5*( c->E[5] + cp->E[7] );

                     c->E[4] = El_avg;
                     c->E[5] = Er_avg;
                    cp->E[6] = El_avg;
                    cp->E[7] = Er_avg;

                    Bl_avg = .5*( c->B[4] + cp->B[6] );
                    Br_avg = .5*( c->B[5] + cp->B[7] );

                     c->B[4] = Bl_avg;
                     c->B[5] = Br_avg;
                    cp->B[6] = Bl_avg;
                    cp->B[7] = Br_avg;
                }

            }
        }
    }

    if(CT_Solver == 0)
        avg_Efields_Extra_Duffell16(theDomain);
    else if(CT_Solver == 1)
        avg_Efields_Extra_GardinerStone05(theDomain);
   
   //E ALONG THE POLE...
   if(theDomain->NgRa == 0 && NUM_EDGES >= 4)
   {
      j=0;
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         double E = 0.0;
         double B = 0.0;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = theCells[jk]+i;
            E += c->E[1];
         }
         E /= Np[jk];
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = theCells[jk]+i;
            c->E[0] = E;
            c->B[0] = B;
         }
      }
   }
   if(theDomain->NgRb == 0 && NUM_EDGES >= 4)
   {
      j=Nr-1;
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         double E = 0.0;
         double B = 0.0;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = theCells[jk]+i;
            E += c->E[0];
         }
         E /= Np[jk];
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = theCells[jk]+i;
            c->E[1] = E;
            c->B[1] = B;
         }
      }
   }
   if(NUM_EDGES == 8 && theDomain->NgZa == 0)
   {
      k=0;
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
         double E = 0.0;
         double B = 0.0;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = theCells[jk]+i;
            E += c->E[5];
         }
         E /= Np[jk];
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = theCells[jk]+i;
            c->E[5] = E;
            c->B[5] = B;
         }
      }
   }
   if(NUM_EDGES == 8 && theDomain->NgZb == 0)
   {
      k=Nz-1;
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
         double E = 0.0;
         double B = 0.0;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = theCells[jk]+i;
            E += c->E[4];
         }
         E /= Np[jk];
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = theCells[jk]+i;
            c->E[4] = E;
            c->B[4] = B;
         }
      }
   }

   if(NUM_EDGES >= 4)
   {
       for( j=0 ; j<Nr ; ++j ){
          for( k=0 ; k<Nz ; ++k ){
             int jk = j+Nr*k;
             for( i=0 ; i<Np[jk] ; ++i ){
                struct cell * c  = theCells[jk]+i;
                int ip = (i+1)%Np[jk];
                struct cell * cp = theCells[jk]+ip;

                cp->E[2] = c->E[0];   
                cp->E[3] = c->E[1];   
                cp->B[2] = c->B[0];   
                cp->B[3] = c->B[1];  

                if( NUM_EDGES == 8 ){
     
                   cp->E[6] = c->E[4];   
                   cp->E[7] = c->E[5];   
                   cp->B[6] = c->B[4];   
                   cp->B[7] = c->B[5];  

                }
             }
          }
       }
   }
}


void avg_Efields_Extra_Duffell16( struct domain * theDomain ) 
{
    // Perform extra averaging for Er and Ez according to Duffell 2016
    int Nf = theDomain->fIndex_r[theDomain->N_ftracks_r];
    struct face * theFaces = theDomain->theFaces_1;

    int n;

    for( n=0 ; n<Nf ; ++n )
    {
        //This loop computes E & B interp and stores it in the face
        //whose forward edge is the interpolated edge
        struct face * f = theFaces+n;
        struct cell * c1;
        struct cell * c2;
        if( f->LRtype == 0 ){
            c1 = f->L;
            c2 = f->R;
        }else{
            c1 = f->R;
            c2 = f->L;
        }
        // c1's front face is trailing c2's ==> c1's front face sets the
        // forward boundary of this face.

        double p1 = c1->piph; //equal to forward phi of this face
        double p2 = c2->piph; //past the forward phi of this faces
        double dp1 = get_dp(p2,p1); //positive, will be in [0, pi]. 
                                    //distance from c2's front face to
                                    //forward phi of this face.
        double dp2 = c2->dphi - dp1; //positive.  distance from c2's back
                                     //face to forward phi of this face.
        if( f->LRtype == 0 )
        {
            double Eavg = ( dp2*c2->E[0] + dp1*c2->E[2] )/(dp1+dp2);
            double Bavg = ( dp2*c2->B[0] + dp1*c2->B[2] )/(dp1+dp2);
            f->E = .5*(f->L->E[1] + Eavg);
            f->B = .5*(f->L->B[1] + Bavg);
        }
        else
        {
            double Eavg = ( dp2*c2->E[1] + dp1*c2->E[3] )/(dp1+dp2);
            double Bavg = ( dp2*c2->B[1] + dp1*c2->B[3] )/(dp1+dp2);
            f->E = .5*(f->R->E[0] + Eavg);
            f->B = .5*(f->R->B[0] + Bavg);
        }
    }

    for( n=0 ; n<Nf ; ++n )
    {
        struct face * f = theFaces+n;
        if( f->LRtype==0 )
        {
            f->L->E[1] = f->E;
            f->L->B[1] = f->B;
        }
        else
        {
            f->R->E[0] = f->E;
            f->R->B[0] = f->B;
        }
        f->E = 0.0;
        f->B = 0.0;
    }

    if( NUM_EDGES == 8 )
    {
        //REPEAT THE ABOVE FOR VERTICALLY-ORIENTED FACES & RADIAL EDGES
        Nf = theDomain->fIndex_z[theDomain->N_ftracks_z];
        theFaces = theDomain->theFaces_2;
        int n;
        for( n=0 ; n<Nf ; ++n )
        {
            struct face * f = theFaces+n;
            struct cell * c1;
            struct cell * c2;
            if( f->LRtype == 0 )
            { 
                c1 = f->L;
                c2 = f->R;
            }
            else
            {
                c1 = f->R;
                c2 = f->L;
            }    
            double p1 = c1->piph;
            double p2 = c2->piph;
            double dp1 = get_dp(p2,p1);
            double dp2 = c2->dphi - dp1; 
            if( f->LRtype == 0 )
            { 
                double Eavg = ( dp2*c2->E[4] + dp1*c2->E[6] )/(dp1+dp2);
                double Bavg = ( dp2*c2->B[4] + dp1*c2->B[6] )/(dp1+dp2);
                f->E = .5*(f->L->E[5] + Eavg);
                f->B = .5*(f->L->B[5] + Bavg);
            }
            else
            {
                double Eavg = ( dp2*c2->E[5] + dp1*c2->E[7] )/(dp1+dp2);
                double Bavg = ( dp2*c2->B[5] + dp1*c2->B[7] )/(dp1+dp2);
                f->E = .5*(f->R->E[4] + Eavg);
                f->B = .5*(f->R->B[4] + Bavg);
            }    
        }

        for( n=0 ; n<Nf ; ++n )
        {
            struct face * f = theFaces+n;
            if( f->LRtype==0 )
            {
                f->L->E[5] = f->E;
                f->L->B[5] = f->B;
            }
            else
            {
                f->R->E[4] = f->E;
                f->R->B[4] = f->B;
            }    
            f->E = 0.0; 
            f->B = 0.0; 
        }
    }
}
void avg_Efields_Extra_GardinerStone05( struct domain * theDomain )
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;

    int *fI_r  = theDomain->fIndex_r;
    struct face *theFaces_r = theDomain->theFaces_1;
    int Nfr = Nr - 1;

    int kmin = NgZa;
    int kmax = Nz - NgZb;
    int jmin = NgRa == 0 ? 0 : NgRa - 1;
    int jmax = NgRb == 0 ? Nr-1 : Nr - NgRb;

    int k;
    for(k=kmin; k<kmax; k++)
    {
        double zm = theDomain->z_kph[k-1];
        double zp = theDomain->z_kph[k];
        double z = get_centroid(zp, zm, 2);

        int j;
        for(j=jmin; j<jmax; j++)
        {
            // j is inner annulus, j+1 the outer.
            double rmm = theDomain->r_jph[j-1];  // innermost face r
            double rf = theDomain->r_jph[j];     // r at face between annuli
            double rpp = theDomain->r_jph[j+1];  // outermote face r

            double rm = get_centroid(rf, rmm, 1); // r of inner annulus
            double rp = get_centroid(rpp, rf, 1); // r of outer annulus

            int JK = j + Nfr * k;

            int f;
            for(f=fI_r[JK]; f<fI_r[JK+1]; f++)
            {
                //Looping over faces.  fp is the next face.
                int fp = f < fI_r[JK+1]-1 ? f+1 : fI_r[JK];

                // cell C is shared by the faces. L & R on the other side.
                struct cell *cC = NULL;
                struct cell *cL = NULL;
                struct cell *cR = NULL;

                double rC = 0;
                double rLR = 0;

                int idx_EB_L = 0;
                int idx_EB_R = 0;

                if(theFaces_r[f].L == theFaces_r[fp].L)
                { 
                    cC = theFaces_r[f].L;
                    cL = theFaces_r[f].R;
                    cR = theFaces_r[fp].R;
                    rC = rm;
                    rLR = rp;
                    idx_EB_L = 0;
                    idx_EB_R = 2;
                }
                else if(theFaces_r[f].R == theFaces_r[fp].R)
                {
                    cC = theFaces_r[f].R;
                    cL = theFaces_r[f].L;
                    cR = theFaces_r[fp].L;
                    rC = rp;
                    rLR = rm;
                    idx_EB_L = 1;
                    idx_EB_R = 3;
                }
                else
                {
                    fprintf(stderr, "Faces don't share a cell!\n");
                }

                double xC[] = {rC,  cC->piph - 0.5*cC->dphi, z};
                double xL[] = {rLR, cL->piph - 0.5*cL->dphi, z};
                double xR[] = {rLR, cR->piph - 0.5*cR->dphi, z};

                double EL[3], EC[3], ER[3];

                prim_to_E(cL->prim, EL, xL);
                prim_to_E(cC->prim, EC, xC);
                prim_to_E(cR->prim, ER, xR);

                double x[] = {rf, cL->piph, z};

                double BrC = cC->prim[BRR];
                double BrL = cL->prim[BRR];
                double BrR = cR->prim[BRR];

                // Re-orient to put all phi's on same branch.
                xL[1] = x[1] + get_signed_dp(xL[1], x[1]);
                xC[1] = x[1] + get_signed_dp(xC[1], x[1]);
                xR[1] = x[1] + get_signed_dp(xR[1], x[1]);

                // Need to build a 2nd order estimate of Ez at this edge.
                // First attempt: treat r and phi as orthogonal 2d coords,
                //    form planar approx of Ez.

                // Using "dx" and "dA" here but these are not really areas.
                // But fine to 2nd order! (I think :p)
                double dxL = xL[0] - xC[0];
                double dyL = xL[1] - xC[1];
                double dxR = xR[0] - xC[0];
                double dyR = xR[1] - xC[1];
                double dx = x[0] - xC[0];
                double dy = x[1] - xC[1];

                double dA = dxL * dyR - dyL * dxR;

                double Ez_cells_avg = EC[2] 
                                + ((dx*dyR - dy*dxR) * (EL[2]-EC[2]) 
                                +  (dxL*dy - dyL*dx) * (ER[2]-EC[2])) / dA;
                double Br_cells_avg = BrC
                                + ((dx*dyR - dy*dxR) * (BrL-BrC) 
                                +  (dxL*dy - dyL*dx) * (BrR-BrC)) / dA;

                cL->E[idx_EB_L] = 2 * cL->E[idx_EB_L] - Ez_cells_avg;
                cR->E[idx_EB_R] = 2 * cR->E[idx_EB_R] - Ez_cells_avg;
                cL->B[idx_EB_L] = 2 * cL->B[idx_EB_L] - Br_cells_avg;
                cR->B[idx_EB_R] = 2 * cR->B[idx_EB_R] - Br_cells_avg;
            }
        }
    }

#if NUM_EDGES == 8

    int *fI_z  = theDomain->fIndex_z;
    struct face *theFaces_z = theDomain->theFaces_2;
    int Nfz = Nz - 1;

    kmin = NgZa == 0 ? 0 : NgZa - 1;
    kmax = NgZb == 0 ? Nz-1 : Nz - NgZb;
    jmin = NgRa;
    jmax = Nr - NgRb;

    for(k=kmin; k<kmax; k++)
    {
        // k is lower annulus, k+1 the upper.
        double zmm = theDomain->z_kph[k-1];
        double zpp = theDomain->z_kph[k+1];
        double zf = theDomain->z_kph[k];

        double zm = get_centroid(zf, zmm, 2); // z of upper annulus
        double zp = get_centroid(zpp, zf, 2); // z of lower annulus

        int j;
        for(j=jmin; j<jmax; j++)
        {
            double rm = theDomain->r_jph[j-1];
            double rp = theDomain->r_jph[j];
            double r = get_centroid(rp, rm, 1);

            int JK = j + Nfz * k;

            int f;
            for(f=fI_z[JK]; f<fI_z[JK+1]; f++)
            {
                //Looping over faces.  fp is the next face.
                int fp = f < fI_z[JK+1]-1 ? f+1 : fI_z[JK];

                // cell C is shared by the faces. L & R on the other side.
                struct cell *cC = NULL;
                struct cell *cL = NULL;
                struct cell *cR = NULL;

                double zC = 0;
                double zLR = 0;

                int idx_EB_L = 0;
                int idx_EB_R = 0;

                if(theFaces_z[f].L == theFaces_z[fp].L)
                { 
                    cC = theFaces_z[f].L;
                    cL = theFaces_z[f].R;
                    cR = theFaces_z[fp].R;
                    zC = zm;
                    zLR = zp;
                    idx_EB_L = 4;
                    idx_EB_R = 6;
                }
                else if(theFaces_z[f].R == theFaces_z[fp].R)
                {
                    cC = theFaces_z[f].R;
                    cL = theFaces_z[f].L;
                    cR = theFaces_z[fp].L;
                    zC = zp;
                    zLR = zm;
                    idx_EB_L = 5;
                    idx_EB_R = 7;
                }
                else
                {
                    fprintf(stderr, "Faces don't share a cell!\n");
                }

                double xC[] = {r, cC->piph - 0.5*cC->dphi, zC};
                double xL[] = {r, cL->piph - 0.5*cL->dphi, zLR};
                double xR[] = {r, cR->piph - 0.5*cR->dphi, zLR};

                double EL[3], EC[3], ER[3];

                prim_to_E(cL->prim, EL, xL);
                prim_to_E(cC->prim, EC, xC);
                prim_to_E(cR->prim, ER, xR);

                double x[] = {r, cL->piph, zf};

                double BzC = cC->prim[BZZ];
                double BzL = cL->prim[BZZ];
                double BzR = cR->prim[BZZ];

                // Re-orient to put all phi's on same branch.
                xL[1] = x[1] + get_signed_dp(xL[1], x[1]);
                xC[1] = x[1] + get_signed_dp(xC[1], x[1]);
                xR[1] = x[1] + get_signed_dp(xR[1], x[1]);

                // Need to build a 2nd order estimate of Ez at this edge.
                // First attempt: treat r and phi as orthogonal 2d coords,
                //    form planar approx of Ez.

                // Using "dx" and "dA" here but these are not really areas.
                // But fine to 2nd order! (I think :p)
                double dxL = xL[2] - xC[2];
                double dyL = xL[1] - xC[1];
                double dxR = xR[2] - xC[2];
                double dyR = xR[1] - xC[1];
                double dx = x[2] - xC[2];
                double dy = x[1] - xC[1];

                double dA = dxL * dyR - dyL * dxR;

                double Er_cells_avg = EC[0] 
                                + ((dx*dyR - dy*dxR) * (EL[0]-EC[0]) 
                                +  (dxL*dy - dyL*dx) * (ER[0]-EC[0])) / dA;
                double Bz_cells_avg = BzC
                                + ((dx*dyR - dy*dxR) * (BzL-BzC) 
                                +  (dxL*dy - dyL*dx) * (BzR-BzC)) / dA;

                cL->E[idx_EB_L] = 2 * cL->E[idx_EB_L] - Er_cells_avg;
                cR->E[idx_EB_R] = 2 * cR->E[idx_EB_R] - Er_cells_avg;
                cL->B[idx_EB_L] = 2 * cL->B[idx_EB_L] - Bz_cells_avg;
                cR->B[idx_EB_R] = 2 * cR->B[idx_EB_R] - Bz_cells_avg;
            }
        }
    }
#endif



    /*
    double phi_max = theDomain->phi_max;

    int I0[Nr*Nz];
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = j + Nr*k;
            I0[jk] = 0;
            for(i=0; i<Np[jk]; i++)
            {
                int im = (i == 0) ? Np[jk]-1 : i - 1;

                double piph = get_dp(theCells[jk][i].piph, 0.0);
                double pimh = get_dp(theCells[jk][im].piph, 0.0);
                if(piph > 0.0 && pimh <= 0.0)
                {
                    I0[jk] = i;
                    break;
                }
            }
        }

    for(k=NgZa; k<Nz-NgZb; k++)
    {
        for(j=0; j<Nr-1; j++)
        {
            int jkL = j + Nr*k;
            int jkR = j+1 + Nr*k;

            int iL = I0[jkL];
            int iR = I0[jkR];

            int count;
            for(count=0; count < Np[jkL]+Np[jkR]; count++)
            {
                
            }
        }
    }
    */




}

void subtract_advective_B_fluxes( struct domain * theDomain ){

   int i,j,k;
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   if (NUM_EDGES >= 4)
   {
       for( k=0 ; k<Nz ; ++k ){
          for( j=0 ; j<Nr ; ++j ){
             int jk = j+Nr*k;
             double xp[3] = {r_jph[j], 0.0, z_kph[k]};
             double xm[3] = {r_jph[j-1], 0.0, z_kph[k-1]};
             double xc[3];
             get_centroid_arr(xp, xm, xc);

             double x[3] = {xm[0], 0.0, xc[2]};
             double hL = get_scale_factor(x, 0);
             x[0] = xp[0];
             double hR = get_scale_factor(x, 0);
             x[0] = xc[0]; x[2] = xm[2];
             double hD = get_scale_factor(x, 0);
             x[2] = xp[2];
             double hU = get_scale_factor(x, 0);

             for( i=0 ; i<Np[jk] ; ++i ){
                struct cell * c  = theCells[jk]+i;

                double wm = c->wiph;
                double wp = c->wiph;

                c->E[0] -= hL*wm * c->B[0];
                c->E[1] -= hR*wp * c->B[1];

                if( NUM_EDGES == 8 ){
                   c->E[4] += hD*wm * c->B[4];
                   c->E[5] += hU*wp * c->B[5];
                }
             }
          }
       }
   }
}


void check_flipped( struct domain * theDomain , int dim ){

   struct face * theFaces;
   int * fI;
   int Nf_t;

   if( dim==0 ){
      theFaces = theDomain->theFaces_1;
      fI = theDomain->fIndex_r;
      Nf_t = theDomain->N_ftracks_r;
   }else{
      theFaces = theDomain->theFaces_2;
      fI = theDomain->fIndex_z;
      Nf_t = theDomain->N_ftracks_z;
   }

   int n;
   int n0 = 0; 
   int jk;
   for( jk=0 ; jk<Nf_t ; ++jk ){
      for( n=0 ; n<fI[jk]-n0 ; ++n ){
         struct face * f  = theFaces + n0 + n; 
         int np = n+1; 
         if( np==fI[jk]-n0 ) np=0; 
         struct face * fp = theFaces + n0 + np;

         double pp,pm;
         if( f->LRtype==0 ){
            pm = f->L->piph;
         }else{
            pm = f->R->piph;
         }
         if( fp->LRtype==0 ){
            pp = fp->L->piph;
         }else{
            pp = fp->R->piph;
         }
         double dp = get_signed_dp( pp , pm );
         if( dp<0. ){ fp->flip_flag=1; }//printf("FLIPPED YO\n");}
      }    
      n0 = fI[jk];
   }

}

void get_phi_pointer( struct face * f , double ** P , double ** RK_P , int dim ){

   if( dim==0){
       if(NUM_FACES >= 3) {
          if( f->LRtype == 0 ){ 
             *P    = &(f->L->Phi[2]);
             *RK_P = &(f->L->RK_Phi[2]);
          }else{
             *P    = &(f->R->Phi[1]);
             *RK_P = &(f->R->RK_Phi[1]);
          }
       }
   }else {
       if(NUM_FACES >= 5) {
          if( f->LRtype == 0 ){ 
             *P    = &(f->L->Phi[4]);
             *RK_P = &(f->L->RK_Phi[4]);
          }else{
             *P    = &(f->R->Phi[3]);
             *RK_P = &(f->R->RK_Phi[3]);
          }  
       }
   }

}

void flip_fluxes( struct domain * theDomain , int dim ){

   struct face * theFaces;
   int * fI;
   int Nf_t;
   if( dim==0 ){
      theFaces = theDomain->theFaces_1;
      fI = theDomain->fIndex_r;
      Nf_t = theDomain->N_ftracks_r;
   }else{
      theFaces = theDomain->theFaces_2;
      fI = theDomain->fIndex_z;
      Nf_t = theDomain->N_ftracks_z;
   }

   int n;
   int n0 = 0;
   int jk;
   for( jk=0 ; jk<Nf_t ; ++jk ){
      for( n=0 ; n<fI[jk]-n0 ; ++n ){
         struct face * f  = theFaces + n0 + n;
         if( f->flip_flag ){

            int np = n+1;
            if( np==fI[jk]-n0 ) np=0;
            struct face * fp = theFaces + n0 + np;
            int nm = n-1;
            if( nm==-1 ) nm = fI[jk]-n0-1;
            struct face * fm = theFaces + n0 + nm;

            double *P = NULL;
            double *Pm, *Pp, *RK_P, *RK_Pm, *RK_Pp;

            get_phi_pointer( f  , &P  , &RK_P  , dim );
            get_phi_pointer( fm , &Pm , &RK_Pm , dim );
            get_phi_pointer( fp , &Pp , &RK_Pp , dim );

            double Phi = *P;
            double RK_Phi = *RK_P;

            *P      = *Pm    + Phi;
            *RK_P   = *RK_Pm + RK_Phi;
            *Pm     = -Phi;
            *RK_Pm  = -RK_Phi;
            *Pp     = *Pp    + Phi;
            *RK_Pp  = *RK_Pp + RK_Phi;

            f->LRtype  = !f->LRtype;
            fm->LRtype = !fm->LRtype;

         }
      }    
      n0 = fI[jk];
   }

}


int phi_switch( double dphi , double Pmax , int mode ){
    // Returns "sign" of dphi, taking periodicity into account
    //mode == 0: return 1 if dphi > 0.0, 0 otherwise
    //mode == 1: return 1 if dphi < 0.0, 0 otherwise

   while( dphi > .5*Pmax ) dphi -= Pmax;
   while( dphi <-.5*Pmax ) dphi += Pmax;
   if( mode == 1 ) dphi = -dphi;

   int LR = 0;
   if( dphi > 0.) LR = 1;

   return( LR );

}

int get_which4( double phi , double phiR , double phiU , double phiUR , int * LR_alt , int * UD_alt , int mode , double Pmax ){

    // Determine whether phi(0), phiR(1), phiU(2), or phiUR(3) is smallest,
    // taking into account periodicity, and return its code (0,1,2,or 3).
    //
    // LR_alt = the LR toggle for the top or bottom, whichever which4 isnt
    // UD_alt = the UD toggle for the left or right, whichever which4 isnt
    //
    // if mode == 1, returns code for largest phi, LR_alt and UD_alt not set.
    
    //comments apply to mode==0 case

   int which4;
   double dphi;
   
   dphi = phi - phiR;
   int LR_D = phi_switch( dphi , Pmax , mode );  // =1 if phi > phiR

   dphi = phiU - phiUR;
   int LR_U = phi_switch( dphi , Pmax , mode );  // =1 if phiU > phiUR

   double phi1 = phi;
   if( LR_D ) phi1 = phiR;  //phi1 is the smaller of phi, phiR
   double phi2 = phiU;
   if( LR_U ) phi2 = phiUR; //phi2 is the smaller of phiU, phiUR

   dphi = phi1-phi2;
   int UD = phi_switch( dphi , Pmax , mode );   // =1 if phi1 > phi2

   if( UD==0 ){
      if( LR_D==0 ) which4 = 0; // phi is smallest
      else which4 = 1;          // phiR is smallest
   }else{
      if( LR_U==0 ) which4 = 2; // phiU is smallest
      else which4 = 3;          // phiUR is smallest
   }

   if( mode == 0 ){
      if( which4 == 0 || which4 == 1 ) *LR_alt = LR_U;
      else                             *LR_alt = LR_D;

      if( which4 == 0 || which4 == 2 ){
         dphi = phiR - phiUR;
         *UD_alt = phi_switch( dphi , Pmax , mode );
      }else{
         dphi = phi - phiU;
         *UD_alt = phi_switch( dphi , Pmax , mode );
      }
   }

   return( which4 );
}


void make_edge_adjust( struct domain * theDomain , double dt ){

#if (NUM_FACES < 5) && (NUM_AZ_EDGES < 4)
    UNUSED(theDomain);
    UNUSED(dt);
    return;
#else
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   double Pmax = theDomain->phi_max;
   int j,k;

   int I0[Nr*Nz];
   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
         int found=0;
         int quad_prev=0;
         
         int i;
         for( i=0 ; i<Np[jk] && !found ; ++i ){
            struct cell * c = theCells[jk]+i;
            //TODO
            // Make this coordinate independent.
            double convert = 2.*M_PI/Pmax;
            double sn = sin(c->piph*convert);
            double cs = cos(c->piph*convert);
            if( sn>0. && cs>0. && quad_prev ){
               quad_prev = 0; 
               found = 1; 
               I0[jk] = i; 
            }else if( sn<0. && cs>0. ){
               quad_prev=1;
            }    
         }    
         if( !found ) I0[jk]=0;
      }    
   }

   for( k=0 ; k<Nz-1 ; ++k ){
      for( j=0 ; j<Nr-1 ; ++j ){
         int jk   = j  +Nr*k;
         int jkR  = j+1+Nr*k;
         int jkU  = j  +Nr*(k+1);
         int jkUR = j+1+Nr*(k+1);

         double xp[3] = {r_jph[j],0.0,z_kph[k]};
         double xm[3] = {r_jph[j],0.0,z_kph[k]};

         int i   = I0[jk  ];
         int iR  = I0[jkR ];
         int iU  = I0[jkU ];
         int iUR = I0[jkUR];
         
         double rL = get_centroid(r_jph[j],   r_jph[j-1], 1);
         double rR = get_centroid(r_jph[j+1], r_jph[j],   1);
         double zD = get_centroid(z_kph[k],   z_kph[k-1], 2);
         double zU = get_centroid(z_kph[k+1], z_kph[k],   2);

         int Ne = Np[jk] + Np[jkR] + Np[jkU] + Np[jkUR];
         int e;
         for( e=0 ; e<Ne ; ++e ){

            struct cell * c   = &(theCells[jk  ][i  ]);
            struct cell * cR  = &(theCells[jkR ][iR ]);
            struct cell * cU  = &(theCells[jkU ][iU ]);
            struct cell * cUR = &(theCells[jkUR][iUR]);

            int LR_alt;
            int UD_alt;
            int which4      = get_which4(   c->piph,
                                            cR->piph,
                                            cU->piph,
                                            cUR->piph, 
                                            &LR_alt , &UD_alt , 0 , Pmax );
            int which4_back = get_which4(   c->piph - c->dphi,
                                            cR->piph - cR->dphi,
                                            cU->piph-cU->dphi, 
                                            cUR->piph-cUR->dphi,
                                            NULL , NULL , 1 , Pmax );

            double * PhiL;
            double * PhiR;
            double * PhiU;
            double * PhiD;

            double E;

            if( which4 == 0 ){
               PhiL = c->Phi+4;
               PhiD = c->Phi+2;
               E = .25*( c->E_phi[3] + c->E_phi[1] );
               if( LR_alt==0 ){
                  PhiU = cU->Phi+2;
                  E += .25*cU->E_phi[1];
               }else{
                  PhiU = cUR->Phi+1;
                  E += .25*cUR->E_phi[0];
               }
               if( UD_alt==0 ){
                  PhiR = cR->Phi+4;
                  E += .25*cR->E_phi[3];
               }else{
                  PhiR = cUR->Phi+3;
                  E += .25*cUR->E_phi[2];
               }
            }else if( which4 == 1 ){
               PhiR = cR->Phi+4;
               PhiD = cR->Phi+1;
               E = .25*( cR->E_phi[3] + cR->E_phi[0] );
               if( LR_alt==0 ){
                  PhiU = cU->Phi+2;
                  E += .25*cU->E_phi[1];
               }else{
                  PhiU = cUR->Phi+1;
                  E += .25*cUR->E_phi[0];
               }
               if( UD_alt==0 ){
                  PhiL = c->Phi+4;
                  E += .25*c->E_phi[3];
               }else{
                  PhiL = cU->Phi+3;
                  E += .25*cU->E_phi[2];
               }
            }else if( which4 == 2 ){
               PhiL = cU->Phi+3;
               PhiU = cU->Phi+2;
               E = .25*( cU->E_phi[2] + cU->E_phi[1] );
               if( LR_alt==0 ){
                  PhiD = c->Phi+2;
                  E += .25*c->E_phi[1];
               }else{
                  PhiD = cR->Phi+1;
                  E += .25*cR->E_phi[0];
               }
               if( UD_alt==0 ){
                  PhiR = cR->Phi+4;
                  E += .25*cR->E_phi[3];
               }else{
                  PhiR = cUR->Phi+3;
                  E += .25*cUR->E_phi[2];
               }
            }else{
               PhiR = cUR->Phi+3;
               PhiU = cUR->Phi+1;
               E = .25*( cUR->E_phi[2] + cUR->E_phi[0] );
               if( LR_alt==0 ){
                  PhiD = c->Phi+2;
                  E += .25*c->E_phi[1];
               }else{
                  PhiD = cR->Phi+1;
                  E += .25*cR->E_phi[0];
               }
               if( UD_alt==0 ){
                  PhiL = c->Phi+4;
                  E += .25*c->E_phi[3];
               }else{
                  PhiL = cU->Phi+3;
                  E += .25*cU->E_phi[2];
               }
            } 

            if( which4 == 0 ){
               xp[1] = c->piph;
            }else if( which4 == 1 ){
               xp[1] = cR->piph;
            }else if( which4 == 2 ){
               xp[1] = cU->piph;
            }else{
               xp[1] = cUR->piph;
            }

            if( which4_back == 0 ){
               xm[1] = c->piph  - c->dphi;
            }else if( which4_back == 1 ){
               xm[1] = cR->piph - cR->dphi;
            }else if( which4_back == 2 ){
               xm[1] = cU->piph - cU->dphi;
            }else{
               xm[1] = cUR->piph- cUR->dphi;
            }

            // Gardiner & Stone adjustment
            if(CT_Solver == 1)
            {
                double Ec = 0.0;
                double dphi = get_dp(xp[1], xm[1]);
                double phi = xp[1] - 0.5*dphi;
                double x[3] = {xp[0], phi, xp[2]};

                int q;
                double prim[NUM_Q], Ecell[3];
                // Cell c
                double dphic = get_signed_dp(phi, c->piph-0.5*c->dphi);
                for(q=0; q<NUM_Q; q++)
                    prim[q] = c->prim[q] + dphic * c->gradp[q];
                prim_to_E(prim, Ecell, x);
                Ec += (x[0]-rL)*(x[2]-zD)*Ecell[1];
                // Cell cU
                dphic = get_signed_dp(phi, cU->piph-0.5*cU->dphi);
                for(q=0; q<NUM_Q; q++)
                    prim[q] = cU->prim[q] + dphic * cU->gradp[q];
                prim_to_E(prim, Ecell, x);
                Ec += (x[0]-rL)*(zU-x[2])*Ecell[1];
                // Cell cR
                dphic = get_signed_dp(phi, cR->piph-0.5*cR->dphi);
                for(q=0; q<NUM_Q; q++)
                    prim[q] = cR->prim[q] + dphic * cR->gradp[q];
                prim_to_E(prim, Ecell, x);
                Ec += (rR-x[0])*(x[2]-zD)*Ecell[1];
                // Cell cUR
                dphic = get_signed_dp(phi, cUR->piph-0.5*cUR->dphi);
                for(q=0; q<NUM_Q; q++)
                    prim[q] = cUR->prim[q] + dphic * cUR->gradp[q];
                prim_to_E(prim, Ecell, x);
                Ec += (rR-x[0])*(zU-x[2])*Ecell[1];
                Ec /= (rR-rL)*(zU-zD);

                E = 2*E-Ec;  //Gardiner & Stone adjustment (their Ez0 scheme)
            }
            

            double dl = get_dL( xp , xm , 0 );
//if( e==0 ) printf("dl = %e which4 = %d, which4_back = %d, phip = %e phim = %e dphi=%e \n",dl,which4,which4_back,xp[1],xm[1],xp[1]-xm[1]);
            add_E_phi( PhiL , PhiR , PhiD , PhiU , E*dl*dt );

            if( which4 == 0 ){
               ++i;
               if( i   == Np[jk  ] ) i  =0;
            }else if( which4 == 1 ){
               ++iR;
               if( iR  == Np[jkR ] ) iR =0;
            }else if( which4 == 2 ){
               ++iU;
               if( iU  == Np[jkU ] ) iU =0;
            }else{
               ++iUR;
               if( iUR == Np[jkUR] ) iUR=0;
            }
         }

      }
   }

#endif
}

