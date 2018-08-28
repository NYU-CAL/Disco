
#include "paul.h"

double get_moment_arm( double * , double * );
double get_dV( double * , double * );

int num_diagnostics( void );
void initializePlanets( struct planet * );
void init_tracerList( struct domain * );
void initializeTracers( struct domain * );

void setICparams( struct domain * );
void setHydroParams( struct domain * );
void setGeometryParams( struct domain * );
void setRiemannParams( struct domain * );
void setGravParams( struct domain * );
void setPlanetParams( struct domain * );
void setTracerParams( struct domain * );
void setHlldParams( struct domain * );
void setOmegaParams( struct domain * );
void setRotFrameParams( struct domain * );
void setMetricParams( struct domain * );
void setFrameParams(struct domain * );
void setDiagParams( struct domain * );
void setNoiseParams( struct domain * );
void setProcessCoords( struct domain * );
void build_tracerList( struct domain * );

int get_num_rzFaces( int , int , int );

static double t_start = 0.0;

void setupDomain( struct domain * theDomain ){

   srand(theDomain->rank);
   rand();

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   theDomain->theCells = (struct cell **) malloc( Nr*Nz*sizeof(struct cell *) );
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      theDomain->theCells[jk] = (struct cell *) malloc( Np[jk]*sizeof(struct cell) );
   }

   setGravParams( theDomain );
   setPlanetParams( theDomain );
   int Npl = theDomain->Npl;
   theDomain->thePlanets = (struct planet *) malloc( Npl*sizeof(struct planet) );
   initializePlanets( theDomain->thePlanets );

   //set tracer params
   setProcessCoords( theDomain );
   setGeometryParams( theDomain );
   setTracerParams( theDomain );

   int num_tools = num_diagnostics();
   theDomain->num_tools = num_tools;
   theDomain->theTools.t_avg = 0.0;
   theDomain->theTools.Qrz = (double *) malloc( Nr*Nz*num_tools*sizeof(double) );
   int i;
   for( i=0 ; i<Nr*Nz*num_tools ; ++i ) theDomain->theTools.Qrz[i] = 0.0;

   double Pmax = theDomain->theParList.phimax;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      double p0 = Pmax*(double)rand()/(double)RAND_MAX;
      double dp = Pmax/(double)Np[jk];
      for( i=0 ; i<Np[jk] ; ++i ){
         double phi = p0+dp*(double)i;
         if( phi > Pmax ) phi -= Pmax;
         theDomain->theCells[jk][i].piph = phi;
         theDomain->theCells[jk][i].dphi = dp;
      }
   }

   //Initialize Tracers
   theDomain->theTracers = (struct tracerList *) malloc( sizeof(struct tracerList) );
   init_tracerList( theDomain );  //initialize each processor's tracer list
   initializeTracers( theDomain );  //initialize tracers on each process
   //initTracers_Rand( theDomain );
  

   theDomain->t       = theDomain->theParList.t_min;
   theDomain->t_init  = theDomain->theParList.t_min;
   theDomain->t_fin   = theDomain->theParList.t_max;
   theDomain->phi_max = theDomain->theParList.phimax;

   theDomain->N_rpt = theDomain->theParList.NumRepts;
   theDomain->N_snp = theDomain->theParList.NumSnaps;
   theDomain->N_chk = theDomain->theParList.NumChecks;
   theDomain->N_xyz = theDomain->theParList.tr_out_num;

   theDomain->count_steps = 0;
   theDomain->final_step  = 0;
   theDomain->check_plz   = 0;

   theDomain->nrpt=-1;
   theDomain->nsnp=-1;
   theDomain->nchk=-1;
   theDomain->nxyz=-1;

   theDomain->theFaces_1 = NULL;
   theDomain->theFaces_2 = NULL;
   theDomain->N_ftracks_r = get_num_rzFaces( Nr , Nz , 1 );
   theDomain->N_ftracks_z = get_num_rzFaces( Nr , Nz , 2 );

   theDomain->fIndex_r = (int *) malloc( (theDomain->N_ftracks_r+1)*sizeof(int) );
   theDomain->fIndex_z = (int *) malloc( (theDomain->N_ftracks_z+1)*sizeof(int) );

   setICparams( theDomain );
   setHydroParams( theDomain );
   setRiemannParams( theDomain );
   setHlldParams( theDomain );
   setOmegaParams( theDomain );
   setRotFrameParams( theDomain );
   setMetricParams( theDomain );
   setFrameParams( theDomain );
   setDiagParams( theDomain );
   setNoiseParams( theDomain );

}

void initial( double * , double * );
void prim2cons( double * , double * , double * , double );
void cons2prim( double * , double * , double * , double );
void restart( struct domain * );
void calc_dp( struct domain * );
void set_wcell( struct domain * );
void adjust_gas( struct planet * , double * , double * , double );
void set_B_fields( struct domain * );
void subtract_omega( double * );
void addNoise(double *prim, double *x);

void setupCells( struct domain * theDomain ){

   int restart_flag = theDomain->theParList.restart_flag;
   if( restart_flag ) restart( theDomain );
   t_start = theDomain->t;

   int noiseType = theDomain->theParList.noiseType;

   calc_dp( theDomain );

   int i,j,k;
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int atmos = theDomain->theParList.include_atmos;

   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            c->wiph = 0.0;
            double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double r = get_moment_arm( xp , xm );
            double dV = get_dV( xp , xm );
            double phi = c->piph-.5*c->dphi;
            double x[3] = { r , phi , .5*(z_kph[k]+z_kph[k-1])};
            if( !restart_flag )
            {
               initial( c->prim , x );
               subtract_omega( c->prim );
               if( atmos ){
                  int p;
                  for( p=0 ; p<Npl ; ++p ){
                     double gam = theDomain->theParList.Adiabatic_Index;
                     adjust_gas( theDomain->thePlanets+p , x , c->prim , gam );
                  }
               }
            }
            if(noiseType != 0)
                addNoise(c->prim, x);            
            prim2cons( c->prim , c->cons , x , dV );
            cons2prim( c->cons , c->prim , x , dV );
            c->real = 1;
            if( xm[0]==0.0 ) //if cell is a middle wedge
                c->origin = 1;
            else
                c->origin = 0;
         }    
      }    
   }
   set_wcell( theDomain );

}


/*
void clear_cell( struct cell * c ){
   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      c->prim[q]   = 0.0;
      c->cons[q]   = 0.0;
      c->RKcons[q] = 0.0;
      c->grad[q]   = 0.0;
      c->gradr[q]  = 0.0;
   }
   c->riph = 0.0;
   c->RKriph = 0.0;
   c->dr = 0.0;
   c->wiph = 0.0;
}
*/

void freeDomain( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      free( theDomain->theCells[jk] );
   }
   free( theDomain->theCells );
   free( theDomain->Np );
   theDomain->r_jph--;
   free( theDomain->r_jph );
   theDomain->z_kph--;
   free( theDomain->z_kph );
   free( theDomain->thePlanets );
   free( theDomain->theTools.Qrz );
   free( theDomain->fIndex_r );
   free( theDomain->fIndex_z );

   int n;
   struct tracer *tr = theDomain->theTracers->head;
   while( tr!=NULL ){
        struct tracer *temp = tr;
        tr = tr->next;
        free(temp);
   }
   free( theDomain->theTracers );
}

void check_dt( struct domain * theDomain , double * dt ){

   double t = theDomain->t;
   double tmax = theDomain->t_fin;
   int final=0;
   int check=0;
   if( t + *dt > tmax ){
      *dt = tmax-t;
      final=1;
   }

   if( theDomain->rank==0 ){
      FILE * abort = NULL;
      abort = fopen("abort","r");
      if( abort ){ final = 1; fclose(abort); }
      FILE * latest = NULL;
      latest = fopen("latest","r");
      if( latest ){ check = 1; fclose(latest); remove("latest");}
   }

   MPI_Allreduce( MPI_IN_PLACE , &final , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , &check , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
   if( final ) theDomain->final_step = 1;
   if( check ) theDomain->check_plz = 1;
   
}

void report( struct domain * );
void snapshot( struct domain * , char * );
void output( struct domain * , char * );
void tracerOutput( struct domain * );
void tracerReport( struct domain * );

int possiblyOutput( struct domain * theDomain , int override ){

   double t = theDomain->t;
   double t_min = theDomain->t_init;
   double t_fin = theDomain->t_fin;
   double Nrpt = theDomain->N_rpt;
   double Nsnp = theDomain->N_snp;
   double Nchk = theDomain->N_chk;
   double Nxyz = theDomain->N_xyz;
   int LogOut = theDomain->theParList.Out_LogTime;
   int step = theDomain->mdStep;
   int n0;
   
   int dt_flag = 0;
   n0 = (int)( t*Nrpt/t_fin );
   if( LogOut ) n0 = (int)( Nrpt*log(t/t_min)/log(t_fin/t_min) );
   if( theDomain->nrpt < n0 || override ){
      theDomain->nrpt = n0;
      //longandshort( &theDomain , &L , &S , &iL , &iS , theDomain.theCells[0] , 0 , 0 );
      report( theDomain );

      if( theDomain->rank==0 ){ 
         printf("t = %.3e   Step = %d\n", t,step);
         dt_flag = 1;
      }
   }
   n0 = (int)( t*Nchk/t_fin );
   if( LogOut ) n0 = (int)( Nchk*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nchk < n0 && Nchk>0) || override || theDomain->check_plz ){
      theDomain->nchk = n0;
      char filename[256];
      if( !override ){
         if( !theDomain->check_plz ){
            if(theDomain->rank==0) printf("Creating Checkpoint #%04d...\n",n0);
            sprintf(filename,"./checkpoint_%04d",n0);
         }else{
            if(theDomain->rank==0) printf("Creating Requested Checkpoint...\n");
            sprintf(filename,"checkpoint_latest");
            theDomain->check_plz = 0;
         }
         output( theDomain , filename );
      }else{
         if(theDomain->rank==0) printf("Creating Final Checkpoint...\n");
         output( theDomain , "output" );
      }
   }
   n0 = (int)( t*Nsnp/t_fin );
   if( LogOut ) n0 = (int)( Nsnp*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nsnp < n0 && Nsnp>0) || override ){
      theDomain->nsnp = n0;
      char filename[256];
      if(!override) sprintf( filename , "snapshot_%04d" , n0 );
      else sprintf( filename , "snapshot" );
      //snapshot( theDomain , filename );
   }

   //Tracer Outputting
   n0 = (int)( (t-t_start)*Nxyz/(t_fin-t_start) );
   if( LogOut ) 
      n0 = (int)( Nxyz*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nxyz < n0 && Nxyz>0) || override ){
      theDomain->nxyz = n0;
      tracerOutput( theDomain );
      //tracerReport( theDomain );
   }

   return dt_flag;
}

int trOutStep( struct domain *theDomain ){

   int step = theDomain->mdStep;
   int printStep = theDomain->theParList.tr_out_step;

   if( (step+printStep)%printStep==0 )
      return 1;
   else
      return 0;
}

int getListSize( struct tracerList * );

void tracerOutput( struct domain *theDomain ){

   char filename[256];
   sprintf(filename, "%s.xyz", "tracers" );

   int Ntr_tot = theDomain->Ntr;
   //int step = theDomain->mdStep;
   int step = theDomain->nxyz;
   int rank = theDomain->rank;
   int size = theDomain->size;
   double t = theDomain->t;

   int rewrite_flag = theDomain->theParList.restart_flag;
   if( step )
       rewrite_flag = 0;
   int outFlag = theDomain->theParList.tr_out_flag;

   MPI_Allreduce( MPI_IN_PLACE, &Ntr_tot, 1, MPI_INT, MPI_SUM, theDomain->theComm );
/*
   if( rank==0 && step==0 ){
	   FILE * pFile = fopen(filename, "w");
      printf("Rewriting Tracer Out-File\n");
	   fclose(pFile);
   }
*/
   if( !rank ){
      if( !step || rewrite_flag ){
        FILE * pFile = fopen(filename, "w");
        fclose( pFile );
      }
   }
   MPI_Barrier( theDomain->theComm);

   int rk;
   int count=0;
   for( rk=0; rk<size; ++rk){
      //MPI_Barrier( theDomain->theComm );
      if( rank==rk ){
         FILE * pFile = fopen(filename, "a");
         if( rank==0 && outFlag==0 ){
            fprintf(pFile, "%d \nAtoms. Timestep: %d\n", Ntr_tot, step);
            //printf(" Ntr_tot from output: %d \n", Ntr_tot );
            // step,time, id,tpe,  x, y, z, r, phi, vr, om, vz
            //fprintf(pFile, "%d %4.2f %d %d %.2f %.2f %.2f %.2f %.2f \n",
            //         step,t, 0,0, 0.0,0.0,0.0, 0.0,0.0 );
         }
         struct tracer *tr = theDomain->theTracers->head;
         while( tr != NULL){
            int id   = tr->ID;
   	        int type = tr->Type;
   	        double r = tr->R;
      	    double phi = tr->Phi;
      	    double z = tr->Z;
      	    double x = r*cos(phi);
         	double y = r*sin(phi);
         	double vr = tr->prim[URR];
      	    double om = tr->prim[UPP];
         	double vz = tr->prim[UZZ];
            double rho = tr->prim[RHO];
            double ppp = tr->prim[PPP];
            //double dPdp = tr->dPdp;
            double C_j  = tr->C_j;
            double dCj  = tr->dPdp;
            //double gravfrac = tr->gravfrac;
            //double Lp = tr->Lp;
            //double Jj = tr->Jj;
            //double Ss = tr->Ss;
      	    fprintf(pFile, "%d %4.2f %d %d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n",    
                           step,t, 
                           id,type, 
                           x,y,z, 
                           vr,om,vz,
                           rho, ppp,
                           C_j, dCj);
            tr = tr->next;
            count++;
         }
         fclose(pFile);
      }
      MPI_Barrier( theDomain->theComm );
   }

}
