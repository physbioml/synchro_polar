#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#define ETA 0  /* left-right simm break term */
#define K  2.0   /* spring */
#define R 1.5
#define N_ROW 30
#define N_STEP 3000*1240
//#define N_STEP 1000

#define N_CONF 20000
#define XMQ 1


#define SW 0.244186047 //0.5*(PASSO/2)  /* switch referred to 1 */
#define D 6
#define VISC 7.5 /* viscosity */
#define OSEENA K/(4*M_PI*VISC)  
#define STOKES K/(6*M_PI*VISC*R)
#define ERRE 1  /* interpolating linear-quadratic  potential*/
#define N_PASSI_CAMPIONAMENTO_PER_CICLO 40000


struct complex{
    
    double Im;
    double Re;
    
};


struct particle{
    
    double theta;
    double phi;
    double x;
    double y;
    
};


/*----- distance ----*/
double compute_distance(double , double, double , double);
/*----- force phase ----*/
double compute_force_phase(double J0, double J1, double J2, double vv, double vvv, double uu, double uuu, double dist);
/*----- force angle ----*/
double compute_force_angle(double K0, double K1, double K2, double vv, double vvv, double uu, double uuu, double dist);
/*----- pos ----*/
double compute_pos(double v, double dt);
/*----- force ----*/
double compute_force(double o,  double dt );


FILE *output, *xt, *phasext, *orderparam, *poincare,*state,*outen, *photo;		  	
int main() 
{
  struct particle data[N_ROW], dataP[N_ROW], dataScratch[N_ROW], dataNew[N_ROW];
  struct complex S_m, S_p, ord_par_theta, ord_par_phi;
  int i,j,k,l,m,t,u,z,nswitch, freq_switch, delta_switch,gt;
  int sum,diff, strobo;
  int  memo[N_ROW], cnt1, cnt2;
  
  //float rando, randomv;
  double fact, fact2, fact_sym, fact2_sym;
    
  char fileoutxyz[150],fileoutxt[150],fileouphase[150] ,fileoutord[150], fileoutpoincare[150],fileoutxtphoto[150];
  double xmq_diag, xmq_non_diag, noise_diag, noise_non_diag,nnu,nnnu,nu[N_ROW],omega, deltas, ordampl, ordphase, ordampl2, ordphase2, devph[N_ROW],  devamp[N_ROW], devphase, devamplit;
  
    double theta[N_ROW], phi[N_ROW], x[N_ROW], y[N_ROW];
    double v0, erre, vvv, vv, uuu, uu, xxx, xx, oo, oo2, dt, yyy, yy, tt, pp, dist, d, rnd;
    double  o[N_ROW],  o2[N_ROW] ;
    double J1, J2, K1, K2, J0, K0;
    double deltax, deltay, deltatheta, deltaphi;
    double kappa;
    double ord_theta, ord_phi, ord_S_p, ord_S_m;
    
    double pippo;
    
  /*  random  gen  */
  gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs0);
  /* sid setting */
  gsl_rng_set(r,120000000);

    strobo =500;
    dt = 0.0001;
    v0  = 1.0;
    erre = 10.0;
    kappa = 10.0;
    
    J0 = 0.0;
    K0 = 0.0;
    
    J1 = 2.;
    J2 = 0.;
  
    K1 = 0.01;
    K2 = 0.;
    
    d=10.0;
    
      /* starting positions */
	for(k=0; k< N_ROW; k++) {
        
		data[k].theta = (2*M_PI*k)/N_ROW + gsl_ran_gaussian(r,XMQ);
		data[k].phi = (2*M_PI*k)/N_ROW + 0.5+ gsl_ran_gaussian(r,XMQ);
        pippo = gsl_ran_gaussian(r,XMQ);
        
        data[k].x = erre*cos((2*M_PI*k)/N_ROW );
        data[k].y = erre*sin((2*M_PI*k)/N_ROW );

		printf(" j %d\t x  %lf\t y %lf \n ", k, data[k].x, data[k].y );
		
    }

      /* external file names */
        sprintf(fileoutord, "%d-particles_ord_par_RK.dat",N_ROW );
        orderparam= fopen (fileoutord, "w");

        /* external file names */
        sprintf(fileoutxt,"%d-particles-pos_RK.dat",N_ROW);
        xt= fopen (fileoutxt, "w");
                   
      cnt1=0;
      cnt2=0;
     
      nswitch=0; 
      for (i=0;i<=N_STEP; i++){/* time-cycle */
          for(k=0; k< N_ROW; k++){
              o[k]=0;	/* prepare oseen*/
              o2[k]=0;
          }
          
          ord_par_theta.Re = 0.;
          ord_par_theta.Im = 0.;
          ord_par_phi.Re = 0.;
          ord_par_phi.Im = 0.;
          
          S_m.Re = 0. ;
          S_m.Im = 0. ;
          S_p.Re = 0. ;
          S_p.Im = 0. ;

          
          for(k=0; k< N_ROW; k++){/* starting k-rowers cycle*/
              vvv =  data[k].theta ;
              uuu =  data[k].phi ;

              for(l=0; l< N_ROW;l++){/* internal cycle: l!=k */
                  if(l > k){
                      
                      vv = data[l].theta;
                      uu  =  data[l].phi;

                      //dist = compute_distance(xx, xxx, yy, yyy); /* truncated distance */
                      dist = d; /* truncated distance */

                      fact = compute_force_phase( J0,  J1,  J2,  vv,  vvv,  uu,  uuu, dist);
                      fact2 = compute_force_angle( K0,  K1,  K2,  vv,  vvv,  uu,  uuu, dist);
                    
                      fact_sym = compute_force_phase( J0,  J1,  J2,  vvv,  vv,  uuu,  uu, dist);
                      fact2_sym = compute_force_angle( K0,  K1,  K2,  vvv,  vv,  uuu,  uu, dist);

                      /*add contributes to oseen non diagonal term*/
                      o[k] += (fact + fact_sym);
                      o2[k] += (fact2 + fact2_sym);
            
                  }/* end if */
	    
              }/* end l-cycle */
	  
              if(o[k]*dt > 1)  printf(" oo dt su part %d e' %f, dt o int troppo grandi\n", k, o[k]*dt), exit(0);
              if(o2[k]*dt > 1)  printf(" oo2 dt su part %d e' %f, dt o int troppo grandi\n", k, o2[k]*dt), exit(0);
              
              xx = data[k].x;
              yy = data[k].y;
              
              tt = data[k].theta;
              pp = data[k].phi;
              
              oo = o[k];
              oo2 = o2[k];

              /* molecular dynamics on positions and angles*/
              
              //deltatheta = compute_force(oo, dt)*dt  ;
              //deltaphi=  compute_force(oo2, dt)*dt ;
              //deltax = compute_pos(v0*cos(data[k].theta), dt)*dt;
              //deltay = compute_pos(v0*sin(data[k].theta), dt)*dt;
              
              dataP[k].theta = compute_angle(oo, dt)*dt  ;
              dataP[k].phi = compute_angle(oo2, dt)*dt  ;
              
              dataP[k].x = compute_pos(v0*cos(data[k].theta), dt)*dt;
              dataP[k].y = compute_pos(v0*sin(data[k].theta), dt)*dt;
   
              dataScratch[k].theta = 0.5*dataP[k].theta + tt;
              dataScratch[k].phi = 0.5*dataP[k].phi + pp;
              dataScratch[k].x = 0.5*dataP[k].x + xx;
              dataScratch[k].y = 0.5*dataP[k].y + yy;

              dataNew[k].theta =  tt + compute_angle(dataScratch[k]., PosForceScratch[k].x )*dt ;

              dataNew[k].x =  xx + compute_pos(v0*cos(dataScratch[k].theta), dt)*dt;
              dataNew[k].y =  yy + compute_pos(v0*sin(dataScratch[k].theta), dt)*dt;
              
              /* increments*/
              xx += deltax ;
              yy += deltay ;
              
              tt += deltatheta ;
              pp += deltaphi ;
              
              if(isnan(xx) == 1) printf("NANNA x!"), exit(0);
              if(isnan(yy) == 1) printf("NANNA y!"), exit(0);
              if(isnan(tt) == 1) printf("NANNA theta!"), exit(0);
              if(isnan(pp) == 1) printf("NANNA phi!"), exit(0);
              
              if(deltatheta > 0.3) printf("deltatheta = %g troppo!\n", deltatheta), exit(0);
              if(deltaphi > 0.3) printf("deltaphi = %g troppo!\n", deltaphi), exit(0);
              
              data[k].x = xx;
              data[k].y = yy;
              data[k].theta = tt;
              data[k].phi = pp;
              
              ord_par_theta.Re += cos(theta[k]);
              ord_par_phi.Re  += cos(phi[k]) ;
              ord_par_theta.Im += sin(theta[k]);
              ord_par_phi.Im  += sin(phi[k]) ;
              
              S_p.Re += cos(theta[k] + phi[k]) ;
              S_p.Im += sin(theta[k] + phi[k]);
              S_m.Re += cos(theta[k] - phi[k]) ;
              S_m.Im += sin(theta[k] - phi[k]) ;

              
          }/* end k-cycle */
	
		  ord_theta = (ord_par_theta.Re*ord_par_theta.Re + ord_par_theta.Im*ord_par_theta.Im)/((double)N_ROW*(double)N_ROW);
		  ord_phi = (ord_par_phi.Re*ord_par_phi.Re + ord_par_phi.Im*ord_par_phi.Im)/((double)N_ROW*(double)N_ROW);
          
          ord_S_p =(S_p.Re*S_p.Re + S_p.Im*S_p.Im)/((double)N_ROW*(double)N_ROW);
          ord_S_m =(S_m.Re*S_m.Re + S_m.Im*S_m.Im)/((double)N_ROW*(double)N_ROW);

          
if(i%strobo==0) fprintf(orderparam, "%lf\t %lf\t %lf\t %lf\t %lf\n", 0.001*dt*((double)i), ord_theta, ord_phi, ord_S_p, ord_S_m);

          if(i%strobo==0){
              fprintf(xt, "%lf\t", 0.001*dt*((double)i));
              for(k=0;k<9;k++) fprintf(xt, "%lf\t  %lf\t",  data[k].x, data[k].y );
              fprintf(xt, "%lf\t  %lf\n",  data[9].x, data[9].y );
          }
          
	      }/* end time-cycle*/
      
      printf("nswitch %d\n",nswitch);	       
    
  gsl_rng_free(r); /* free rnd generator */
  return (0);
}


/************** functions *************/

/*----- distance ----*/
double compute_distance(double xxx, double xx, double yyy, double yy){
    
    double dist;
    dist  = (xx - xxx)*(xx - xxx) + (yy - yyy)*(yy - yyy) ;
    
    return sqrt(dist);
    
}

/*----- force phase ----*/
double compute_force_phase(double J0, double J1, double J2, double vv, double vvv, double uu, double uuu, double dist){
    double fact;

    fact  = J0*sin(vv-vvv)  + J1*sin(uu-uuu)*sin(vv-vvv)  +  J2*cos(uu-uuu)*sin(vv-vvv);
    fact = fact/(dist*dist);

    return fact;

}

/*----- force angle ----*/

double compute_force_angle(double K0, double K1, double K2, double vv, double vvv, double uu, double uuu, double dist){

    double fact;

    fact  = K0*sin(uu-uuu) + K1*sin(vv-vvv)*sin(uu-uuu) + K2*cos(vv-vvv)*sin(uu-uuu);
    fact = fact/(dist*dist);

    return fact;

}



/*----- force ----*/
double compute_angle(double o,  double dt ){
    
    double angle;
    
    angle = o;
    if(angle*dt > 0.3) printf("force = %g troppo!\n", angle), exit(0);
    if(isnan(angle) == 1) printf("force NANNA!"), exit(0);
    
    return angle;
}



/*----- pos ----*/
double compute_pos(double v, double dt){
    
    double pos;
    
    pos  = v;
    if(pos*dt > 0.3) printf("position = %g troppo!\n", pos), exit(0);
    if(isnan(pos) == 1) printf("position NANNA!"), exit(0);
    
    return pos;
}

/*------ phase  -----*/
double compute_phase(double s, double x){
    
    double phase;
    
    phase = atan2(s, x);
    
    return phase;
}

/*------ amplitude  -----*/
double compute_amplitude(double s, double x){
    
    double amplitude;
    
    amplitude = x*x + s*s;
    amplitude = sqrt(amplitude);
    
    return amplitude;
}
