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

#define N_CONF 20000
#define XMQ 1


#define SW 0.244186047 //0.5*(PASSO/2)  /* switch referred to 1 */
#define D 6
#define VISC 7.5 /* viscosity */
#define OSEENA K/(4*M_PI*VISC)  
#define STOKES K/(6*M_PI*VISC*R)
#define ERRE 1  /* interpolating linear-quadratic  potential*/
#define N_PASSI_CAMPIONAMENTO_PER_CICLO 40000

FILE *output, *xt, *phasext, *orderparam, *poincare,*state,*outen, *photo;		  	
int main() 
{
  int i,j,k,l,m,t,u,z,nswitch, freq_switch, delta_switch,gt;                    
  int sum,diff, strobo;
  int  memo[N_ROW], cnt1, cnt2;
  
 // double  f[N_ROW], sigma[N_ROW], ss,phase[N_ROW],amplitude[N_ROW], mu;
  //double xtracer[2],vtracer[2];
  //double   fact_imm, dist, distq, dist_imm, d_bordoq,rumore_nondiag,rumore_diag;
  //double vtracerx,vtracery,xtracerx,xtracery;
    
  //double  alpha, beta;
  //double teta;
  //double rnd,var_rumore, d;
  //double eps,nx,ny,rho,sigmax,sigmay,ciccio,nnx,nny;
  float rando, randomv;
  char fileoutxyz[150],fileoutxt[150],fileouphase[150] ,fileoutord[150], fileoutpoincare[150],fileoutxtphoto[150];
  double xmq_diag, xmq_non_diag, noise_diag, noise_non_diag,nnu,nnnu,nu[N_ROW],omega, deltas, ordampl, ordphase, ordampl2, ordphase2, devph[N_ROW],  devamp[N_ROW], devphase, devamplit;
  
    double theta[N_ROW], phi[N_ROW], x[N_ROW], y[N_ROW];
    double v0, erre, vvv, vv, uuu, uu, xxx, xx, oo, oo2, dt, fact, fact2, yyy, yy, tt, pp, dist, d, rnd;
    double  o[N_ROW],  o2[N_ROW] ;
    double J1, J2, K1, K2, J0, K0;
    double deltax, deltay, deltatheta, deltaphi;
    double kappa;
    double Re_ord_par_theta, Im_ord_par_theta, Re_ord_par_phi, Im_ord_par_phi, Re_S_p, Re_S_m,
    Im_S_p, Im_S_m;
    
    double ord_par_theta, ord_par_phi, S_p, S_m;

  /*  random  gen  */
  gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs0);
 
  /* sid setting */
  gsl_rng_set(r,120000000); 
    
  // for(gt=0;gt<5;gt++){ /* feedback frequency cycle*/
   // delta_switch=delta_switch*2;
   // printf("oseen%f\t stokes %f\n",OSEENA,STOKES);
    
    //d=15;
    //mu=0.01;
	
   // for(u=0;u<6;u++){/* distance cycle */
     // mu = mu/(double)10;    
     // printf("mu %lf\n",mu);
      //dt=2*M_PI*sqrt(6*M_PI*VISC*R)/(10*(double)N_PASSI_CAMPIONAMENTO_PER_CICLO);
	//printf("dt %lf\n",dt);

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
    double pippo;
    
      /* starting positions */
	for(k=0; k< N_ROW; k++) {
		theta[k] = (2*M_PI*k)/N_ROW + gsl_ran_gaussian(r,XMQ);
		phi[k] = (2*M_PI*k)/N_ROW + 0.5+ gsl_ran_gaussian(r,XMQ);
        pippo = gsl_ran_gaussian(r,XMQ);
        
        x[k] = erre*cos((2*M_PI*k)/N_ROW );
        y[k] = erre*sin((2*M_PI*k)/N_ROW );

		// = gsl_rng_uniform(r);
		// = gsl_rng_uniform(r);
		printf(" j %d\t x  %lf\t y %lf \n ", k, x[k], y[k] );
		
    }

      /* external file names */
        sprintf(fileoutord, "%d-particles_ord_par.dat",N_ROW );
        orderparam= fopen (fileoutord, "w");

        /* external file names */
        sprintf(fileoutxt,"%d-particles-pos.dat",N_ROW);
        xt= fopen (fileoutxt, "w");
                   
      cnt1=0;
      cnt2=0;
     
      nswitch=0; 
      for (i=0;i<=N_STEP; i++){/* time-cycle */
	
	//	gsl_ran_bivariate_gaussian(r,sigmax,sigmay,rho, &nx, &ny);
	//nu[0]=nx;
	//nu[1]=ny;
	//sigma[0]= sin(omega*(double)i*dt);
	//sigma[1]= sin(omega*(double)i*dt + phi);
	/* prepare oseen*/
          
          for(k=0; k< N_ROW; k++){
              o[k]=0;
              o2[k]=0;
          }
          
	for(k=0; k< N_ROW; k++){/* starting k-rowers cycle*/
	  vvv =  theta[k] ;
       uuu =  phi[k] ;

	  for(l=0; l< N_ROW;l++){/* internal cycle: l!=k */	    
	    if(l != k){ 
	      //xx = sqrt(x[l]*x[l] + y[l]*y[l]);
	      vv =theta[l];
	      uu  = phi[l];
	      /* force wih model parameters  F = (r*f -sigma*sw)   */
	      fact  = J0*sin(vv-vvv)  + J1*sin(uu-uuu)*sin(vv-vvv)  +  J2*cos(uu-uuu)*sin(vv-vvv)  ;
          fact2  = K0*sin(uu-uuu) + K1*sin(vv-vvv)*sin(uu-uuu) + K2*cos(vv-vvv)*sin(uu-uuu)  ;
			/* this is the asynchrony*/
	      //fact +=  (1/(double)3)*mu*xx*xx*xx;      
	      /* (k,l) distance */	      
	      //dist  = ((xx + d*(double)l) - (xxx + d*(double)k));
	      dist = d; /* truncated distance */
	      fact = (fact)/fabs(dist*dist);
          fact2 = (fact2)/fabs(dist*dist);

            /*add contributes to oseen non diagonal term*/
	      o[k] += (fact); 	    
          o2[k] += (fact2);
            
	    }/* end l-cycle */
	    
	  }/* end if */	
	  
	  if(o[k]*dt > 1)  printf(" oo dt su part %d e' %f, dt o int troppo grandi\n", k, o[k]*dt), exit(0);
     if(o2[k]*dt > 1)  printf(" oo2 dt su part %d e' %f, dt o int troppo grandi\n", k, o2[k]*dt), exit(0);
 
	}/* end k-cycle */
	
	/*cycle again*/
          Re_ord_par_theta = 0;
          Re_ord_par_phi = 0;
          Im_ord_par_theta = 0;
          Im_ord_par_phi = 0;
          Re_S_p = 0;
          Im_S_p = 0;
          Re_S_m = 0;
          Im_S_m = 0;

	for(j=0;j<N_ROW;j++){
	  /*move state, position, oseen in dummy variables*/
      xxx = sqrt(x[k]*x[k] + y[k]*y[k]) ;

      xx = x[j];
      yy = y[j];

     tt = theta[j];
     pp = phi[j];

	  oo = o[j];
      oo2 = o2[j];
	  //oo=0;
	  //nnnu = nu[j];
        
	  /* molecular dynamics on positions and angles*/
	  deltatheta = oo*dt   ;
	  deltaphi=  oo2*dt ;
	  
     deltax = v0*(cos(theta[j]))*dt;
     deltay = v0*(sin(theta[j]))*dt;

        /* increments*/

     xx += deltax ;
     yy += deltay ;
 
        tt += deltatheta ;
        pp += deltaphi ;

        if(isnan(xx) == 1) printf("NANNA x!"), exit(0);
        if(isnan(yy) == 1) printf("NANNA y!"), exit(0);
        if(isnan(tt) == 1) printf("NANNA theta!"), exit(0);
        if(isnan(pp) == 1) printf("NANNA phi!"), exit(0);
        
	  if(deltax > 0.3) printf("deltax = %g troppo!\n", deltax), exit(0);
      if(deltay > 0.3) printf("deltay = %g troppo!\n", deltay), exit(0);
      if(deltatheta > 0.3) printf("deltatheta = %g troppo!\n", deltatheta), exit(0);
      if(deltaphi > 0.3) printf("deltaphi = %g troppo!\n", deltaphi), exit(0);

      x[j] = xx;
	  y[j] = yy;
      theta[j] = tt;
      phi[j] = pp;
        
	  
        //if(i > 100*40000 && i%strobo==0){
//		  		if(i > 2998*40000) {
//			fprintf(xt, "%d\t %lf\t %lf\n", j, 0.001*dt*((double)i), x[j]);
//			fprintf(phasext, "%d\t %lf\t %lf\n", j, 0.001*dt*((double)i), phase[j]);
//					}

        //Re_ord_par_theta, Im_ord_par_theta, Re_ord_par_phi, Im_ord_par_phi, S_p, S_m
        
        Re_ord_par_theta += cos(theta[j]);
        Re_ord_par_phi  += cos(phi[j]) ;
        Im_ord_par_theta += sin(theta[j]);
        Im_ord_par_phi  += sin(phi[j]) ;
        Re_S_p += cos(theta[j] + phi[j]) ;
        Im_S_p += sin(theta[j] + phi[j]);
        Re_S_m += cos(theta[j] - phi[j]) ;
        Im_S_m += sin(theta[j] - phi[j]) ;
        
	}	/* end j-cycle */
          
//		if(i > 2998*40000) {
//		  fprintf(xt, "\n");
//		  fprintf(phasext, "\n");
//		  	  	  }
	 
		  ord_par_theta = (Re_ord_par_theta*Re_ord_par_theta + Im_ord_par_theta*Im_ord_par_theta)/((double)N_ROW*(double)N_ROW);
		  ord_par_phi = (Re_ord_par_phi*Re_ord_par_phi + Im_ord_par_phi*Im_ord_par_phi)/((double)N_ROW*(double)N_ROW);
          S_p =(Re_S_p*Re_S_p + Im_S_p*Im_S_p)/((double)N_ROW*(double)N_ROW);
          S_m =(Re_S_m*Re_S_m + Im_S_m*Im_S_m)/((double)N_ROW*(double)N_ROW);

          
if(i%strobo==0)
    fprintf(orderparam, "%lf\t %lf\t %lf\t %lf\t %lf\n", 0.001*dt*((double)i), ord_par_theta, ord_par_phi, S_p, S_m);

if(i%strobo==0) fprintf(xt, "%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n", 0.001*dt*((double)i), x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3], x[4], y[4], x[5], y[5], x[6], y[6], x[7], y[7], x[8], y[8], x[9], y[9]);

          
	/* FILE PRINTING */      
	/* xt print */
	/* print according to the frame rate */
	      }/* end time-cycle*/
      
      printf("nswitch %d\n",nswitch);	       
      
      /*free & fclose*/    
      //fclose (xt);
	// fclose (phasext);
  //    } /* uu-cycle */
      //} /* uuu-cycle */
	
  gsl_rng_free(r); /* free rnd generator */
  return (0);
}
