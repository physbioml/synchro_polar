#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#define ETA 0  /* left-right simm break term */
#define K  2.0   /* spring */
#define R 1.5
#define N_ROW 100
#define N_STEP 3000*40000

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
  int  memo[N_ROW],uu, cnt1, cnt2;
  
  double x[N_ROW], f[N_ROW], sigma[N_ROW], ss,phase[N_ROW],amplitude[N_ROW], mu;	              
  double v[N_ROW], o[N_ROW] ;
  double xtracer[2],vtracer[2];
  double  fact, fact_imm, dist, distq, dist_imm, d_bordoq,rumore_nondiag,rumore_diag;
  double xx, xxx, vv,vvv, oo, dt;
  double vtracerx,vtracery,xtracerx,xtracery;
  double deltax;                        
  double  alpha, beta;
  double teta;  
  double rnd,var_rumore, d;
  double eps,nx,ny,rho,sigmax,sigmay,ciccio,nnx,nny;
  float rando, randomv;
  char fileoutxyz[150],fileoutxt[150],fileouphase[150] ,fileoutord[150], fileoutpoincare[150],fileoutxtphoto[150];
  double xmq_diag, xmq_non_diag, noise_diag, noise_non_diag,nnu,nnnu,nu[N_ROW],omega, phi[N_ROW], deltas, ordampl, ordphase, ordampl2, ordphase2, devph[N_ROW],  devamp[N_ROW], devphase, devamplit;
  
  /*  random  gen  */
  gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs0);
 
  /* sid setting */
  gsl_rng_set(r,120000000); 
  
  eps=1;
  delta_switch=500;
  
  // for(gt=0;gt<5;gt++){ /* feedback frequency cycle*/
    delta_switch=delta_switch*2;
    printf("oseen%f\t stokes %f\n",OSEENA,STOKES);
    
    d=15;
    mu=0.01;
	strobo =5000;
	
   // for(u=0;u<6;u++){/* distance cycle */
     // mu = mu/(double)10;    
      printf("mu %lf\n",mu);
      dt=2*M_PI*sqrt(6*M_PI*VISC*R)/(10*(double)N_PASSI_CAMPIONAMENTO_PER_CICLO);
	printf("dt %lf\n",dt);
      /* bivariate gaussian */
   
      /* starting positions */
	
	
      x[0]=0 ;
      x[1]=-.85342*SW;
      sigma[0] =1;
      sigma[1] =-1;
      omega=0.001;
      phi[0]=1.6;
      phi[1]=0;
      //alpha=0;
     alpha = -0.01;

	
	for(k=0; k< N_ROW; k++) {
		x[k] = gsl_ran_gaussian(r,XMQ);
		sigma[k] = gsl_ran_gaussian(r,XMQ);
		//x[k] = gsl_rng_uniform(r);
		//sigma[k] = gsl_rng_uniform(r);
		printf(" j %d\t x0  %lf\t sigma0 %lf \n ", k, x[k], sigma[k] );
		//phi[k] = gsl_ran_gaussian(r,XMQ); 
	}

      /* external file names */
	sprintf(fileoutxt,"numerics/%drow_amplitude_d%lf_k%0.2lf_eta%0.1f_sw%0.3lf_eps%0.2lf_R%0.1lf_mu%lf_alpha%lf_v%d.dat",N_ROW, d,K,VISC,SW, eps,R,mu,alpha, N_STEP);    
      xt= fopen (fileoutxt, "w"); 
      
      sprintf(fileouphase,"numerics/%drow_phase_d%lf_k%0.2lf_eta%0.1f_sw%0.3lf_eps%0.2lf_R%0.1lf_mu%lf_alpha%lf_v%d.dat",N_ROW, d,K,VISC,SW, eps,R,mu,alpha, N_STEP);    
	phasext= fopen (fileouphase, "w"); 

	sprintf(fileoutord,"numerics/%drow_ord_par_d%lf_k%0.2lf_eta%0.1f_sw%0.3lf_eps%0.2lf_R%0.1lf_mu%lf_alpha%lf_v%d.dat",N_ROW, d,K,VISC,SW, eps,R,mu,alpha, N_STEP);    
      orderparam= fopen (fileoutord, "w"); 

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
	for(k=0; k< N_ROW; k++) o[k]=0;
	
	for(k=0; k< N_ROW; k++){/* starting k-rowers cycle*/
	  xxx = x[k] ;
	  vvv = v[k];
	  
	  for(l=0; l< N_ROW;l++){/* internal cycle: l!=k */	    
	    if(l != k){ 
	      xx = x[l];
	      vv = v[l];
	      ss = sigma[l];	    	    
	      /* force wih model parameters  F = (r*f -sigma*sw)   */
	      fact  = (ss);
			/* this is the asynchrony*/
	      //fact +=  (1/(double)3)*mu*xx*xx*xx;      
	      /* (k,l) distance */	      
	      dist  = ((xx + d*(double)l) - (xxx + d*(double)k));
	      //dist = d; /* truncated distance */
	      fact = (fact)/fabs(dist*dist);
	      fact *= OSEENA  ;
	      /*add contributes to oseen non diagonal term*/
	      o[k] += (fact); 	    
	      
	    }/* end l-cycle */
	    
	  }/* end if */	
	  
	  if(o[k]*dt > 1)  printf(" oo dt su part %d e' %f, dt o int troppo grandi\n"
				  ,k, o[k]*dt), exit(0); 
	  
	}/* end k-cycle */
	
	
	/*cycle again*/
	ordampl =0;
	ordphase = 0;
	ordampl2 =0;
	ordphase2 = 0;
	for(j=0;j<N_ROW;j++){
	  /*move state, position, oseen in dummy variables*/
	  ss = sigma[j];
	  xx = x[j];
	  oo = o[j];
	  //oo=0;
	  //nnnu = nu[j];
	  /* molecular dynamics on positions*/
	  vv =(ss);
	  //vv = vv*STOKES;
	  vv=vv + oo;
	  deltax = vv*dt;
	  deltas = -xx+ mu*ss*(1-xx*xx) + alpha*xx*xx*xx;
	  deltas= deltas*dt;
	   ss += deltas;
	  deltax = vv*dt; 
	  xx += deltax ;	  
	 
/* 	  /\*TEST*\/ */
/* 	  if(xx > SW || xx < -SW) { */
/* 	    ss= -ss; */
/* /\* 	    if(j==1 && ss==1) fprintf(xt, "%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",0.001*dt*((double)i),x[1]-x[0], x[0], sigma[0], x[1],sigma[1] );  *\/ */
/*  	  } */
	  if(deltax > 0.3) printf("deltax = %g troppo!\n", deltax), exit(0);
	  x[j] = xx;
	  v[j] = vv;
	  sigma[j] =ss;
	  if(isnan(xx) == 1) printf("NANNA!"), exit(0);
	  
	  phase[j]=atan(sigma[j]/x[j]);
	  //phase[j] *= 180/M_PI;
	
	  amplitude[j] = x[j]*x[j] + sigma[j]*sigma[j];
	  amplitude[j] = sqrt(amplitude[j]);
		//i > 40*40000 && 
	    	//if(i > 100*40000 && i%strobo==0){
		  		if(i > 2998*40000) {
			fprintf(xt, "%d\t %lf\t %lf\n", j, 0.001*dt*((double)i), x[j]);
			fprintf(phasext, "%d\t %lf\t %lf\n", j, 0.001*dt*((double)i), phase[j]);
					}

				ordampl = ordampl + amplitude[j];
				ordphase = ordphase + phase[j];
				ordampl2 = ordampl2 + cos(phase[j]);
				ordphase2 = ordphase2 + sin(phase[j]);
	}	/* end j-cycle */ 
	//	if(i > 100*40000 && i%strobo==0) {
		if(i > 2998*40000) {
		  fprintf(xt, "\n");
		  fprintf(phasext, "\n");
		  	  	  }
	 
		  ordampl = ordampl/(double)N_ROW;
		  ordphase = ordphase/(double)N_ROW;
		  ordampl2 = ordampl2/(double)N_ROW;
		  ordphase2 = ordphase2/(double)N_ROW;
		  devphase =0;
		  devamplit =0;

	for(j=0;j<N_ROW;j++){
	  devph[j]= (sin(phase[j])- ordphase2)*(sin(phase[j])- ordphase2) ;
	  devamp[j] = (cos(phase[j])- ordampl2)*(cos(phase[j])- ordampl2) ;
	  devphase += devph[j];
	  devamplit += devamp[j]; 
	  // or 
	  //devph[j]= (phase[j]- ordphase)*(phase[j]- ordphase) ;
	  //devamp[j] = ( amplitude[j]- ordampl )*( amplitude[j]- ordampl) ;
	  //devphase += devph[j];
	  //devamplit += devamp[j]; 

	}
	
   devphase /= (double)N_ROW;
devamplit  /= (double)N_ROW; 
if(i%strobo==0)
   fprintf(orderparam, "%lf\t %lf\t %lf\t %lf\t %lf\t  %lf\t %lf\t %lf\t  %lf\n", 0.001*dt*((double)i), ordampl, ordphase,ordampl2, ordphase2, devamplit, devphase, phase[50], amplitude[50]);

	/* FILE PRINTING */      
	/* xt print */
	/* print according to the frame rate */
	      }/* end time-cycle*/
      
      printf("nswitch %d\n",nswitch);	       
      
      /*free & fclose*/    
      fclose (xt);  
	fclose (phasext);  
  //    } /* uu-cycle */
      //} /* uuu-cycle */
	
  gsl_rng_free(r); /* free rnd generator */
  return (0);
}
