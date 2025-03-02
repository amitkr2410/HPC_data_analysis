/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This file is for lattice generation with the RHMC algorithm */

#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"


#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#endif

#ifdef HAVE_QPHIX
#include "../include/generic_qphix.h"
#endif

#ifdef MILC_GLOBAL_DEBUG
#include "debug.h"
#endif /* MILC_GLOBAL_DEBUG */

/* For information */
#define NULL_FP -1

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int MyFunction( int argc, char **argv )
{
  int i,MeasurementCount,traj_done, naik_index;
  int prompt;
  int s_iters=0, iters=0;
  double dtime, dclock();

  //Plaquette and Field-strength variable
  double SS_Plaq=0.0, ST_Plaq=0.0;
  double Current_Plaq=0.0, Sum_Plaq=0.0, Average_Plaq=0.0;
  double SumPBP=0.0, AveragePBP=0.0;
  complex CurrentPolyakovLoop, SumPolyakovLoop, AveragePolyakovLoop;
  complex CurrentTraceF3iF3iMinusF4iF4i, SumTraceF3iF3iMinusF4iF4i, AverageTraceF3iF3iMinusF4iF4i;
  complex CurrentTraceF4iF3iPlusF3iF4i, SumTraceF4iF3iPlusF3iF4i, AverageTraceF4iF3iPlusF3iF4i;
  //Initialize variable to zero
  CurrentPolyakovLoop=cmplx(0.0,0.0); SumPolyakovLoop=cmplx(0.0,0.0); AveragePolyakovLoop=cmplx(0.0,0.0);
  CurrentTraceF3iF3iMinusF4iF4i =cmplx(0.0,0.0); CurrentTraceF4iF3iPlusF3iF4i =cmplx(0.0,0.0);
  SumTraceF3iF3iMinusF4iF4i =cmplx(0.0,0.0); AverageTraceF3iF3iMinusF4iF4i =cmplx(0.0,0.0);
  SumTraceF4iF3iPlusF3iF4i  =cmplx(0.0,0.0); AverageTraceF4iF3iPlusF3iF4i =cmplx(0.0,0.0);

  //FileName to save observables
  FILE *fplaquette, *ftracefmunu;
  char FileNamePlaq[1000], FileNameFmunu[1000];

  // Initialization 
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1); 
  g_sync();
  
  /* set up lattice parameters */
  prompt = setup();

  printf("Amit MyFFApp/control.c prompt= %d \n",prompt);
  printf("Amit MyFFApp/control.c before while(readin(prompt)==0) called \n");
  /* loop over input sets */
  while( readin(prompt) == 0)
    {
      sprintf(FileNamePlaq,"OutPut/DataPlaquette%d.txt", (int)(beta*1000));
      sprintf(FileNameFmunu,"OutPut/DataTraceFmunu%d.txt",(int)(beta*1000));
      fplaquette = fopen(FileNamePlaq,"w");
      ftracefmunu = fopen(FileNameFmunu,"w");

      fprintf(fplaquette,"#Beta=%e, Nt=%d, Ns=%d^3 \n", beta, nt, nx);
      fprintf(fplaquette,"#Iters \t Current_Plaq \t Average_Plaq \t CurrentPolyakovLoop.real \t CurrentPolyakovLoop.imag \t AveragePolyakovLoop.real \t AveragePolyakovLoop.imag \t PBP \t AveragePBP\n");
      fprintf(ftracefmunu,"#Beta=%e, Nt=%d, Ns=%d^3 \n", beta, nt, nx);
      fprintf(ftracefmunu,"#Iter \t TraceF3iF3iMinusF4iF4i.real \t TraceF3iF3iMinusF4iF4i.imag \t AvgTraceF3iF3iMinusF4iF4i.real \t AvgTraceF3iF3iMinusF4iF4i.imag \t TraceF4iF3iPlusF3iF4i.real \t TraceF4iF3iPlusF3iF4i.imag \t AvgTraceF4iF3iPlusF3iF4i.real \t AvgTraceF4iF3iPlusF3iF4i.imag \n");
      /* perform warmup trajectories */
      #ifdef MILC_GLOBAL_DEBUG
      global_current_time_step = 0;
      #endif /* MILC_GLOBAL_DEBUG */
    
      // Start application timer
      dtime = -dclock();
      printf(" Amit MyFFApp/control.c inside while(readin(prompt)==0) \n");
      for( traj_done=0; traj_done < warms; traj_done++ )
      	{
	  rephase(OFF);
	  SS_Plaq=0.0; ST_Plaq=0.0;
	  d_plaquette(&SS_Plaq, &ST_Plaq);
	  printf("Amit MyFFApp/control.c Plaquette = (%e,%e)\n",SS_Plaq, ST_Plaq);	  
	  rephase(ON);
      	  update();	  
        }
      
      node0_printf("default MyFFApp/control.c WARMUPS COMPLETED\n"); fflush(stdout);
      
      /* perform measuring trajectories, reunitarizing and measuring 	*/
      MeasurementCount=0;		/* number of measurements 		*/
      
      
      for( traj_done=0; traj_done < trajecs; traj_done++ )
	{ 
          #ifdef MILC_GLOBAL_DEBUG
          #ifdef HISQ_REUNITARIZATION_DEBUG
	  {
	    int isite, idir;
	    site *s;
	    FORALLSITES(isite,s)
	      {
		for( idir=XUP;idir<=TUP;idir++ )
		  {
		    lattice[isite].on_step_Y[idir] = 0;
		    lattice[isite].on_step_W[idir] = 0;
		    lattice[isite].on_step_V[idir] = 0;
		  }
	      }
	  }
          #endif /* HISQ_REUNITARIZATION_DEBUG */
          #endif /* MILC_GLOBAL_DEBUG */
	  
	  /* do the Measurement and then the trajectories */

	  rephase(OFF);
          printf(" Amit MyFFApp/control.c s_iters=update() called at iters = %d \n", iters);
	  /* measure every "propinterval" trajectories */
	  //	 if( (traj_done%propinterval)==(propinterval-1) )
	  //{
	      iters = 1 + iters;
	      MeasurementCount = MeasurementCount + 1;
	      if(iters ==200) 
		{
		  MeasurementCount = 1; Sum_Plaq=0.0;
		  SumPolyakovLoop = cmplx(0.0,0.0);
		  SumTraceF3iF3iMinusF4iF4i = cmplx(0.0,0.0);  SumTraceF4iF3iPlusF3iF4i = cmplx(0.0,0.0);
		  SumPBP = 0.0;
		}
	      /* call gauge_variable fermion_variable measuring routines */
	      //rephase(OFF);	      
	      /* Compute plaquette and display output */
	      SS_Plaq=0.0; ST_Plaq=0.0;
	      d_plaquette(&SS_Plaq, &ST_Plaq);
	      Current_Plaq = ((SS_Plaq+ST_Plaq)/2.0);
	      Sum_Plaq = Sum_Plaq + Current_Plaq;
	      Average_Plaq = Sum_Plaq/MeasurementCount;
	      printf("Amit MyFFApp/control.c Plaquette=(%e,%e), AveragePlaq=%e \n",SS_Plaq, ST_Plaq, Average_Plaq);

	      /* Calculate trace of polyakov loop */
	      CurrentPolyakovLoop=cmplx(0.0,0.0); 
	      CurrentPolyakovLoop = ploop();
	      CADD(SumPolyakovLoop, CurrentPolyakovLoop, SumPolyakovLoop);
	      CDIVREAL(SumPolyakovLoop, MeasurementCount, AveragePolyakovLoop);
	      printf("Amit MyFFApp/control.c PLoop=(%e,%e), AvgPLoop=(%e,%e)\n", CurrentPolyakovLoop.real, CurrentPolyakovLoop.imag, AveragePolyakovLoop.real, AveragePolyakovLoop.imag);

	      /* Calculate trace of fmunu and output */
	      
	      CurrentTraceF3iF3iMinusF4iF4i = cmplx(0.0,0.0);  CurrentTraceF4iF3iPlusF3iF4i = cmplx(0.0,0.0);  	      
	      fmunu_fmunu(&CurrentTraceF3iF3iMinusF4iF4i, &CurrentTraceF4iF3iPlusF3iF4i);
	      CADD(SumTraceF3iF3iMinusF4iF4i, CurrentTraceF3iF3iMinusF4iF4i, SumTraceF3iF3iMinusF4iF4i);
	      CADD(SumTraceF4iF3iPlusF3iF4i, CurrentTraceF4iF3iPlusF3iF4i, SumTraceF4iF3iPlusF3iF4i);
	      CDIVREAL(SumTraceF3iF3iMinusF4iF4i, MeasurementCount, AverageTraceF3iF3iMinusF4iF4i);
	      CDIVREAL(SumTraceF4iF3iPlusF3iF4i, MeasurementCount, AverageTraceF4iF3iPlusF3iF4i);	      
	      printf("Amit MyFFApp/control.c TraceF3iF3iMinusF4iF4i=(%e,%e), AvgTrace=(%e,%e) \n",CurrentTraceF3iF3iMinusF4iF4i.real, CurrentTraceF3iF3iMinusF4iF4i.imag, AverageTraceF3iF3iMinusF4iF4i.real, AverageTraceF3iF3iMinusF4iF4i.imag);
	      printf("Amit MyFFApp/control.c TraceF4iF3iPlusF3iF4i=(%e,%e), AvgTrace=(%e,%e) \n",CurrentTraceF4iF3iPlusF3iF4i.real, CurrentTraceF4iF3iPlusF3iF4i.imag, AverageTraceF4iF3iPlusF3iF4i.real, AverageTraceF4iF3iPlusF3iF4i.imag); 
	      //fflush(stdout);
	      
	      fprintf(fplaquette,"%d \t %e \t %e \t %e \t %e \t %e \t %e ", iters, Current_Plaq, Average_Plaq, CurrentPolyakovLoop.real, CurrentPolyakovLoop.imag, AveragePolyakovLoop.real, AveragePolyakovLoop.imag);
	      fprintf(ftracefmunu,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", iters, CurrentTraceF3iF3iMinusF4iF4i.real, CurrentTraceF3iF3iMinusF4iF4i.imag, AverageTraceF3iF3iMinusF4iF4i.real, AverageTraceF3iF3iMinusF4iF4i.imag, CurrentTraceF4iF3iPlusF3iF4i.real, CurrentTraceF4iF3iPlusF3iF4i.imag, AverageTraceF4iF3iPlusF3iF4i.real, AverageTraceF4iF3iPlusF3iF4i.imag );
	      
	 rephase(ON);
	 /* Compute chiral condensate pbp, etc */
	 /* Make fermion links if not already done */
	 restore_fermion_links_from_site(fn_links, par_buf.prec_pbp);
	 for(i = 0; i < par_buf.num_pbp_masses; i++)
	   {
             #if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
	     naik_index = par_buf.ksp_pbp[i].naik_term_epsilon_index;
             #else
	     naik_index = 0;
             #endif
	     f_meas_imp_field( par_buf.npbp_reps, &par_buf.qic_pbp[i], par_buf.ksp_pbp[i].mass, naik_index, fn_links);
	   }
	 SumPBP = SumPBP + PBP;
	 AveragePBP = SumPBP/MeasurementCount;
	 printf("Amit MyFFApp/control.c PBP=%e, AveragePBP=%e \n", PBP, AveragePBP);
	 fprintf(fplaquette,"%e \t %e \n",PBP, AveragePBP);
	 //fflush(stdout);
	 
	 if(traj_done < trajecs - 1)
	   {
	     s_iters=update();
	   }
	}	   	/* end loop over trajectories */
       
      if(this_node==1)
	{
	     printf("\n\n default MyFFApp/control.c RUNNING COMPLETED, This node is %d\n\n",this_node); 
	     fflush(stdout);	     
	}

      dtime += dclock();
      if(this_node==1)
	{
	  printf("\n\n Default MyFFApp/control.c Time = %e seconds\n",dtime);
	  printf("Default MyFFApp/control.c total_iters = %d \n\n",iters);
	}
      fflush(stdout);
      
      /* save lattice if requested */
      if( saveflag != FORGET )
	{
	  rephase( OFF );
	  save_lattice( saveflag, savefile, stringLFN );
	  rephase( ON );
	}
      
      /* Destroy fermion links (created in readin() */
      
#if FERM_ACTION == HISQ
      destroy_fermion_links_hisq(fn_links);
#elif FERM_ACTION == HYPISQ
      destroy_fermion_links_hypisq(fn_links);
#else
      destroy_fermion_links(fn_links);
#endif
      fn_links = NULL;
	   }
  
  fclose(fplaquette);
  fclose(ftracefmunu);
  normal_exit(0);
  return 0;
}

