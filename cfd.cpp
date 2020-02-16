#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "boundary.h"
#include "jacobi.h"
#include "cfdio.h"

int main(int argc, char **argv)
{
  int printfreq=1000; //output frequency
  double error, bnorm;
  double tolerance=0.0001; //tolerance for convergence. <=0 means do not check
  int test_mode = 0;

  //main arrays
  double *psi, *zet;
  //temporary versions of main arrays
  double *psitmp, *zettmp;

  //command line arguments
  int scalefactor, numiter;

  double re = -1.0; // Reynold's number - must be less than 3.7

  //simulation sizes
  int bbase=10;
  int hbase=15;
  int wbase=5;
  int mbase=32;
  int nbase=32;

  int checkerr = 0;

  int m,n,b,h,w;
  int iter;
  int i,j;

  double tstart, tstop, ttot, titer;

  //do we stop because of tolerance?
  if (tolerance > 0) {checkerr=1;}

  //check command line parameters and parse them

  if (argc <3|| argc >4)
  {
    printf("Running in TEST mode\n");
    scalefactor = 10;
    numiter = 5000;
    re = -1.0;
    test_mode = 1;
  }
  else
  {
    scalefactor=atoi(argv[1]);
    numiter=atoi(argv[2]);
  }

  if (argc == 4)
  {
    re=atof(argv[3]);
  }


  if(!checkerr)
  {
    printf("Scale Factor = %i, iterations = %i\n",scalefactor, numiter);
  }
  else
  {
    printf("Scale Factor = %i, iterations = %i, tolerance= %g\n",scalefactor,numiter,tolerance);
  }


  printf("Reynolds number = %f\n",re);

  //Calculate b, h & w and m & n
  b = bbase*scalefactor;
  h = hbase*scalefactor;
  w = wbase*scalefactor;
  m = mbase*scalefactor;
  n = nbase*scalefactor;

  re = re / (double)scalefactor;

  printf("Running CFD on %d x %d grid in serial\n",m,n);

  //allocate arrays

  psi    = new double[(m+2)*(n+2)];
  psitmp = new double[(m+2)*(n+2)];

  //zero the psi array
  for (i=0;i<m+2;i++)
  {
    for(j=0;j<n+2;j++)
    {
      psi[i*(n+2)+j]=0.0;
    }
  }

  //allocate arrays
  zet =   new double[(m+2)*(n+2)];;
  zettmp = new double[(m+2)*(n+2)];;

  //zero the zeta array

  for (i=0;i<m+2;i++)
  {
    for(j=0;j<n+2;j++)
    {
      zet[i*(n+2)+j]=0.0;
    }
  }

  //set the psi boundary conditions

  boundarypsi(psi,m,n,b,h,w);

  //compute normalisation factor for error

  bnorm=0.0;

  for (i=0;i<m+2;i++)
  {
    for (j=0;j<n+2;j++)
    {
      bnorm += psi[i*(n+2)+j]*psi[i*(n+2)+j];
    }
  }


  //update zeta BCs that depend on psi
  boundaryzet(zet,psi,m,n);

  //update normalisation

  for (i=0;i<m+2;i++)
  {
    for (j=0;j<n+2;j++)
    {
      bnorm += zet[i*(n+2)+j]*zet[i*(n+2)+j];
    }
  }


  bnorm=sqrt(bnorm);

  //begin iterative Jacobi loop

  printf("\nStarting main loop...\n\n");

  tstart=gettime();

  for(iter=1;iter<=numiter;iter++)
  {

    //calculate psi for next iteration

    jacobistepvort(zettmp,psitmp,zet,psi,m,n,re);

    //calculate current error if required

    if (checkerr || iter == numiter)
    {
      error = deltasq(psitmp,psi,m,n);

      error += deltasq(zettmp,zet,m,n);

      error=sqrt(error);
      error=error/bnorm;
    }

    //quit early if we have reached required tolerance

    if (checkerr)
    {
      if (error < tolerance)
      {
        printf("Converged on iteration %d\n",iter);
        break;
      }
    }

    //copy back

    for(i=1;i<=m;i++)
    {
      for(j=1;j<=n;j++)
      {
        psi[i*(n+2)+j]=psitmp[i*(n+2)+j];
      }
    }


    for(i=1;i<=m;i++)
    {
      for(j=1;j<=n;j++)
      {
        zet[i*(n+2)+j]=zettmp[i*(n+2)+j];
      }
    }

  	//update zeta BCs that depend on psi
    boundaryzet(zet,psi,m,n);

    //print loop information

    if(iter%printfreq == 0)
    {
      if (!checkerr)
      {
        printf("Completed iteration %d\n",iter);
      }
      else
      {
        printf("Completed iteration %d, error = %g\n",iter,error);
      }
    }
  }

  if (iter > numiter) iter=numiter;

  tstop=gettime();

  ttot=tstop-tstart;
  titer=ttot/(double)iter;


  //print out some stats

  printf("\n... finished\n");
  printf("After %d iterations, the error is %g\n",iter,error);
  printf("Time for %d iterations was %g seconds\n",iter,ttot);
  printf("Each iteration took %g seconds\n",titer);

  if (test_mode == 1)
  {
    if (iter == 5000 && error > 0.0002 && error < 0.0003)
    {
      printf("TEST PASSED\n");
    }
    else
    {
      printf("TEST FAILED\n");
    }
  }

    //output results

  writedatafiles(psi,m,n, scalefactor);

  writeplotfile(m,n,scalefactor);

  //free un-needed arrays
  delete[] psi;
  delete[] psitmp;

  delete[] zet;
  delete[] zettmp;

  printf("... finished\n");

  return 0;
}
