#include <stdio.h>
#include "jacobi.h"

void jacobistepvort(double *zetnew, double *psinew,
  double *zet, double *psi,
  int m, int n, double re) {

  for(int i=1;i<m+1;i++) {
    for(int j=1;j<n+1;j++)	{
     psinew[i*(n+2)+j]=0.25*(  psi[(i-1)*(n+2)+j]+psi[(i+1)*(n+2)+j]+psi[i*(n+2)+j-1]+psi[i*(n+2)+j+1]
      - zet[i*(n+2)+j] );
   }
 }

 for(int i=1;i<m+1;i++) {
  for(int j=1;j<n+1;j++)	{
   zetnew[i*(n+2)+j]=0.25*(zet[(i-1)*(n+2)+j]+zet[(i+1)*(n+2)+j]+zet[i*(n+2)+j-1]+zet[i*(n+2)+j+1])
   - re/16.0*(
     (  psi[i*(n+2)+j+1]-psi[i*(n+2)+j-1])*(zet[(i+1)*(n+2)+j]-zet[(i-1)*(n+2)+j])
     - (psi[(i+1)*(n+2)+j]-psi[(i-1)*(n+2)+j])*(zet[i*(n+2)+j+1]-zet[i*(n+2)+j-1])
     );
 }
}
}

double deltasq(double *newarr, double *oldarr, int m, int n)
{
  int i, j;

  double dsq=0.0;
  double tmp;

  for(int i=1;i<m+1;i++) {
    for(int j=1;j<n+1;j++)	{
     tmp = newarr[i*(n+2)+j]-oldarr[i*(n+2)+j];
     dsq += tmp*tmp;
   }
 }

 return dsq;
}
