
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>


/* Returns maximum absolute difference of empirical and uniform distribution */
double empirical(double *a, int na)
{

  int i;
  int j=na;
  double *b;
  double e;

  if ((b=malloc(sizeof(double)*1))==0) {printf("Error, could not allocate memory");}

  e=fabs(a[na-1]-1);

  for (i=na-2; i>0; i--){

    if (a[i]!=a[i+1]){
      j=i+1;
      b[0]=fabs( a[i] - ( (double)j / ((double)na) ) );

      if (b[0]>e){e=b[0];}
    }
  }

  free(b);

  return e;
}


/* SEP */
void sep(double *xin, int *nxin, double *lambda, int *xout, double *funout)
{

  int *ix;
  int i;
  int j=0;
  int k;
  int count=0;
  int nyin=0;
  int randnum;
  double *yin;
  double *objfunc;

  if ((ix=malloc(sizeof(int)*(*nxin)))==0) {printf("Error, could not allocate memory");}
  if ((objfunc=malloc(sizeof(double)*2))==0) {printf("Error, could not allocate memory");}

  for (i=0; i<*nxin; i++){ ix[i]=1; }

  objfunc[0]=empirical(xin,*nxin);

  for (k=1; k<10000;)
    {
      randnum=(int)( ((double) *nxin)*rand()/(RAND_MAX+1.0) );
      
      ix[randnum]=fabs(1-ix[randnum]);

      for (i=0; i<*nxin; i++){ 
	if (ix[i]==1){nyin=nyin+1;} 
      }
      
      if ((yin=malloc(sizeof(double)*nyin))==0) {printf("Error, could not allocate memory");}
      
      for (i=0; i<*nxin; i++){ if (ix[i]==1){yin[j]=xin[i]; j=j+1;} }

      objfunc[1]=empirical(yin,nyin) + (*lambda)*(*nxin - (double)nyin)*log(*nxin - (double)nyin)/(*nxin);
            
      if (objfunc[1]<objfunc[0])
	{
	  objfunc[0]=objfunc[1];
	  count=0;
	  k=k+1;

	  funout[0]=empirical(yin,nyin);
	}
      else
	{
	  ix[randnum]=fabs(1-ix[randnum]);
	  count=count+1;
	}
      
      nyin=0;
      j=0;
      free(yin);
      if (count>2*(*nxin)){ break; }
    }


  for (i=0; i<*nxin; i++)
    {
      xout[i]=ix[i];
    }

  free(ix);
  free(objfunc);
 
}
