
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>

int compare7(const void *x, const void *y)
{
  double *a, *b;
  a=(double*)x;
  b=(double*)y;
  return ((*a>*b)-(*a<*b));
}

/* Sort in decreasing order */
int comparedecr7(const void *x, const void *y)
{
  double *a, *b;
  a=(double*)x;
  b=(double*)y;
  return ((*a<*b)-(*a>*b));
}


void corsingle(double *vector, double *matrix, int *ngene, int *nsample, double *e)
{
  double *ex0, *ex1, *ex20, *ex21, *exboth;
  int i, j;

  if ((ex0=calloc(1,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex1=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex20=calloc(1,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex21=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((exboth=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}

  /* compute first and second moments */
  for (i=0; i<*nsample; i++){
    ex0[0]  += vector[i];
    ex20[0] += vector[i]*vector[i];
  }

  for (j=0; j<*ngene; j++){
    for (i=0; i<*nsample; i++){
      ex1[j]    += matrix[j*(*nsample)+i];
      ex21[j]   += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];
      exboth[j] += matrix[j*(*nsample)+i]*vector[i];
    }
    /* correlation coefficient */
    e[j]=(exboth[j]-ex0[0]*ex1[j]/(*nsample))/sqrt((ex20[0]-ex0[0]*ex0[0]/(*nsample))*(ex21[j]-ex1[j]*ex1[j]/(*nsample)));
  }
  
  free(ex0);
  free(ex1);
  free(ex20);
  free(ex21);
  free(exboth);
}



void corperm(double *vecperm, int *nperm, double *matrix, int *ngene, int *nsample, double *sobs, double *e, double *f)
{
  double *ex0, *ex1, *ex20, *ex21, *exboth, *stat;
  int i, j, k, j0, j1, j2, jj;

  if ((ex0=calloc(1,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex1=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex20=calloc(1,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex21=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((exboth=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((stat=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}

  for (k=0; k<*nperm; k++){

    ex0[0] = 0;
    ex20[0] = 0;

    for (j=0; j<*ngene; j++){
      ex1[j]=0;
      ex21[j]=0;
      exboth[j]=0;
      stat[j]=0;
    }

    /* compute first and second moments */
    for (i=0; i<*nsample; i++){
      ex0[0]  += vecperm[k*(*nsample)+i];
      ex20[0] += vecperm[k*(*nsample)+i]*vecperm[k*(*nsample)+i];
    }
   
    /* compute mixed term */
    for (j=0; j<*ngene; j++){
      for (i=0; i<*nsample; i++){
	ex1[j]    += matrix[j*(*nsample)+i];
	ex21[j]   += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];
	exboth[j] += matrix[j*(*nsample)+i]*vecperm[k*(*nsample)+i];
      }

      /* correlation coefficient */
      stat[j]=(exboth[j]-ex0[0]*ex1[j]/(*nsample))/sqrt((ex20[0]-ex0[0]*ex0[0]/(*nsample))*(ex21[j]-ex1[j]*ex1[j]/(*nsample)));
    }
    
    qsort((void*)stat,*ngene,sizeof(double),compare7);
    for (j=0; j<*ngene; j++){
      e[j]+=stat[j];
    }

    for (j=0; j<*ngene; j++){
      stat[j]=fabs(stat[j]);
    }
    qsort((void*)stat,*ngene,sizeof(double),comparedecr7);

    /* Compute p-values by comparing stat and sobs */
    /* Find rank of sobs in stat by dividing intervals into halfs (bin search) */
    for (j=0; j<*ngene; j++){
      jj=0;
      j0=0;
      j2=(*ngene)-1;
      j1=j0+ceil((j2-j0)/2);

        while ((j0!=j1)&&(j1!=j2)){
	  if (fabs(sobs[j])<stat[j1]){
	    j0=j1;
	    j1=j0+ceil((j2-j0)/2);
	  }
	  if (fabs(sobs[j])>stat[j1]){
	    j2=j1;
	    j1=j0+floor((j2-j0)/2);
	  }
	  if (fabs(sobs[j])==stat[j1]){
	    j2=j1;
	  }
	}
	if (fabs(sobs[j])>stat[j1]){
	  jj=j0;
	}
	if (fabs(sobs[j])<stat[j1]){
	  jj=j2;
	}
	if (fabs(sobs[j])==stat[j1]){	  
	  jj=j1+1;
	}
	if (fabs(sobs[j])==stat[(*ngene)-1]){	  
	  jj=(*ngene);
	}
		
	f[j]+=jj;
    }   
  }

  for (j=0; j<*ngene; j++){ 
    e[j]=e[j]/(*nperm);
    f[j]=f[j]/((*nperm)*(*ngene));
  }
    
  free(ex0);
  free(ex1);
  free(ex20);
  free(ex21);
  free(exboth);
  free(stat);
}



void corci(double *vecperm, int *nperm, double *matrix, int *ngene, int *nsample, double *sobs, double *e)
{
  double *ex0, *ex1, *ex20, *ex21, *exboth, *stat;
  int i, j, k;

  if ((ex0=calloc(1,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex1=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex20=calloc(1,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex21=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((exboth=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((stat=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}

  for (k=0; k<*nperm; k++){

    ex0[0] = 0;
    ex20[0] = 0;

    for (j=0; j<*ngene; j++){
      ex1[j]=0;
      ex21[j]=0;
      exboth[j]=0;
      stat[j]=0;
    }

    /* compute first and second moments */
    for (i=0; i<*nsample; i++){
      ex0[0]  += vecperm[k*(*nsample)+i];
      ex20[0] += vecperm[k*(*nsample)+i]*vecperm[k*(*nsample)+i];
    }
   
    /* compute mixed term */
    for (j=0; j<*ngene; j++){
      for (i=0; i<*nsample; i++){
	ex1[j]    += matrix[j*(*nsample)+i];
	ex21[j]   += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];
	exboth[j] += matrix[j*(*nsample)+i]*vecperm[k*(*nsample)+i];
      }

      /* correlation coefficient */
      stat[j]=(exboth[j]-ex0[0]*ex1[j]/(*nsample))/sqrt((ex20[0]-ex0[0]*ex0[0]/(*nsample))*(ex21[j]-ex1[j]*ex1[j]/(*nsample)));
    }
        
    /* Maximum absolute difference between stat and sobs */
    qsort((void*)stat,*ngene,sizeof(double),compare7);
    qsort((void*)sobs,*ngene,sizeof(double),compare7);

    for (j=0; j<*ngene; j++){
      stat[j]-=sobs[j];
      stat[j]=fabs(stat[j]);
    }

    qsort((void*)stat,*ngene,sizeof(double),comparedecr7);

    e[k]=stat[0];     
  }

    
  free(ex0);
  free(ex1);
  free(ex20);
  free(ex21);
  free(exboth);
  free(stat);
}


