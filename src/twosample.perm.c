
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <stdlib.h>

int compare2(const void *x, const void *y)
{
  double *a, *b;
  a=(double*)x;
  b=(double*)y;
  return ((*a>*b)-(*a<*b));
}


void unpairedperm(int *id, int *nperm, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, int *which1, int *which0, double *s0, double *e, double *f)
{
  double *ex1, *ex0, *ex21, *ex20, *r, *s, *stat, *dstat;
  int i, j, k, *test, *indx;

  ex1=Calloc(*ngene,double);
  ex0=Calloc(*ngene,double);
  ex21=Calloc(*ngene,double);
  ex20=Calloc(*ngene,double);
  r=Calloc(*ngene,double);
  s=Calloc(*ngene,double);
  stat=Calloc(*ngene,double);
  test=Calloc(1,int);
  indx=Calloc((*nperm)*(*ngene),int);
  dstat=Calloc((*nperm)*(*ngene),double);

  for (k=0; k<*nperm; k++){

    for (j=0; j<*ngene; j++){
      ex1[j]=0;
      ex0[j]=0;
      ex21[j]=0;
      ex20[j]=0;
      r[j]=0;
      s[j]=0;
      stat[j]=0;
    }

    /* compute first and second moment to calculate mean and variance */
    for (j=0; j<*ngene; j++){

      for (i=0; i<*nsample; i++){
	test[0] = id[k*(*nsample)+i];

	if (test[0]==0){ex0[j] += matrix[j*(*nsample)+i];}
	if (test[0]==1){ex1[j] += matrix[j*(*nsample)+i];}
      }
      for (i=0; i<*nsample; i++){
	test[0] = id[k*(*nsample)+i];

	if (test[0]==0){ex20[j] += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];}
	if (test[0]==1){ex21[j] += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];}
      }

    }

    for (j=0; j<*ngene; j++){

      ex0[j]=ex0[j]/(*n0);
      ex1[j]=ex1[j]/(*n1);
      ex20[j]=ex20[j]/(*n0);
      ex21[j]=ex21[j]/(*n1);
      
      /* difference in means */
      r[j]=ex1[j] - ex0[j];
      
      /* pooled variance */
      s[j]=(*n1)*(ex21[j]-ex1[j]*ex1[j]) + (*n0)*(ex20[j]-ex0[j]*ex0[j]);
      s[j]=sqrt( s[j]*((double)1/(*n1)+(double)1/(*n0))/((*nsample)-2) );
      

      /* t test statistic */
      if (*meth==1){
	stat[j]=r[j]/s[j];
      }
      
      /* Z test statistic  */
      if (*meth==2){            
	stat[j]=r[j]/(s[j] + *s0);
      }

      /* fold change equivalent */
      if (*meth==3){
	stat[j]=r[j];
      }
      
    }
    

    for (j=0; j<*ngene; j++){
      dstat[j+(*ngene)*k]=fabs(stat[j]);
    }

    qsort((void*)stat,*ngene,sizeof(double),compare2);
    for (j=0; j<*ngene; j++){
      e[j]+=stat[j];
    }

  }


  /* Compute p-values */    
  for (j=0; j<(*nperm)*(*ngene); j++){
    indx[j]=j;
  }

  rsort_with_index((double*)dstat,(int*)indx,(int)(*nperm)*(*ngene));
  
  /* First ngene values correspond to the original labeling */
  for (j=0; j<(*nperm)*(*ngene); j++){      		
    if (indx[j]<*ngene){
      f[indx[j]] = (*nperm)*(*ngene)-j;
    }
  }   

  for (j=0; j<*ngene; j++){ 
    e[j]=e[j]/(*nperm);
    f[j]=f[j]/((*nperm)*(*ngene));
  }
  
  Free(ex1);
  Free(ex0);
  Free(ex21);
  Free(ex20);
  Free(r);
  Free(s);
  Free(stat);
  Free(test);
  Free(dstat);
  Free(indx);
}




void pairedperm(int *id, int *nperm, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, int *which1, int *which0, double *s0, double *e, double *f)
{
  double *r, *s, *ex2, *stat, *dstat;
  double *diff;
  int i, j, k, *indx;

  diff=Calloc(*n1,double);
  r=Calloc(*ngene,double);
  s=Calloc(*ngene,double);
  ex2=Calloc(*ngene,double);
  stat=Calloc(*ngene,double);
  indx=Calloc((*nperm)*(*ngene),int);
  dstat=Calloc((*nperm)*(*ngene),double);


  for (k=0; k<*nperm; k++){

    for (j=0; j<*ngene; j++){
      r[j]=0;
      s[j]=0;
      ex2[j]=0;
      stat[j]=0;
    }
    for (j=0; j<*n1; j++){
      diff[j]=0;
    }

    for (j=0; j<*ngene; j++){

      /* compute pairwise differences and second moment for variance calculation */
      for (i=0; i<*n0; i++){

	diff[i]=matrix[j*(*nsample)+which1[i]]-matrix[j*(*nsample)+which0[i]];
	
	/* especially important for permuted values: where are the original pairs and in which order do we have to substract them? */
	if (id[k*(*nsample)+which0[i]]==1){diff[i]= -diff[i];}
	
	r[j] += diff[i];
	ex2[j] += diff[i]*diff[i];	
      }
    }
      
    for (j=0; j<*ngene; j++){
    
      r[j]=r[j]/(*n1);
      ex2[j]=ex2[j]/(*n1);
      
      s[j]=(*n1)*(ex2[j]-r[j]*r[j]);
      s[j]=sqrt( s[j]/((*n1)*((*n1)-1)) );
    
      /* t test statistic */
      if (*meth==1){
	stat[j]=r[j]/s[j];
      }

      /* Z test statistic */
      if (*meth==2){
	stat[j]=r[j]/(s[j] + *s0);
      }

      /* fold change equivalent */
      if (*meth==3){
	stat[j]=r[j];
      }
      
    }
    

    for (j=0; j<*ngene; j++){
      dstat[j+(*ngene)*k]=fabs(stat[j]);
    }

    qsort((void*)stat,*ngene,sizeof(double),compare2);
    for (j=0; j<*ngene; j++){
      e[j]+=stat[j];
    }
  }


  /* Compute p-values */    
  for (j=0; j<(*nperm)*(*ngene); j++){
    indx[j]=j;
  }

  rsort_with_index((double*)dstat,(int*)indx,(int)(*nperm)*(*ngene));
  
  /* First ngene values correspond to the original labeling */
  for (j=0; j<(*nperm)*(*ngene); j++){      		
    if (indx[j]<*ngene){
      f[indx[j]] = (*nperm)*(*ngene)-j;
    }
  }   
  
  for (j=0; j<*ngene; j++){ 
    e[j]=e[j]/(*nperm);
    f[j]=f[j]/((*nperm)*(*ngene));
  }

  Free(diff);
  Free(r);
  Free(s);
  Free(ex2);
  Free(stat);
  Free(dstat);
  Free(indx);
}

