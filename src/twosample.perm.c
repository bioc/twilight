
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>

int compare2(const void *x, const void *y)
{
  double *a, *b;
  a=(double*)x;
  b=(double*)y;
  return ((*a>*b)-(*a<*b));
}


void unpairedperm(int *id, int *nperm, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, double *sobs, int *which1, int *which0, double *s0, double *e, double *f)
{
  double *ex1, *ex0, *ex21, *ex20, *r, *s, *ssort, *stat;
  int i, j, k, *test;
    
  if ((ex1=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex0=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex21=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex20=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((r=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((s=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ssort=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((stat=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((test=calloc(1,sizeof(int)))==0) {printf("Error, could not allocate memory");}


  for (k=0; k<*nperm; k++){

    for (j=0; j<*ngene; j++){
      ex1[j]=0;
      ex0[j]=0;
      ex21[j]=0;
      ex20[j]=0;
      r[j]=0;
      s[j]=0;
      ssort[j]=0;
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
      
      /* fold change equivalent */
      if (*meth==3){
	stat[j]=r[j];
      }
      
      ssort[j]=s[j];
    }
    

    /* Z test statistic, needs calculation of median(s) */
    if (*meth==2){
            
      if (*s0==0){
	qsort((void*)ssort,*ngene,sizeof(double),compare2);
	
	if (fmod(*ngene,2)==1){
	  *s0=ssort[(*ngene-1)/2];
	}
	if (fmod(*ngene,2)==0){
	  *s0=(ssort[(*ngene)/2]+ssort[(*ngene)/2-1])/2;
	}
      }

      for (j=0; j<*ngene; j++){
	stat[j]=r[j]/(s[j] + *s0);
      }
    }

    /* Compute p-values by comparing stat and sobs */
    for (j=0; j<*ngene; j++){
      if ( (fabs(stat[j])>=fabs(sobs[j])) | (fabs(sobs[j])-fabs(stat[j])<10e-10) ){
	f[j]+=1;
      }
    }

    qsort((void*)stat,*ngene,sizeof(double),compare2);
    for (j=0; j<*ngene; j++){
      e[j]+=stat[j];
    }

  }

  for (j=0; j<*ngene; j++){ 
    e[j]=e[j]/(*nperm);
    f[j]=f[j]/(*nperm); 
  }
  
  free(ex1);
  free(ex0);
  free(ex21);
  free(ex20);
  free(r);
  free(s);
  free(ssort);
  free(stat);
  free(test);
}




void pairedperm(int *id, int *nperm, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, double *sobs, int *which1, int *which0, double *s0, double *e, double *f)
{
  double *r, *s, *ssort, *ex2, *stat;
  double *diff;
  int i, j, k;

  if ((diff=calloc(*n1,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((r=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((s=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ssort=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((ex2=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}
  if ((stat=calloc(*ngene,sizeof(double)))==0) {printf("Error, could not allocate memory");}

  for (k=0; k<*nperm; k++){

    for (j=0; j<*ngene; j++){
      r[j]=0;
      s[j]=0;
      ssort[j]=0;
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

      /* fold change equivalent */
      if (*meth==3){
	stat[j]=r[j];
      }
      
      ssort[j]=s[j];
    }
    

    /* Z test statistic, needs calculation of median(s) */
    if (*meth==2){
      
      if (*s0==0){
	qsort((void*)ssort,*ngene,sizeof(double),compare2);
	
	if (fmod(*ngene,2)==1){
	  *s0=ssort[(*ngene-1)/2];
	}
	if (fmod(*ngene,2)==0){
	  *s0=(ssort[(*ngene)/2]+ssort[(*ngene)/2-1])/2;
	}
      }

      for (j=0; j<*ngene; j++){
	stat[j]=r[j]/(s[j] + *s0);
      }
    }


    /* Compute p-values by comparing stat and sobs */
    for (j=0; j<*ngene; j++){
      if ( (fabs(stat[j])>=fabs(sobs[j])) | (fabs(sobs[j])-fabs(stat[j])<10e-10) ){
	f[j]+=1;
      }
    }
    
    qsort((void*)stat,*ngene,sizeof(double),compare2);
    for (j=0; j<*ngene; j++){
      e[j]+=stat[j];
    }
    
  }

  for (j=0; j<*ngene; j++){ 
    e[j]=e[j]/(*nperm);
    f[j]=f[j]/(*nperm); 
  }

  free(diff);
  free(r);
  free(s);
  free(ssort);
  free(ex2);
  free(stat);
}

