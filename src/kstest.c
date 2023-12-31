
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


int compare11(const void *x, const void *y)
{
  double *a, *b;
  a=(double*)x;
  b=(double*)y;
  return ((*a>*b)-(*a<*b));
}


/* returns indeces of *a in "first" mode */
int* indexx(int n, double *a)
{
  int *indx;
  int i;

  indx=Calloc(n, int);

  for (i=0; i<n; i++){
    a[i]=fabs(a[i]);
    indx[i]=i;
  }
  
  rsort_with_index((double*)a,(int*)indx,(int)n);

  return indx;
  /* indx will be freed in kolmogoroff with free(indx) */
}


/* Returns maximum absolute difference of empirical and uniform distribution */
/* Kolmogoroff-Smirnoff-Test of uniformity */
/* Caution: Input *a are the absolute observed values, not p-values! */
double* kolmogoroff(double *a, int nrowa, int ncola)
{

  int i;
  int k;
  double *b;
  double *asort;
  double *pval;
  int n=nrowa*ncola;
  int *indx;
  double *eout;


  b=Calloc(1, double);
  asort=Calloc(nrowa, double);
  pval=Calloc(n, double);
  eout=Calloc(ncola, double);

  indx=indexx(n,a);

  for (i=0; i<n; i++){
    pval[indx[i]]=(double)(n-i)/((double)n);   /* returns pooled p-value */
  }

  for (k=0; k<ncola; k++){
    
    for (i=0; i<nrowa; i++){
      asort[i]=pval[i*ncola+k];
    }
    
    qsort((void*)asort,nrowa,sizeof(double),compare11);    

    eout[k]=0;
    
    for (i=1; i<nrowa; i++){
      
      if (asort[i]!=asort[i-1]){
	b[0]=fabs( asort[i] - ( (double) i/ ((double) nrowa) ) );
	if (b[0]>eout[k]){eout[k]=b[0];}
	b[0]=fabs( asort[i] - ( ((double) i+1)/ ((double) nrowa) ) );
	if (b[0]>eout[k]){eout[k]=b[0];}
      }
    }    
  }

  Free(b); 
  Free(asort);  
  Free(pval);  
  Free(indx);  

  return eout;

  /* eout will be freed in unpairedKSTEST with free(e) */
}



/* Computes test statistics, pooled p-values and KS statistics */
void unpairedKSTEST(int *id, int *nperm, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, double *s0, double *f)
{
  double *ex1, *ex0, *ex21, *ex20, *r, *s, *ssort, *stat, *statall, *e;
  int i, j, k, *test;
  int nrowa=*ngene, ncola=*nperm;

  ex1=Calloc(*ngene,double);
  ex0=Calloc(*ngene,double);
  ex21=Calloc(*ngene,double);
  ex20=Calloc(*ngene,double);
  r=Calloc(*ngene,double);
  s=Calloc(*ngene,double);
  ssort=Calloc(*ngene,double);
  stat=Calloc(*ngene,double);
  statall=Calloc((*ngene)*(*nperm),double);
  test=Calloc(1,int);


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
	qsort((void*)ssort,*ngene,sizeof(double),compare11);

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


    /* puts the ngene*nperm matrix of scores into a vector to be passed over to kolmogoroff */
    for (j=0; j<*ngene; j++){
      statall[j*(*nperm)+k]=fabs(stat[j]);
    }
  }

  e=kolmogoroff(statall,nrowa,ncola);

  for (k=0; k<*nperm; k++){
    f[k]=e[k];
  }

  Free(ex1);
  Free(ex0);
  Free(ex21);
  Free(ex20);
  Free(r);
  Free(s);
  Free(ssort);
  Free(stat);
  Free(statall);
  Free(test);
  Free(e);
}



/* Computes test statistics, pooled p-values and KS statistics */
void pairedKSTEST(int *id, int *nperm, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, int *which1, int *which0, double *s0, double *f)
{
  double *r, *s, *ssort, *ex2, *stat, *statall, *e;
  double *diff;
  int i, j, k;
  int nrowa=*ngene, ncola=*nperm;

  diff=Calloc(*n1,double);
  r=Calloc(*ngene,double);
  s=Calloc(*ngene,double);
  ssort=Calloc(*ngene,double);
  ex2=Calloc(*ngene,double);
  stat=Calloc(*ngene,double);
  statall=Calloc((*ngene)*(*nperm),double);


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
	qsort((void*)ssort,*ngene,sizeof(double),compare11);
	
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

    /* puts the ngene*nperm matrix of scores into a vector to be passed over to kolmogoroff */
    for (j=0; j<*ngene; j++){
      statall[j*(*nperm)+k]=fabs(stat[j]);
    }
  }

  e=kolmogoroff(statall,nrowa,ncola);

  for (k=0; k<*nperm; k++){
    f[k]=e[k];
  }
  

  Free(diff);
  Free(r);
  Free(s);
  Free(ssort);
  Free(ex2);
  Free(stat);
  Free(statall);
  Free(e);
}






void correlationKSTEST(double *vecperm, int *nperm, double *matrix, int *ngene, int *nsample, double *f)
{
  double *ex0, *ex1, *ex20, *ex21, *exboth, *stat, *statall, *e;
  int i, j, k;
  int nrowa=*ngene, ncola=*nperm;

  ex0=Calloc(1,double);
  ex1=Calloc(*ngene,double);
  ex20=Calloc(1,double);
  ex21=Calloc(*ngene,double);
  exboth=Calloc(*ngene,double);
  stat=Calloc(*ngene,double);
  statall=Calloc((*ngene)*(*nperm),double);


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
    
    /* puts the ngene*nperm matrix of scores into a vector to be passed over to kolmogoroff */
    for (j=0; j<*ngene; j++){
      statall[j*(*nperm)+k]=fabs(stat[j]);
    }
  }

  e=kolmogoroff(statall,nrowa,ncola);

  for (k=0; k<*nperm; k++){
    f[k]=e[k];
  }

    
  Free(ex0);
  Free(ex1);
  Free(ex20);
  Free(ex21);
  Free(exboth);
  Free(stat);
  Free(statall);
  Free(e);
}
