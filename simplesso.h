// PREALPHA version

#include <stdlib.h>
#include <stdio.h>

double scalar(double *x,double *y, int n);
int firstnegative(double *y,int n);
double *rapport(double **An, double *W, double *x,double *bn,int n,int m);
void sostituisci(int *index,int h,int k,int n);
double **creaAn(double **A,int *index, int m,int n);
double *creatmpb(double *b,int *index,int n);
double *creabn(double *b,int *index, int m,int n);
double *prodotto(double **Inverse, double *b,int n);
int min(double *rap, int n);

double *Simplesso(int m, int n, double **A, double *b, double *c,int *index)
{
	int i;
	//calcolo inversa(indicizzata per colonne),  An (indicizzata per righe), tempb, bn
	double **An=creaAn(A,index,m,n);
	double *tempb=creatmpb(b,index,n);
	double *bn=creabn(b,index,m,n);
	double **Inverse=NULL;
	
	//calcolo h
	double *y=(double *) malloc (n*sizeof(double));
	for(i=0;i<n;i++)
		y[i]=scalar(c,Inverse[i],n);
	int h=firstnegative(y,n);
	free(y);	
	
	//calcolo x
	
	double *x=prodotto(Inverse,tempb,n);
	
	if (h==-1)
		return x;
	
	
	//rapporti
	double *rap=rapport(An,Inverse[h],x,bn,n,m-n);
	int k= min(rap, m-n);
	if (k==-1)
	{
		printf("errore pochi vincoli\n");
		return NULL;
	}
	//finale
	free(Inverse);free(tempb);free(bn);free(rap);free(x);free(An);
	sostituisci(index,h,k,n);
	return Simplesso(m,n,A,b,c,index);
}
	
double scalar(double *x,double *y, int n)
{
	double sum=0;
	int i;
	for(i=0; i<n; i++)
		sum += (x[i]*y[i]);
	return sum;
}

int firstnegative(double *y,int n)
{
	int i=0;
	while (y[i]>=0 && i<n)
		i++;
	if (i==n)
		return -1;
	else return i;
}

void sostituisci(int *index,int h,int k,int n)
{
	if (k>=h)
	{
		int i=0;
		while(index[i]<h)
			i++;
		i++;
		while(index[i]<k)
		{
			index[i-1]=index[i];
			i++;
		}
		index[i-1]=k;
		
	}
	else 
	{
		int i=n-1;
		while(index[i]>h)
			i--;
		i--;
		while(index[i]>k)
			index[i+1]=index[i];
		index[i]=k;
		
	}
		
}

double **creaAn(double **A,int *index, int m,int n)
{
	int i,j,k;
	j=k=0;
	n=m-n;
	double **An= (double **) malloc (n*sizeof(double *));
	for (i=0;i<m;i++)
		if (i!=index[k])
			An[j++]=A[i];
	else k++;
	return An;
}

double *creabn(double *b,int *index, int m,int n)
{
	int i,j,k;
	j=k=0;
	n=m-n;
	double *bn= (double *) malloc (n*sizeof(double));
	for (i=0;i<m;i++)
		if (i!=index[k])
			bn[j++]=b[i];
	else k++;
	return bn;
}

double *creatmpb(double *b,int *index,int n)
{
	double *tmpb= (double *) malloc (n*sizeof(double));
	int i;
	for (i=0; i<n; i++)
		tmpb[i]=b[index[i]];
	return tmpb;
}

double *prodotto(double **Inverse, double *b,int n)
{
	double *x=(double *) malloc (n*sizeof(double));
	int i,j;
	for (i=0; i<n;i++)
		for(j=0;j<n;j++)
			x[i]+=Inverse[j][i] * b[j];
	return x;
}

double *rapport(double **An, double *W, double *x,double *bn,int n,int m)
{
	double *rap = (double *) malloc (m*sizeof(double));
	int i;
	for(i=0;i<m; i++)
		rap[i]=-1*scalar(An[i],W,n);
	for (i=0;i<m;i++)
		if (rap[i]>0)
			rap[i]= (bn[i]-scalar(An[i],x,n))/rap[i];
	return rap;
}

int min(double *rap, int n)
{
	int i,indmin;
	indmin = -1;
	for (i=0;i<n;i++)
		if (rap[i]>=0 && rap[i]<rap[indmin])
			indmin=i;
	return indmin;
}
