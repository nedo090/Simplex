// PREALPHA version

#include <stdlib.h>
#include <stdio.h>

double scalar(double *x,double *y, int n);
int firstnegative(double *y,int n);//testata
double *rapport(double **An, double *W, double *x,double *bn,int n,int m);
void sostituisci(int *index,int h,int k,int n);//testata
double **creaAn(double **A,int *index, int m,int n);
double *creatmpb(double *b,int *index,int n);//testata
double *creabn(double *b,int *index, int m,int n);//testata
double *prodotto(double **Inverse, double *b,int n);
int min(double *rap, int n);//testata

double *Simplesso(int m, int n, double **A, double *b, double *c,int *index)
{
	int i;
	//calcolo inversa(indicizzata per colonne),  An (indicizzata per righe), tempb, bn
	double **An=creaAn(A,index,m,n);
	double *tempb=creatmpb(b,index,n);
	double *bn=creabn(b,index,m,n);
	double **Inverse=NULL;
	
	//calcolo indice uscente h
	double *y=(double *) malloc (n*sizeof(double));
	for(i=0;i<n;i++)
		y[i]=scalar(c,Inverse[i],n);
	int h=firstnegative(y,n);
	free(y);	
	
	//calcolo x
	
	double *x=prodotto(Inverse,tempb,n);
	// in caso di vertice ottimo restituisce x
	if (h==-1)
		return x;
	
	
	//calcolo rapporti, risultato infinito e indice entrante k
	double *rap=rapport(An,Inverse[h],x,bn,n,m-n);
	int k= min(rap, m-n);
	if (k==-1)
	{
		printf("errore pochi vincoli\n");
		return NULL;
	}
	//libera memoria, aggiorna gli indici e reitera il simplesso
	free(Inverse);free(tempb);free(bn);free(rap);free(x);free(An);
	sostituisci(index,h,k,n);
	return Simplesso(m,n,A,b,c,index);
}
//prodotto scalare	
double scalar(double *x,double *y, int n)
{
	double sum=0;
	int i;
	for(i=0; i<n; i++)
		sum += (x[i]*y[i]);
	return sum;
}

//cerca primo elemento negativo
int firstnegative(double *y,int n)
{
	int i=0;
	while (y[i]>=0 && i<n)
		i++;
	if (i==n)
		return -1;
	else return i;
}

//scambia gli indici tenendo il vettore ordinato
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
		while(index[i]>k && i!=-1)
		{
			index[i+1]=index[i];
			i--;
		}
		//if (i==0)
		//	index[0]=k;
		 index[i+1]=k;
		
	}
		
}
//crea matrice indici non di base prendendo i puntatori direttamente dalla matrice originale
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

// crea vettore b indici non di base
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

//crea vettore b di base
double *creatmpb(double *b,int *index,int n)
{
	double *tmpb= (double *) malloc (n*sizeof(double));
	int i;
	for (i=0; i<n; i++)
		tmpb[i]=b[index[i]];
	return tmpb;
}


// calcola la x eseguendo il prodotto righe per colonna
//Inverse è indicizzata per colonna
double *prodotto(double **Inverse, double *b,int n)
{
	double *x=(double *) malloc (n*sizeof(double));
	int i,j;
	for (i=0; i<n;i++)
		for(j=0;j<n;j++)
			x[i]+=Inverse[j][i] * b[j];
	return x;
}

//calcola i rapporti, prima i denominatori, poi i numeratori validi
double *rapport(double **An, double *W, double *x,double *bn,int n,int m)
{
	double *rap = (double *) malloc (m*sizeof(double));
	int i;
	for(i=0;i<m; i++)
		rap[i]=-1*scalar(An[i],W,n);
	for (i=0;i<m;i++)
		if (rap[i]>0)
			rap[i]= (bn[i]-scalar(An[i],x,n))/rap[i];
		else rap[i]= -1.0;
	return rap;
}

//calcola im minimo indice maggiore di 0
int min(double *rap, int n)
{
	int i,indmin;
	i=0;
	while (rap[i]<0)
		i++;
	indmin = i;
	for (;i<n;i++)
		if (rap[i]>=0 && rap[i]<rap[indmin])
			indmin=i;
	return indmin;
}
