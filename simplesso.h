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
	//calcolo inversa(indicizzata per colonne) fittizia nella versione alpha
	//double **Inverse=
	double Inverse [][]={1,2,3,4}; n=2; //fittizi

	
	//calcolo indice uscente h
	double *y=(double *) malloc (n*sizeof(double));
	for(i=0;i<n;i++)
		y[i]=scalar(c,Inverse[i],n);
	int h=firstnegative(y,n);
	free(y);	
	
	//calcolo x e vettore b di base tempb
	double *tempb=creatmpb(b,index,n);
	double *x=prodotto(Inverse,tempb,n);
	free(tempb);
	// in caso di vertice ottimo restituisce x
	if (h==-1)
		{
			free(Inverse);free(x);
			return x;
		}
	
	
	//calcolo rapporti, risultato infinito e indice entrante k
	double **An=creaAn(A,index,m,n); //matrice A non di base
	double *bn=creabn(b,index,m,n);
	double *rap=rapport(An,Inverse[h],x,bn,n,m-n);
	int k= min(rap, m-n);
	free(Inverse);free(bn);free(rap);free(x);free(An);
	if (k==-1)
	{
		printf("errore pochi vincoli\n");
		return NULL;
	}
	//aggiorna gli indici e reitera il simplesso
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
	while (i<n && y[i]>=0)
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
		while(i!=-1 && index[i]>k)
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
	j=k=i=0;
	double **An= (double **) malloc ((m-n)*sizeof(double *));
	while(k<n)
		{ 
			for(;i<index[k];i++,j++)
				An[j]=A[i];
				k++;i++;
		}
		for(;i<m;i++,j++)
			An[j]=A[i];
				
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
	double *x=(double *) calloc (n,sizeof(double));
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
