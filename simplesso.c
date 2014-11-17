// Beta version

//test mode
#define limit 2
#define verbose
#define singleinstance
//
#include <stdlib.h>
#include "inverse.c"

#ifdef verbose
#include <stdio.h>
#endif

double scalar(double *x,double *y, int n);
int firstnegative(double *y,int n);
double *rapport(double **An, double *W, double *x,double *bn,int n,int m);
void sostituisci(int *index,int h,int k,int n);
double **creaAn(double **A,int *index, int m,int n);
double *creatmpb(double *b,int *index,int n);
double *creabn(double *b,int *index, int m,int n);
double *prodotto(double **Inverse, double *b,int n);
int min(double *rap, int n);
double **selectbase(double **A,int *index,int n);
int normalize(int k, int *index, int n);

double *Simplesso(int m, int n, double **A, double *b, double *c,int *index)
{
	int i;
	#ifdef limit
	n=limit; //limitazione testing version
	#endif
	//calcolo inversa(indicizzata per colonne) 
	double **base=selectbase(A,index,n);
	double **Inverse=inversion(base);
	free(base);

	
	//calcolo indice uscente h
	double *y=(double *) malloc (n*sizeof(double));
	for(i=0;i<n;i++)
		y[i]=scalar(c,Inverse[i],n);
	int h=firstnegative(y,n);
	#ifdef vervose
	printf("\n y= ");
	for(i=0;i<n;i++)
		printf("%f",y[i]);
	printf("\n");
	#endif
	free(y);	
	
	//calcolo x e vettore b di base tempb
	double *tempb=creatmpb(b,index,n);
	double *x=prodotto(Inverse,tempb,n);
	free(tempb);
	// in caso di vertice ottimo restituisce x
	#ifdef vervose
	printf("\n x= ");
	for(i=0;i<n;i++)
		printf("%f",x[i]);
	printf("\n");
	#endif
	
	if (h==-1)
		{
			free(Inverse);free(x);
			return x;
		}
	h=index[h]; //assegna il valore effettivo dell'indice uscente
	
	//calcolo rapporti, risultato infinito e indice entrante k
	double **An=creaAn(A,index,m,n); //matrice A non di base
	double *bn=creabn(b,index,m,n);
	double *rap=rapport(An,Inverse[h],x,bn,n,m-n);
	int k = min(rap, m-n);
	free(Inverse);free(bn);free(rap);free(x);free(An);
	if (k==-1)
	{
		#ifdef verbose
		printf("errore pochi vincoli\n");
		#endif
		return NULL;
	}
	k = normalize(k,index,n); //normalizza indice entrante col valore effettivo
	//aggiorna gli indici e reitera il simplesso
	sostituisci(index,h,k,n);
	
	#ifdef verbose
	printf("\nindice entrante=%d\nindice uscente=%d\nnuova bese:",k,h);
	for(i=0;i<n;i++)
		printf("%d ",index[i]);
	printf("\n");
	#endif
	
	#ifdef singleinstance
	return NULL;
	#endif 
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
	#ifdef verbose
	printf("rapporti:");
	#endif
	double *rap = (double *) malloc (m*sizeof(double));
	int i;
	for(i=0;i<m; i++)
		rap[i]=-1*scalar(An[i],W,n);
	for (i=0;i<m;i++)
		if (rap[i]>0)
			{
				rap[i]= (bn[i]-scalar(An[i],x,n))/rap[i];
				#ifdef verbose
				printf("%f ",rap[i]);
				#endif
			}
		else rap[i]= -1.0;
	return rap;
}

//calcola im minimo indice maggiore di 0
int min(double *rap, int n)
{
	int i,indmin;
	i=0;
	while (i<n && rap[i]<0)
		i++;
	
	if (i==n)
		return -1;

	indmin = i;
	for (;i<n;i++)
		if (rap[i]>=0 && rap[i]<rap[indmin])
			indmin=i;
	return indmin;
}

//seleziona i vettori della base dalla matrice A
double **selectbase(double **A,int *index,int n)
{
	double **base = (double **)malloc(n*sizeof(double *));
	int i;
	for(i=0;i<n;i++)
	base[i]= A[index[i]];
	return base;
}

//normalizza indice entrante tenendo conto delle righe di base sottratte alla matrice iniziale
int normalize(int k, int *index, int n)
{
	int i=0;
	while (i<n && k>=index[i])
	{
		if (k<index[i])
			i++;
		else{
				i++;k++;
		   	}
		
	}
	return k;
}
