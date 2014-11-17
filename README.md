Simplex
=======

C implementation of Simplex Algorithm

Beta

simplesso primale

parametri dimensioni matrice A m e n, matrice A, vettori b e c, vettore indici di base

tutti i valori, compresi indici di base, cominciano da 0

restituisce vettore x o NULL per risultato infinito

limitato a n=2 (in attesa di funzione inversa per n>2)

modalità e limitazioni settabili nelle prime righe del file

singleinstance mode esegue una sola iterazione dell'algoritmo

verbose mode stampa durante l'esecuzione i valori di x, y, indici uscente ed entrante, nuova base e rapporti

verbose  e singleinstance mode attivi di deault
