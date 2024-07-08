/* 
	-------- ESECUZIONE ALGORITMO ------------
	
	COMPILARE TRAMITE COMANDO: make

   	SE CORRECTNESS = 0 il comando per eseguire l'algoritmo è:
   	- ./select <nomeFile con i numeri> <i-esima statistica d'ordine>
   	- Es: ./select "array.txt" 6
   
   	SE CORRECTNESS = 1 il comando per eseguire l'algoritmo è:
   	- ./select <nomeFile con i numeri> <i-esima statistica d'ordine> <nomoFile di output>
   	- Es: ./select "array.txt" 6 "output.txt"
	
	ESECUZIONE TEST DI CORRETTEZZA:
	IMPOSTARE CORRECTNESS = 1 NEL FILE select.c E 
	IMPOSTARE CORRECTNESS = 1 E BENCHMARK = 0 NEL FILE main.c 
	COMPILARE TRAMITE COMANDO make 
	LANCIARE DA TERMINALE IL COMANDO: python checkResult.py
	
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>

#define LINE_DIMENSION 100
#define CORRECTNESS 0

int readFile(char *fileName, int **A, int *size);
int ceilingDivision(int a, int b);
void swap(int A[], int i, int j);
void insertionSort(int A[], int p, int r);
int partition_around(int A[], int p, int r, int x);
int select(int A[], int p, int r, int i);
void write_double_to_file_check(const char *filename, int value, int ith);


int main (int argc, char* argv[]){

    int *X = NULL; //Array con i valori
    int ith = 0;
    int N = 0;

	#if CORRECTNESS
    	if (argc < 4 || argc >= 5)
    	{
        
        	printf("Insert a valid number of arguments: <exec> <fileName> <ith> <outputFileName>\n");
        	return 0;
    	}
    #else
    	if (argc < 3 || argc >= 4)
    	{
        
        	printf("Insert a valid number of arguments: <exec> <fileName> <ith>\n");
        	return 0;
    	}
    #endif


    if(readFile(argv[1],&X,&N)){
        printf("Insert a valid file \n");
        return 0;
    }
    

    ith = atoi(argv[2]);
    if (ith <= 0 || ith > N)
    {
        printf("Insert a valid ith order statistics \n");
        return 0;
    }

    int result = select(X,0,N-1,ith);
    printf("Result select: %d\n",result);
    #if CORRECTNESS
   		write_double_to_file_check(argv[3], result, ith);
    #endif

    return 0;
}


int readFile(char *fileName, int **A, int *size) {

    FILE* file = fopen(fileName, "r");
    if (!file) {
        printf("Impossibile aprire il file %s in lettura.\n", fileName);
        return 1;
    }

    char line[LINE_DIMENSION];
    int dimension = 0;
    int* safe_pointer;

    while(fgets(line, sizeof(line), file)){
        char *token = strtok(line, "\n"); 
        while (token != NULL) {
            dimension++;
            safe_pointer = realloc(*A,dimension*sizeof(int));
            if (safe_pointer == NULL)
            {
                free(*A);
                return 1;
            } else {
                *A = safe_pointer;
            }
            (*A)[dimension-1] = atoi(token);
            //printf("%d\n", A[dimension-1]);
            fflush(stdout);
            token = strtok(NULL, "\n");
        }
    }
    *size = dimension;
    fclose(file);
    return 0;
}

int ceilingDivision(int a, int b) {
    return (a + b - 1) / b;
}

void swap(int A[], int i, int j) {
    int temp = A[i];
    A[i] = A[j];
    A[j] = temp;
}

void insertionSort(int A[], int p, int r) {
    for (int i = p + 1; i <= r; i++) {
        int key = A[i];
        int j = i - 1;
        while (j >= p && A[j] > key) {
            A[j + 1] = A[j];
            j--;
        }
        A[j + 1] = key;
    }
}

int partition_around(int A[], int p, int r, int x) {
    for (int i = p; i <= r; i++) {
        if (A[i] == x) {
            swap(A, i, r);
            break;
        }
    }
    int pivot = A[r];
    int i = p - 1;
    for (int j = p; j < r; j++) {
        if (A[j] <= pivot) {
            i++;
            swap(A, i, j);
        }
    }
    swap(A, i + 1, r);
    return i + 1;
}

int select(int A[], int p, int r, int i) {
    while ((r - p + 1) % 5 != 0) {
        for (int j = p + 1; j <= r; j++) {
            if (A[p] > A[j]) {
                swap(A, p, j);
            }
        }
        if (i == 1) {
            return A[p];
        }
        p = p + 1;
        i = i - 1;
    }

    int g = (r - p + 1) / 5;
    for (int j = 0; j < g; j++) {
        insertionSort(A, p + j * 5, p + j * 5 + 4);
        swap(A, p + j, p + j * 5 + 2);
    }

    int medianOfMedians = select(A, p, p + g - 1, ceilingDivision(g, 2));
    int q = partition_around(A, p, r, medianOfMedians);
    int k = q - p + 1;

    if (i == k) {
        return A[q];
    } else if (i < k) {
        return select(A, p, q - 1, i);
    } else {
        return select(A, q + 1, r, i - k);
    }
}

void write_double_to_file_check(const char *filename, int value, int ith) {
    // Apertura del file in modalità append ("a")
    // "a" apre il file per la scrittura, mantenendo il contenuto esistente, o crea un nuovo file se non esiste
    FILE *file = fopen(filename, "a");
    
    // Verifica se l'apertura del file è andata a buon fine
    if (file == NULL) {
        fprintf(stderr, "Errore: Impossibile aprire il file %s\n", filename);
        exit(1);
    }
    
    // Scrittura dei valori interi nel file
    fprintf(file, "%d %d\n", value, ith);
    
    // Chiusura del file
    if (fclose(file) != 0) {
        fprintf(stderr, "Errore: Impossibile chiudere il file %s\n", filename);
        exit(1);
    }

}