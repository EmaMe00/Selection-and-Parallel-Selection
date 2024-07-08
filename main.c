/* 
	-------- ESECUZIONE ALGORITMO ------------
	
	COMPILARE TRAMITE COMANDO: make
    ASSICURARSI CHE LA MACRO C SIA IMPOSTATA A UN VALORE MAGGIORE UGUALE DI 2, DI DEFAULT SI CONSIGLIA DI LASCIARE 2

   	SE BENCHMARK = 0 E CORRECTNESS = 0 il comando per eseguire l'algoritmo è:
   	- mpirun -np <numero processori> ./main <nomeFile con i numeri> <i-esima statistica d'ordine> <0 per array disordinato, 1 per array ordinato>
   	- Es: mpirun -np 6 ./main "array.txt" 10 0
   
   	SE BENCHMARK = 1 E CORRECTNESS = 0 il comando per eseguire l'algoritmo è:
   	- mpirun -np <numero processori> ./main <nomeFile con i numeri> <i-esima statistica d'ordine> <0 per array disordinato, 1 per array ordinato> <nomoFile di output>
   	- Es: mpirun -np 6 ./main "array.txt" 10 0 "output.txt"
   
    SE BENCHMARK = 0 E CORRECTNESS = 1 il comando per eseguire l'algoritmo è:
   	- mpirun -np <numero processori> ./main <nomeFile con i numeri> <i-esima statistica d'ordine> <0 per array disordinato, 1 per array ordinato> <nomoFile di output>
   	- Es: mpirun -np 6 ./main "array.txt" 10 0 "output.txt"
   
   	NON DEVONO MAI ESSERE IMPOSTATI BENCHMARK = 1 E CORRECTNESS = 1
	
	ESECUZIONE TEST BENCHMARK:
	IMPOSTARE CORRECTNESS = 0 E BENCHMARK = 1 NEL FILE main.c 
	COMPILARE TRAMITE COMANDO make 
	LANCIARE DA TERMINALE IL COMANDO: python benchmark.py
	
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
#include <mpi.h>

#define LINE_DIMENSION 100 		//Macro utilizzata per la lettura da file dei dati
#define BENCHMARK 0  			//Macro utilizzata per effettuare i benchmark
#define CORRECTNESS 0			//Macro utilizzata per effettuare un check di correttezza dell'applicazione parallela
#define C 2                     //Valore sperimentale per il controllo del numero di iterazioni. C >= 2 per il corretto funzionamento dell'algoritmo anche con p = 1
#define MAXBASE 1000000			//Soglia massima per la base utilizzata da Radix Sort

int readFile(char *fileName, int **A, int *size);  //Lettura del file di input con i numeri da analizzare									
int select(int A[], int p, int r, int i); //Algoritmo di selezione sequenziale
void swap(int A[], int i, int j); 
int partition_around(int A[], int p, int r, int x); //Algoritmo per partizionare gli elementi attorno ad un pivot
int ceilingDivision(int a, int b); //Funzione per arrotondare
int partition(int A[], int p, int q); //Procedura classica di partition
void insertionSort(int A[], int p, int r); 
void printArray(int A[], int size, int printElement, char string[], int rank); //Funzione di debug per stampare un'array di interi
void printArrayF(float A[], int size, int printElement, char string[]); //Funzione di debug per stampare un'array di float
int weighted3median(int X[], float Y[], float S, float E, int p, int q, float W); //Funzione per determinare la lower weighted3median
int weighted3median_upper(int X[], float Y[], float S, float E, int p, int q, float W); //Funzione per determinare la uppert weighted3median

//Le funzioni con la W lavorano contemporaneamente sia sui dati che sui loro pesi, in modo che gli spostamenti degli elementi siano sincronizzati mantenendo la correlazione con gli indici
int selectW(int A[], float Y[], int p, int r, int i); 
void swapW(int A[], float Y[], int i, int j);
int partition_aroundW(int A[], float Y[], int p, int r, int x);
void insertionSortW(int A[], float Y[], int p, int r);

//Funzioni di partizione (il conteggio relativo agli elementi maggiori mino o uguali ai pivot viene effettuato al di fuori delle funzioni, quindi queste le ricevono come parametro)
void partition_2pivt(int A[], int n, int x, int y, int* lxi, int* exi, int* bxyi, int* eyi, int* gyi); 
int* five_way_partition(int A[], int n, int i, int x, int y,int* new_value, int* new_i, int* new_ni, int* new_N,int LX, int EX, int BXY, int EY, int GY, int lxi, int exi, int bxyi, int eyi, int gyi);
int* three_way_partition(int A[], int n, int i, int x, int *output_value, int *output_i, int *output_ni, int *output_N, int L, int G, int N, int li, int gi);

//Radix Sort che lavora in base 10 (DEBUG)
int getMax(int arr[], int n);
void countingSort(int arr[], int n, int exp); 
void radixSort(int arr[], int n);

//Funzioni di debug e di benchmark per scrivere su file
void write_double_to_file(const char *filename, double value);
void write_double_to_file_check(const char *filename, int value, int ith);

//Radix Sort che lavora su una base arbitraria (verrà utilizzato in base n in modo da garantire tempo lineare)
void radix_sort_b(int arr[], int n, int base);
void counting_sort_b(int arr[], int n, int base, int digit_place);
int get_digit_b(int number, int base, int digit_place);

int main(int argc, char* argv[]){

    MPI_Init(&argc, &argv);

    int N = 0; 					//Variabile che conterrà il numero di elementi totali presenti nei p processi in ogni iterazione, all'inizio N è il numero di valori presenti in X
    int *X = NULL; 				//Array contenente i valori

    int rank; 					//id del processo
    int p; 						//numero di processori
    int ith;				 	//i-esima statistica d'ordine
    int ith_for_print = 0;		//Variabile di debug
    int sorted = 0;			    //Utilizzo o meno di Radix Sort

    int* send_BUF = NULL;		//Buffer di trasferimento
    int* recv_BUF = NULL;		//Buffer di ricezione

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    //Il processo 0 è incaricato di acquisire i paramatri
    if (!rank)
    {
    	
        #if BENCHMARK || CORRECTNESS
        if (argc < 5 || argc >= 6)
        {
            printf("Insert a valid number of arguments: <exec> <fileName> <ith> <outputFileName> \n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            MPI_Finalize();
            return 0;
        }
        #else
        if (argc < 4 || argc >= 5)
        {
            
            printf("Insert a valid number of arguments: <exec> <fileName> <ith>\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            MPI_Finalize();
            return 0;
        }
        #endif

        sorted = atoi(argv[3]);
        if (sorted < 0 || sorted > 1)
        {
            printf("Insert 1 for sorted array or 0 for unsorted array \n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            MPI_Finalize();
            return 0;
        }


        if(readFile(argv[1],&X,&N)){
            printf("Insert a valid file \n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            MPI_Finalize();
            return 0;
        }
        
        if (N < p)
        {
            printf("N cannot be < p \n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            MPI_Finalize();
            return 0;
        }
        
        ith = atoi(argv[2]);
        ith_for_print = ith;
        if (ith <= 0 || ith > N)
        {
            printf("Insert a valid ith order statistics \n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            MPI_Finalize();
            return 0;
        }
        
        if (N/(C*p) < 2)
        {
            printf("constant C inappropriate, decrease it or p\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            MPI_Finalize();
            return 0;
        }
        
        if (C == 1 && p == 1)
        {
            printf("C and p cannot be 1, increase C or p \n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            MPI_Finalize();
            return 0;
        }
        

    }

    MPI_Bcast(&ith, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sorted, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int elements_per_proc = N/p; 
    int *sub_array = malloc(elements_per_proc * sizeof(int)); 
    if (sub_array == NULL) {
        fprintf(stderr, "Errore: Impossibile allocare memoria per sub_array\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        MPI_Finalize();
        return 1; 
    }

    // Scatter dell'array da processo 0 a tutti i processi
    MPI_Scatter(X, elements_per_proc, MPI_INT, sub_array, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

    // Il processo zero mantiene più valori rispetto agli altri se necessario
    if (!rank)
    {
        int extra = N % p; //Gli elementi che eccedono, se il numero di processori non divisibile per p vengono assegnati al processo 0
        if (extra > 0)
        {
            elements_per_proc = elements_per_proc + extra;
            int *tmp_pointer = realloc(sub_array, elements_per_proc * sizeof(int));
            if (tmp_pointer != NULL)
            {
                sub_array = tmp_pointer;
            }else {
                printf("Realloc error! \n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            int j = 0;
            for (int i = ((N/p)*p); i < N; i++)
            {
                sub_array[(N/p)+j] = X[i];
                j++;
            }
        }
        free(X);
    }
	
	//DEBUG
    // Tutti i processi stampano l'array ricevuto
    /*
    fflush(stdout);
    printf("Process %d ", rank);
    fflush(stdout);
    printArray(sub_array,elements_per_proc,1,"array");
    */

    MPI_Barrier (MPI_COMM_WORLD);
    double seconds = - MPI_Wtime();

    int *A = NULL;
    int new_value = 0, new_N = 0, new_i = 0, new_ni = 0;

    int n = N;
    int base = N;
	
	//Prima di iniziare la computazione ogni processore ordina il suo array
    if (sorted)
    {
    	if(base < 2 || base > MAXBASE){
    		base = MAXBASE;
    	}
		MPI_Bcast(&base, 1, MPI_INT, 0, MPI_COMM_WORLD);
        radix_sort_b(sub_array, elements_per_proc, base);
        //DEBUG
        //radixSort(sub_array, elements_per_proc);
    }
    
    
    while (N > (n/(p*C))){
    
        int m_i = 0;
        
        if(elements_per_proc > 0){
        	m_i = select(sub_array,0,elements_per_proc-1,ith);
        }
        
		//DEBUG
        //Stampo le mediane scelte dai singoli processi
        //printf("Process %d median: %d\n", rank, m_i);
        //fflush(stdout);

        send_BUF = malloc(2*sizeof(int));
        recv_BUF = malloc(2*p*sizeof(int));
        if (send_BUF == NULL || recv_BUF == NULL) {
            fprintf(stderr, "Errore: Impossibile allocare memoria\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            MPI_Finalize();
            return 1; // Esce dal programma con un codice di errore
        }
        
        send_BUF[0] = m_i;
        send_BUF[1] = elements_per_proc;

        // Esegue MPI_Allgather per raccogliere i dati da tutti i processi
        MPI_Allgather(send_BUF, 2, MPI_INT, recv_BUF, 2, MPI_INT, MPI_COMM_WORLD);

        int m_i_vec[p];
        float n_i_vec[p];

        fflush(stdout);
        for (int i = 0; i < p; i++)
        {
            m_i_vec[i] = recv_BUF[2*i];
            n_i_vec[i] = (float)recv_BUF[(2*i)+1];
        }
        free(send_BUF),free(recv_BUF);

		//DEBUG
        /* 
        if (!rank)
        {
            printArray(m_i_vec,p,1,"m_i");
            printArrayF(n_i_vec,p,1,"n_i");
        }
        */    
        

        int x,y;
        x = weighted3median(m_i_vec,n_i_vec,0,0,0,p-1,(float)N);
        y = weighted3median_upper(m_i_vec,n_i_vec,0,0,0,p-1,(float)N);
		//DEBUG
        /* 
        if (!rank)
        {
            printf("WM_L WM_U: %d %d\n", x,y);
            fflush(stdout);
        }
        */
        
        if (x != y)
        {
            int lxi = 0, exi = 0, bxyi = 0, eyi = 0, gyi = 0;
            if (!sorted)
            {
                partition_2pivt(sub_array, elements_per_proc, x, y, &lxi, &exi, &bxyi, &eyi, &gyi);
            }else{
                for (int i = 0; i < elements_per_proc; i++) {
                    if (sub_array[i] < x) {
                        lxi++;
                    } else if (sub_array[i] == x) {
                        exi++;
                    } else if (sub_array[i] > x && sub_array[i] < y) {
                        bxyi++;
                    } else if (sub_array[i] == y) {
                        eyi++;
                    } else if (sub_array[i] > y) {
                        gyi++;
                    }
                }
            }

            send_BUF = malloc(5*sizeof(int));
            recv_BUF = malloc(5*p*sizeof(int));
            if (send_BUF == NULL || recv_BUF == NULL) {
                fprintf(stderr, "Errore: Impossibile allocare memoria\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                MPI_Finalize();
                return 1; // Esce dal programma con un codice di errore
            }
            send_BUF[0] = lxi, send_BUF[1] = exi, send_BUF[2] = bxyi ,send_BUF[3] = eyi, send_BUF[4] = gyi;
            
            // Esegue MPI_Allgather per raccogliere i dati da tutti i processi
            MPI_Allgather(send_BUF, 5, MPI_INT, recv_BUF, 5, MPI_INT, MPI_COMM_WORLD);

            int LX = 0, EX = 0, BXY = 0, EY = 0, GY = 0;
            for (int i = 0; i < p; i++)
            {
                LX = LX + recv_BUF[(5*i)];
                EX = EX + recv_BUF[(5*i)+1];
                BXY = BXY + recv_BUF[(5*i)+2];
                EY = EY + recv_BUF[(5*i)+3];
                GY = GY + recv_BUF[(5*i)+4];
            }
            
            free(send_BUF),free(recv_BUF);
            //DEBUG
            /* 
            if (!rank)
            {
                printf("LX EX BXY EY GY: %d %d %d %d %d\n",LX,EX,BXY,EY,GY);
                fflush(stdout);
            }
            */
            A = five_way_partition(sub_array, elements_per_proc, ith, x, y, &new_value, &new_i, &new_ni, &new_N, LX, EX, BXY, EY, GY, lxi, exi, bxyi, eyi, gyi);
        } else {

            int L = 0, G = 0;

            if (!sorted)
            {
                partition_around(sub_array,0,elements_per_proc-1,x);
            }
            int li = 0, gi = 0;
            for (int j = 0; j < elements_per_proc; j++) {
                if (sub_array[j] < x) {
                    li++;
                } else if (sub_array[j] > x) {
                    gi++;
                }
            }
    
            send_BUF = malloc(2*sizeof(int));
            recv_BUF = malloc(2*p*sizeof(int));
            if (send_BUF == NULL || recv_BUF == NULL) {
                fprintf(stderr, "Errore: Impossibile allocare memoria\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                MPI_Finalize();
                return 1; // Esce dal programma con un codice di errore
            }
            send_BUF[0] = li, send_BUF[1] = gi;
            // Esegue MPI_Allgather per raccogliere i dati da tutti i processi
            MPI_Allgather(send_BUF, 2, MPI_INT, recv_BUF, 2, MPI_INT, MPI_COMM_WORLD);

            for (int i = 0; i < p; i++)
            {
                L = L + recv_BUF[2*i];
                G = G + recv_BUF[(2*i)+1];
            }
            free(send_BUF),free(recv_BUF);
		
			//DEBUG
            //printf("L G: %d %d\n",L,G);
            A = three_way_partition(sub_array, elements_per_proc, ith, x, &new_value, &new_i, &new_ni, &new_N, L, G, N, li, gi);
        }

        if (new_value != INT_MAX) {
            MPI_Barrier (MPI_COMM_WORLD);
            seconds += MPI_Wtime();
            if (!rank)
            {
                fflush(stdout);
                printf("\nResult P: %d\n", new_value);
                fflush(stdout);
                printf ("N = %d, ith =%d, Processes = %d, Time = %12.6f sec ",n, ith_for_print, p, seconds);
                #if CORRECTNESS
                    write_double_to_file_check(argv[4], new_value, ith_for_print);
                #endif
                #if BENCHMARK
                    write_double_to_file(argv[4], seconds);
                #endif
            }
            free(sub_array);
            free(A);
            MPI_Finalize();
            return 0;
        }

        ith = new_i;
        elements_per_proc = new_ni;
        if(sub_array != NULL){
        	free(sub_array);
        }
        sub_array = A;
        A = NULL;
        if(new_N == N){
        	break;
        }else{
        	N = new_N;
        }

        //DEBUG
        //printArray(sub_array, elements_per_proc, 1, "sub_array",rank);
        //printf("rank:%d  N:%d N/p:%d elements_per_proc:%d\n",rank, N,(n/(p)),elements_per_proc);
        
        //printf("N = %d\n", N);
		//printf("n/p*C = %d\n", n/(p*C));
    }
	//DEBUG
    //printf("\nIl processo %d è uscito dal while! \n",rank);
    
    //DEBUG
    //printf("new_N: %d\n",new_N);
    //printf("new_ni: %d\n",new_ni);
    
    int *final_recv_data = malloc(sizeof(int)*new_N);
    // Prepara i buffer per MPI_Gatherv
    int* recv_counts = NULL;
    int* displs = NULL;

    if (rank == 0) {
        recv_counts = malloc(p * sizeof(int));
        displs = malloc(p * sizeof(int));
    }

    // Raccoglie le dimensioni degli array inviati da ciascun processo
    MPI_Gather(&new_ni, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < p; i++) {
            displs[i] = displs[i - 1] + recv_counts[i - 1];
        }
    }

    // Utilizziamo MPI_Gatherv per inviare gli array al processo principale
    MPI_Gatherv(sub_array, elements_per_proc, MPI_INT, final_recv_data, recv_counts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        
        //DEBUG
        //printArray(final_recv_data, new_N, 1, "Array finale",rank);
        //printf("New i: %d\n",new_i);
        //fflush(stdout);
        
        int result = select(final_recv_data, 0, new_N - 1, new_i);
        fflush(stdout);
        printf("\nResult: %d\n", result);
        #if CORRECTNESS
            write_double_to_file_check(argv[4], result, ith_for_print);
        #endif
        fflush(stdout);

        free(recv_counts);
        free(displs);
        free(final_recv_data);
    }
    if (elements_per_proc > 0)
    {
        free(A), free(sub_array);
    }
    
    MPI_Barrier (MPI_COMM_WORLD);
    seconds += MPI_Wtime();

    if (!rank) {
        printf ("N = %d, ith =%d, Processes = %d, Time = %12.6f sec ",n, ith_for_print, p, seconds);
        fflush(stdout);
        #if BENCHMARK
            write_double_to_file(argv[4], seconds);
        #endif
   }
    

    MPI_Finalize();
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
            fflush(stdout);
            token = strtok(NULL, "\n");
        }
    }
    *size = dimension;
    fclose(file);
    printf("File read successfully! \n");
    return 0;
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

void swap(int A[], int i, int j) {
    int temp = A[i];
    A[i] = A[j];
    A[j] = temp;
}

int ceilingDivision(int a, int b) {
    return (a + b - 1) / b;
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

int partition(int A[], int p, int q) {
    int x = A[p]; // primo elemento come pivot
    int i = p;

    for (int j = p + 1; j <= q; j++) {
        if (A[j] <= x) {
            i = i + 1;
            swap(A, i, j);
        }
    }
    swap(A, p, i);
    return i;
}

void printArray(int A[], int size, int printElement, char string[], int rank){
    if (printElement)
    {
        printf("Rank: %d - %s (dim %d): ",rank,string,size);
        fflush(stdout);
        for (int i = 0; i < size; i++)
        {
            printf("%d ",A[i]);
            fflush(stdout);
        }
        printf("\n");
        fflush(stdout);
    }else {
        printf("Rank: %d - %s (dim %d): ",rank,string,size);
        fflush(stdout);
    }
}

void printArrayF(float A[], int size, int printElement, char string[]){
    if (printElement)
    {
        printf("%s (dim %d): ",string,size);
        fflush(stdout);
        for (int i = 0; i < size; i++)
        {
            printf("%.2f ",A[i]);
            fflush(stdout);
        }
        printf("\n");
        fflush(stdout);
    }else {
        printf("%s : %d\n",string,size);
        fflush(stdout);
    }
}

int weighted3median(int X[], float Y[], float S, float E, int p, int q, float W) {
    int length = q - p + 1;

    if (length == 1) {
        return X[p]; // Caso base: se c'è un solo elemento, lo restituisce
    }

    int median = (length % 2) ? (length / 2) + 1 : length / 2;
    int index_median = selectW(X, Y, p, q, median); // Trova l'indice del mediano tramite selectW
    int x_i = X[index_median]; // Valore mediano

    float A = 0.0f, B = 0.0f;
    float one_third_W = (1.0f / 3.0f) * W;
    float two_thirds_W = (2.0f / 3.0f) * W;

    // Calcola A (somma dei pesi degli elementi minori di x_i) e B (somma dei pesi degli elementi maggiori di x_i)
    for (int j = p; j <= q; j++) {
        if (X[j] < x_i) {
            A += Y[j];
        } else if (X[j] > x_i) {
            B += Y[j];
        }
    }

    // Controlla le condizioni di peso
    if (A + S <= one_third_W && B + E <= two_thirds_W) {
        return x_i; // Se soddisfa le condizioni, restituisce x_i come risultato
    }

    // Se A + S è troppo grande, riduce l'intervallo a sinistra di x_i
    if (A + S > one_third_W) {
        E = B + E + Y[index_median]; // Aggiorna E sommando il peso dell'elemento mediano
        return weighted3median(X, Y, S, E, p, index_median - 1, W); // Ricorsione sulla parte sinistra
    }

    // Se B + E è troppo grande, riduce l'intervallo a destra di x_i
    if (B + E > two_thirds_W) {
        S = A + S + Y[index_median]; // Aggiorna S sommando il peso dell'elemento mediano
        return weighted3median(X, Y, S, E, index_median + 1, q, W); // Ricorsione sulla parte destra
    }

    return x_i; // Caso di default, restituisce x_i
}

int weighted3median_upper(int X[], float Y[], float S, float E, int p, int q, float W) {
    int length = q - p + 1;

    if (length == 1) {
        return X[p];
    }

    int median = (length % 2) ? (length / 2) + 1 : length / 2;
    int index_median = selectW(X, Y, p, q, median);
    int x_i = X[index_median];

    float A = 0.0f, B = 0.0f;
    float two_thirds_W = (2.0f / 3.0f) * W;
    float one_third_W = (1.0f / 3.0f) * W;

    for (int j = p; j <= q; j++) {
        if (X[j] < x_i) {
            A += Y[j];
        } else if (X[j] > x_i) {
            B += Y[j];
        }
    }

    if (A + S <= two_thirds_W && B + E <= one_third_W) {
        return x_i;
    }

    if (A + S > two_thirds_W) {
        E = B + E + Y[index_median];
        return weighted3median_upper(X, Y, S, E, p, index_median - 1, W);
    }

    if (B + E > one_third_W) {
        S = A + S + Y[index_median];
        return weighted3median_upper(X, Y, S, E, index_median + 1, q, W);
    }

    return x_i;
}

int selectW(int A[], float Y[], int p, int r, int i) {
    while ((r - p + 1) % 5 != 0) {
        for (int j = p + 1; j <= r; j++) {
            if (A[p] > A[j]) {
                swapW(A, Y, p, j);
            }
        }
        if (i == 1) {
            return p;
        }
        p = p + 1;
        i = i - 1;
    }
    
    int g = (r - p + 1) / 5;
    for (int j = 0; j < g; j++) {
        insertionSortW(A, Y, p + j * 5, p + j * 5 + 4);
        swapW(A, Y, p + j, p + j * 5 + 2);
    }
    
    int medianOfMedians = selectW(A, Y, p, p + g - 1, ceilingDivision(g, 2));
    int q = partition_aroundW(A, Y, p, r, medianOfMedians);
    int k = q - p + 1;
    
    if (i == k) {
        return q;
    } else if (i < k) {
        return selectW(A, Y, p, q - 1, i);
    } else {
        return selectW(A, Y, q + 1, r, i - k);
    }
}

void swapW(int A[], float Y[], int i, int j) {
    int tempA = A[i];
    A[i] = A[j];
    A[j] = tempA;
    
    float tempY = Y[i];
    Y[i] = Y[j];
    Y[j] = tempY;
}

int partition_aroundW(int A[], float Y[], int p, int r, int x) {
    int i;
    for (i = p; i < r; i++) {
        if (A[i] == x) {
            break;
        }
    }
    swapW(A, Y, i, r);
    int pivot = A[r];
    int j = p - 1;
    for (int i = p; i <= r - 1; i++) {
        if (A[i] <= pivot) {
            j++;
            swapW(A, Y, j, i);
        }
    }
    swapW(A, Y, j + 1, r);
    return j + 1;
}

void insertionSortW(int A[], float Y[], int p, int r) {
    for (int i = p + 1; i <= r; i++) {
        int keyA = A[i];
        float keyY = Y[i];
        int j = i - 1;
        while (j >= p && A[j] > keyA) {
            A[j + 1] = A[j];
            Y[j + 1] = Y[j];
            j = j - 1;
        }
        A[j + 1] = keyA;
        Y[j + 1] = keyY;
    }
}

void partition_2pivt(int A[], int n, int x, int y, int* lxi, int* exi, int* bxyi, int* eyi, int* gyi) {
    if (n > 0)
    {
        *lxi = *exi = *bxyi = *eyi = *gyi = 0;
        int *partitioned = malloc(sizeof(int)*n);
        int index = 0;

        for (int i = 0; i < n; i++) {
            if (A[i] < x) {
                partitioned[index++] = A[i];
                (*lxi)++;
            } else if (A[i] == x) {
                partitioned[index++] = A[i];
                (*exi)++;
            } else if (A[i] > x && A[i] < y) {
                partitioned[index++] = A[i];
                (*bxyi)++;
            } else if (A[i] == y) {
                partitioned[index++] = A[i];
                (*eyi)++;
            } else if (A[i] > y) {
                partitioned[index++] = A[i];
                (*gyi)++;
            }
        }

        for (int i = 0; i < n; i++) {
            A[i] = partitioned[i];
        }
        free(partitioned);
    }
}

int* five_way_partition(int A[], int n, int i, int x, int y,int* new_value, int* new_i, int* new_ni, int* new_N,int LX, int EX, int BXY, int EY, int GY, int lxi, int exi, int bxyi, int eyi, int gyi) {
    
    int *B = NULL;
    int B_counter = 0;

    if (i <= LX) {
        B = malloc(sizeof(int)*lxi);
        //discard all of the elements ≥ x;
        for (int j = 0; j < n; j++) {
            if (A[j] < x) {
                B[B_counter] = A[j];
                B_counter++;
            }
        }
        *new_i = i;
        *new_ni = lxi;
        *new_N = LX;
        *new_value = INT_MAX;
        return B;
    } else if (i <= LX + EX && i > LX) {
        *new_value = x;
        return B;
    } else if (i <= LX + EX + BXY && i > LX + EX) {
        B = malloc(sizeof(int)*bxyi);
        // Rimuovi gli elementi <= x e >= y
        for (int j = 0; j < n; j++) {
            if (A[j] > x && A[j] < y) {
                B[B_counter] = A[j];
                B_counter++;
            }
        }
        *new_i = i - (LX + EX);
        *new_ni = bxyi;
        *new_N = BXY;
        *new_value = INT_MAX;
        return B;
    } else if (i <= LX + EX + BXY + EY && i > LX + EX + BXY) {
        *new_value = y;
        return B;
    } else {
        B = malloc(sizeof(int)*gyi);
        // Rimuovi gli elementi <= y
        for (int j = 0; j < n; j++) {
            if (A[j] > y) {
                B[B_counter] = A[j];
                B_counter++;
            }
        }
        *new_i = i - (LX + EX + BXY + EY);
        *new_ni = gyi;
        *new_N = GY;
        *new_value = INT_MAX;
        return B;
    }

    *new_value = INT_MAX;
    return B;
}

int* three_way_partition(int A[], int n, int i, int x, int *output_value, int *output_i, int *output_ni, int *output_N, int L, int G, int N, int li, int gi) {
    
    int *B = NULL;
    int B_counter = 0;

    if (i <= L) {
        B = malloc(sizeof(int)*li);
        // Scarta tutti gli elementi >= x
        for (int j = 0; j < n; j++) {
            if (A[j] < x) {
                B[B_counter] = A[j];
                B_counter++;
            }
        }
        *output_value = INT_MAX; // ∞
        *output_i = i;
        *output_ni = li;
        *output_N = L;
        return B;
    } else if (L < i && i <= N - G) {
        // Ritorna x, indicando che i è già stato trovato
        *output_value = x;
        *output_i = INT_MAX; // ∞
        *output_ni = INT_MAX; // ∞
        *output_N = INT_MAX; // ∞
        return B;
    } else {
        B = malloc(sizeof(int)*gi);
        // Scarta tutti gli elementi <= x
        for (int j = 0; j < n; j++) {
            if (A[j] > x) {
                B[B_counter] = A[j];
                B_counter++;
            }
        }
        i = i - (N - G);
        *output_value = INT_MAX; // ∞
        *output_i = i;
        *output_ni = gi;
        *output_N = G;
        return B;
    }
    return B;
}

int getMax(int arr[], int n) {
    int max = arr[0];
    for (int i = 1; i < n; i++) {
        if (arr[i] > max) {
            max = arr[i];
        }
    }
    return max;
}

void countingSort(int arr[], int n, int exp) {
    int *output = malloc(sizeof(int)*n); // array di output
    int count[10] = {0}; // array di conteggio inizializzato a 0

    // Memorizza la frequenza delle cifre in count[]
    for (int i = 0; i < n; i++) {
        count[(arr[i] / exp) % 10]++;
    }

    // Modifica count[] in modo che contenga le posizioni effettive di queste cifre in output[]
    for (int i = 1; i < 10; i++) {
        count[i] += count[i - 1];
    }

    // Costruisce l'array di output
    for (int i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }

    // Copia l'array di output nell'array originale arr[]
    for (int i = 0; i < n; i++) {
        arr[i] = output[i];
    }
    free(output);
}

void radixSort(int arr[], int n) {
    // Trova il numero massimo per conoscere il numero di cifre
    int m = getMax(arr, n);

    // Esegui counting sort per ogni cifra. Inizia con la cifra delle unità (exp=1), poi decine (exp=10), centinaia, ecc.
    for (int exp = 1; m / exp > 0; exp *= 10) {
        countingSort(arr, n, exp);
    }
}

void write_double_to_file_check(const char *filename, int value, int ith) {

    FILE *file = fopen(filename, "a");
    
    if (file == NULL) {
        fprintf(stderr, "Errore: Impossibile aprire il file %s\n", filename);
        exit(1);
    }
    
    fprintf(file, "%d %d\n", value, ith);
    

    if (fclose(file) != 0) {
        fprintf(stderr, "Errore: Impossibile chiudere il file %s\n", filename);
        exit(1);
    }

}

void write_double_to_file(const char *filename, double value) {

    FILE *file = fopen(filename, "a");
    

    if (file == NULL) {
        fprintf(stderr, "Errore: Impossibile aprire il file %s\n", filename);
        exit(1);
    }
    

    fprintf(file, "%lf\n", value);
    

    if (fclose(file) != 0) {
        fprintf(stderr, "Errore: Impossibile chiudere il file %s\n", filename);
        exit(1);
    }
}

int get_digit_b(int number, int base, int digit_place) {
    return (number / (int)pow(base, digit_place)) % base;
}

void counting_sort_b(int arr[], int n, int base, int digit_place) {
    // Array di conteggio per le cifre (da 0 a base-1)
    int *count = (int *)calloc(base, sizeof(int));
    if (count == NULL) {
        fprintf(stderr, "Errore di allocazione di memoria\n");
        return;
    }

    // Conta le occorrenze di ogni cifra all'interno dell'array
    for (int i = 0; i < n; i++) {
        int digit = get_digit_b(arr[i], base, digit_place);
        count[digit]++;
    }

    // Aggiorna il conteggio per mantenere la posizione corretta di ogni cifra
    for (int i = 1; i < base; i++) {
        count[i] += count[i - 1];
    }

    // Costruzione dell'array di output in base alla cifra corrente
    int *output = (int *)malloc(n * sizeof(int));
    if (output == NULL) {
        fprintf(stderr, "Errore di allocazione di memoria\n");
        free(count);
        return;
    }

    for (int i = n - 1; i >= 0; i--) {
        int digit = get_digit_b(arr[i], base, digit_place);
        output[count[digit] - 1] = arr[i];
        count[digit]--;
    }

    // Copia l'array di output nell'array principale
    for (int i = 0; i < n; i++) {
        arr[i] = output[i];
    }

    free(output);
    free(count);
}

void radix_sort_b(int arr[], int n, int base) {
    // Trova il numero massimo nell'array
    int max = arr[0];
    for (int i = 1; i < n; i++) {
        if (arr[i] > max) {
            max = arr[i];
        }
    }

    // Trova il numero massimo di cifre in base 'base'
    int max_digits = 0;
    if (max > 0) {
        max_digits = (int)(log(max) / log(base)) + 1;
    }

    // Applica Radix Sort per ogni cifra, da destra a sinistra
    for (int digit_place = 0; digit_place < max_digits; digit_place++) {
        counting_sort_b(arr, n, base, digit_place);
    }
}