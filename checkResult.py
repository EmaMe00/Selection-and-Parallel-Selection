import subprocess
import random
import os
import sys

# Funzione per confrontare i due file
def confronta_file(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        for line1, line2 in zip(f1, f2):
            if line1.strip() != line2.strip():
                return False
    return True

# Funzione per eliminare il file se esiste
def elimina_file_se_esiste(filename):
    if os.path.exists(filename):
        os.remove(filename)

# Funzione per contare le righe di un file
def conta_righe(filename):
    with open(filename, 'r') as f:
        righe = sum(1 for line in f)
    return righe

# Determina il numero di righe nel file array.txt
nome_file_array = 'array.txt'
valori = conta_righe(nome_file_array)
valori_controllati = 100 #MODIFICARE SOLO LA PROPORSIONE RELATIVA ALLA QUANTITÀ DI DATI DA CAMPIONARE
valori_controllati = int(valori_controllati)
print("Campiono " + str(valori_controllati) + " valori.")

# Genera valori casuali per ith
ith = [str(random.randint(1, valori - 1)) for _ in range(valori_controllati)]

for sort in range (0,2):

    for processori in range(1, 9):

        # Nomi dei file di output
        filename_mpi = f'checkP_proc{processori}_value{valori}_sort{sort}.txt'
        filename_select = f'checkS_proc{processori}_value{valori}_sort{sort}.txt'

        # Elimina i file di output se esistono
        elimina_file_se_esiste(filename_mpi)
        elimina_file_se_esiste(filename_select)

        count = 1
        # Esegui i comandi per ogni valore di ith
        for index, valore in enumerate(ith):
            if count % 50 == 0:
                print("Iterazione: " + str(count))
            count += 1
            mpi_command = ['mpirun', '-np', str(processori), './main', 'array.txt', str(valore), str(sort), filename_mpi]
            select_command = ['./select', 'array.txt', str(valore), filename_select]

            # Esegui il comando mpirun
            try:
                subprocess.run(mpi_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                #subprocess.run(mpi_command, check=True)
            except subprocess.CalledProcessError:
                print(f"Errore durante l'esecuzione del comando:")
                print(mpi_command)
                exit(1)

            # Esegui il comando select
            try:
                subprocess.run(select_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            except subprocess.CalledProcessError:
                print(f"Errore durante l'esecuzione del comando:")
                print(select_command)
                exit(1)

        print("Esecuzione dei comandi completata per: ")
        print("Processori: " + str(processori) + ". Valori totali file: " + str(valori) + ". Valori controllati: " + str(valori_controllati) + ". Sort: " + str(sort))

        # Verifica se i file sono uguali
        if confronta_file(filename_mpi, filename_select):
            print("I file sono uguali.")
            print()
            elimina_file_se_esiste(filename_mpi)
            elimina_file_se_esiste(filename_select)
        else:
            print("Errore: i file non coincidono. Interrompo l'esecuzione.")
            print()
            elimina_file_se_esiste(filename_mpi)
            elimina_file_se_esiste(filename_select)
            sys.exit(1)
