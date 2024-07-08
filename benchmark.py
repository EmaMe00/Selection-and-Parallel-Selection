import subprocess
import random
import os

# Funzione per contare le righe di un file
def conta_righe(filename):
    with open(filename, 'r') as f:
        righe = sum(1 for line in f)
    return righe

def calcola_media_da_file(nome_file):
    try:
        with open(nome_file, 'r') as file:
            # Leggi tutte le righe dal file
            righe = file.readlines()
            
            # Converti i valori in float e calcola la somma
            valori = [float(riga.strip()) for riga in righe]
            somma = sum(valori)
            
            # Calcola la media
            media = somma / len(valori)
            
            return media
    
    except FileNotFoundError:
        print(f"Errore: il file '{nome_file}' non è stato trovato.")
        return None
    except ValueError:
        print("Errore: il file contiene valori non numerici.")
        return None
    except Exception as e:
        print(f"Errore imprevisto: {e}")
        return None

def scrivi_su_file(nome_file, testo):
    try:
        with open(nome_file, 'a+') as file:
            file.write(testo + '\n')
    
    except Exception as e:
        print(f"Errore durante la scrittura nel file '{nome_file}': {e}")

# Funzione per eliminare il file se esiste
def elimina_file_se_esiste(filename):
    if os.path.exists(filename):
        os.remove(filename)


# Determina il numero di righe nel file array.txt
nome_file_array = 'array.txt'

problem_size = [10000,100000,1000000,10000000,100000000,1000000000]

for probl_siz in problem_size:
    
    command = ["python3.11", "generator.py", str(probl_siz)]
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Errore durante l'esecuzione del comando: {e}")
        exit(1)
        
    valori = conta_righe(nome_file_array)
    valori_controllati = 100 #MODIFICARE SOLO LA PROPORSIONE RELATIVA ALLA QUANTITÀ DI ITERAZIONI
    valori_controllati = int(valori_controllati)
    print("Campiono " + str(valori_controllati) + " valori.")
    # Genera valori casuali per ith
    ith = [str(random.randint(1, valori - 1)) for _ in range(valori_controllati)]

    for sort in range(0,2):

        for processori in range (1,9):

            filename = f'bench_proc{processori}_value{valori}_sort{sort}.txt'
            count = 1
            for index, valore in enumerate(ith):
                if count % 50 == 0:
                    print("Iterazione: " + str(count))
                count += 1
                filename = f'bench_proc{processori}_value{valori}_sort{sort}.txt'
                mpi_command = ['mpirun', '-np', str(processori), './main', 'array.txt', valore, str(sort), filename]

                try:
                    subprocess.run(mpi_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    #subprocess.run(mpi_command, check=True)
                except subprocess.CalledProcessError:
                    print(f"Errore durante l'esecuzione {index + 1}.")
                    exit(1)
            
        
            media = calcola_media_da_file(filename)
            media = round(media, 6)
            testo = "Processori: " + str(processori) + ". Valori totali file: " + str(valori) + ". Sort: " + str(sort) + ". TEMPO MEDIO: " + str(media)
            scrivi_su_file("resoconto_benchmark.txt", testo)
            print("Esecuzione dei comandi completata per: ")
            print(testo)
            print()
            elimina_file_se_esiste(filename)
                
