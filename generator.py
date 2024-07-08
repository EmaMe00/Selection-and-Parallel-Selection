import os
import random
import argparse
from tqdm import tqdm

def generate_random_numbers(file_name, count):
    lower_bound = 1
    upper_bound = count
    
    # Se il file esiste, lo elimina
    if os.path.exists(file_name):
        os.remove(file_name)
    
    # Genera numeri casuali e li scrive nel file uno alla volta
    with open(file_name, 'w') as file:
        for _ in tqdm(range(count), desc="Generazione dei numeri"):
            number = random.randint(lower_bound, upper_bound)
            file.write(f"{number}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genera un file con numeri casuali.')
    parser.add_argument('max_value', type=int, help='Il valore massimo dei numeri generati.')

    args = parser.parse_args()
    
    generate_random_numbers("array.txt", args.max_value)
    print(f"File 'array.txt' generato con {args.max_value} numeri casuali compresi tra 1 e {args.max_value}.")
