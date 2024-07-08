# Definisci i compilatori
MPICC = mpicc
CC = gcc

# Definisci le opzioni del compilatore
CFLAGS = -Wall -O2

# Nome dei file eseguibili
MPI_TARGET = main
CC_TARGET = select

# Lista dei file sorgente
MPISRCS = main.c
CCSRCS = select.c

# Regola predefinita
all: $(MPI_TARGET) $(CC_TARGET)

# Regola per creare l'eseguibile MPI
$(MPI_TARGET): $(MPISRCS)
	$(MPICC) $(CFLAGS) -o $(MPI_TARGET) $(MPISRCS)

# Regola per creare l'eseguibile C normale
$(CC_TARGET): $(CCSRCS)
	$(CC) $(CFLAGS) -o $(CC_TARGET) $(CCSRCS)

# Regola per pulire i file generati dalla compilazione
clean:
	rm -f $(MPI_TARGET) $(CC_TARGET) *.o

# Regola per rigenerare tutto da zero
rebuild: clean all
