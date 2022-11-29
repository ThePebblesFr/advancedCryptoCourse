CC = gcc
EXEC = main
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)
CFLAGS = -g -Wall -Wextra -lm -lgmp

$(EXEC): $(OBJ)
	$(CC) -o $(EXEC) $^ $(CFLAGS)
	
%.o : %.c
	$(CC) -c $< -o $@ $(CFLAGS)


clean : 
	rm -f *.o core

mrproper : clean
	rm -f executable
