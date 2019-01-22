CFLAGS=-Wall -D EULER -I.
CC=gcc
LIBS=-lm
DEPS = dg1D.h
SOURCES=main.c meshInit.c calcFlux.c output.c shape.c project.c
OBJ=$(SOURCES:.c=.o)
EXEC=dgrk
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o *~

