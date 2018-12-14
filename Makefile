CC=gcc
CFLAGS=-std=c99 -Wall -Wextra -O3 -fopenmp
LIBS=-lm -lz -lcblas -llapacke
DEPS = data.h io.h par.h traj.h util.h
OBJS = data.o main.o io.o par.o traj.o util.o
EXEC = qhe

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	rm -rf $(EXEC) $(OBJS)
