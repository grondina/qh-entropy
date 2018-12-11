CC=gcc
CFLAGS=-std=c99 -Wall -Wextra -O3
LIBS=-lm -lz -lcblas -llapacke
DEPS = data.h io.h traj.h ref.h util.h
OBJS = data.o main.o io.o traj.o ref.o util.o
EXEC = qhe

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	rm -rf $(EXEC) $(OBJS)
