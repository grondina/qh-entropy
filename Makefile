CC=gcc
CFLAGS=-std=c99 -Wall -Wextra -O3
LIBS=-lm -lcblas -llapacke
DEPS = data.h io.h
OBJS = data.o main.o io.o
EXEC = qh-entropy

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	rm -rf $(EXEC) $(OBJS)
