CC=gcc
CFLAGS=-Wall -Wextra -O3
LIBS=-lm -lcblas -llapacke
DEPS = io.h
OBJS = main.o io.o
EXEC = qh-entropy

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	rm -rf $(EXEC) $(OBJS)
