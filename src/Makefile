CC=gcc
CFLAGS=-O2 -g

SRCS =  bendit.c expoint14.c trinuc_hellserv.c util.c
OBJS = $(SRCS:.c=.o)

all:	clean compile
	echo OK!!!!!!!!!

clean:
	rm -rf ../bin/bendIt *~ *.o *.so
	
compile: $(OBJS)
	$(CC) *.o -o ../bin/bendIt -lm
	rm -rf *.o
	chmod 755 ../bin/bendIt

