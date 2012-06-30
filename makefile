LIB = libbcutils.so

all: default

default:
	gcc -c -fPIC -shared -rdynamic -fmax-errors=1 -std=gnu99 -g -O0 -Wall -m64 Array.c Hashmap.c bcutils.c
	gcc -o $(LIB) -fPIC -shared -rdynamic Array.o Hashmap.o bcutils.o -lm

clean:
	rm -f *.o *~ $(LIB) $(PROG)

install:
	install $(LIB) /usr/lib
	mkdir -p /usr/include/bcutils
	cp *.h /usr/include/bcutils
	chmod 644 /usr/include/bcutils/*
	ldconfig
