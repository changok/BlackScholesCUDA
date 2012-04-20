include ../../Makefile.include

CPPFLAGS += -I../include
LDFLAGS += -ldcmt 
OBJS = example1.o example2.o example3.o log.o timer.o

all: example1 example2 example3 

%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

example1:   example1.o ../include/dc.h 
	$(CC) $(CFLAGS) $(CPPFLAGS) -o example1 example1.o $(LDFLAGS)

example2:   example2.o ../include/dc.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -o example2 example2.o $(LDFLAGS)

example3:   example3.o ../include/dc.h log.o timer.o
	$(CC) $(CFLAGS) $(CPPFLAGS) -o example3 example3.o log.o timer.o $(LDFLAGS)

log.o: log.c log.h

timer.o: timer.c timer.h

clean :
	/bin/rm -f example1 example2 example3 $(OBJS)


