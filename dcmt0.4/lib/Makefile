include ../../Makefile.include

OBJS = serialize.o init.o check32.o prescr.o mt19937.o eqdeg.o sieve.o genmtrand.o random_seed.o 
OUT = libdcmt.$(DYLIB_SUFFIX)

all: lib

lib : $(OBJS)
	$(LINKER) $(CFLAGS) $(LD_SHARED_FLAG) -o $(OUT) $(OBJS)

#	ar -cru libdcmt.a $(OBJS)

%.o : %.c
	$(CC) -c $(CFLAGS) $(PICFLAG) $(CPPFLAGS) $< -o $@

init.o : init.c dci.h 

check32.o : check32.c dci.h 

prescr.o : prescr.c dci.h 

mt19937.o : mt19937.c

eqdeg.o : eqdeg.c dci.h 

random_seed.o: random_seed.c ../include/dc.h

sieve.o : sieve.c dci.h 

genmtrand.o : genmtrand.c dci.h 

serialize.o : serialize.c ../include/dc.h ../include/serialize.h

clean :
	/bin/rm -f $(OBJS) $(OUT)

oclean :
	/bin/rm -f $(OBJS)

