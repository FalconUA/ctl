TARGET = demo
CC = g++
CFLAGS = -std=c++11 
HEADERS = cryptlib.hpp ctlfac.hpp ctllog.hpp ctlmod.hpp typespec.hpp other.h
OBJECTS = \
	demo.o \
	primitive_root_modulo.o \
	fac_fermat.o fac_lenstra.o fac_pollard.o fac_qsieve.o fac_lenstra.o \
	log_primitive.o log_shanks.o log_pollard.o log_index.o log_polhig_hellman.o log_adleman.o \

demo: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	rm -r *.o *.gch $(TARGET)
