EXEC   = landy-szalay

OPTIMIZE  =  -O2  -m64 


OPTIMIZE += $(OPTS)

OBJS   = main.o kdtree2.o write_shock_catalogues.o read_athena_header.o 

CC     = g++

INCL   = kdtree2.hpp write_shock_catalogues.hpp shock_data_types.hpp read_athena_header.hpp read_athena_tracers.h

LIBS   = -lm


CFLAGS  = $(OPTIMIZE) #-I/usr/local/include
CXXFLAGS = #-stdlib=libstdc++
CPPFLAGS  = $(OPTIMIZE) # -I/usr/local/include/ -I../ -stdlib=libstdc++
LDFLAGS = #-m64 -stdlib=libstdc++ #-shared-intel #-shared-libgcc 

$(EXEC): $(OBJS) 
	 $(CC) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)   

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

