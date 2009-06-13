EXEC   = P-GroupFinder
#HDF_DIR = /sw
HDF_DIR = /usr/local/hdf5/1.8.2s

#OPT1   = -DFOF_ONLY -DUSE_HDF5 # does only FOF group finding
OPT1   = -DUSE_HDF5 # does only FOF group finding

OPTIMIZE = -O3

OPTIONS = $(OPTIMIZE) $(OPT1) -g


OBJS = main.o io_input.o io_enzo.o io_hdf.o io_hdf5.o cmpfunc.o allvars.o \
	allocate.o ngbtree.o forcetree.o subfind.o density.o unbind.o \
	find_subgroups.o order_by_potential.o properties.o nrsrc/nrutil.o \
	nrsrc/sort2_int.o nrsrc/sort_int.o nrsrc/iindexx.o nrsrc/selectb.o \
	nrsrc/indexx.o nrsrc/sort2_flt_int.o 

CC     = mpicc

INCL   = parameters.h allvars.h proto.h ngbtree.h forcetree.h \
	 nrsrc/nrutil.h

CFLAGS = $(OPTIONS) -I/usr/include/ -I$(HDF_DIR)/include

LIBS = -L$(HDF_DIR)/lib 		\
	-lm -lhdf5 -lz -lsz
#LIBS = -L$(HDF_DIR)/lib 		\
#	-limf -lm -lmfhdf -ldf -ljpeg -lsz -lz	\
#	-lmpi

$(EXEC): $(OBJS) 
	 $(CC) $(OBJS) $(LIBS) -o $(EXEC)  

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

