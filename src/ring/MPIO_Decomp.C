// Parallel particle sorting and distribution for ParallelParticleIO
 
// Robert Harkness
// April 2008
 
 
#include <mpi.h>
#include <hdf5.h>
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
 
#define MAXCPU 16384
#define FAIL 1
#define MAX_LINE_LENGTH          512
#define MAX_GRID_TAG_SIZE         16
#define MAX_TASK_TAG_SIZE         16
#define GRID_TAG_FORMAT        "4.4"
#define TASK_TAG_FORMAT        "4.4"
 
#include "macros_and_parameters.h"
 
// HDF5 prototypes
 
#include "mpio_extern_hdf5.h"
 
int Enzo_Dims_create(int nnodes, int ndims, int *dims); 
 
 
Eint32 main(Eint32 argc, char *argv[])
{
 
  hid_t       file_id, dset_id, attr_id;
  hid_t       mem_dsp_id, file_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id;
  hid_t       int_file_type_id, int_mem_type_id;
  hid_t       file_acc_template;
  hid_t       xfer_prop_list;
 
  hsize_t     dims[3];
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  hsize_t     mem_stride, mem_count, attr_count;
  hsize_t     slab_stride[3], slab_count[3];
  hsize_t     Slab_Rank;
  hsize_t     Slab_Dims[2];
 
  hssize_t    mem_offset;
  hssize_t    slab_offset[3];
 
  int i, j, k, m, n;
  int a, b, c;
  int ii, jj, kk;
  int ic, ipc, jpc, ppc, iblk;
  int hits;
  int fln, ext, lex;
 
  int thisnode, prevnode, nextnode, ltype, rtype, stype, ier;
 
  int dim, rank, ngrids, gridcounter;

  int mode;
  char *argin;
  char s1[5];
  char proc_name[MAX_LINE_LENGTH];
 
  int enzo_layout[3] = {0,0,0};
 
  int Starts[3];
  int Ends[3];
  int TopGridDims[3];
 
  double DomainLeftEdge[3] = {0.0, 0.0, 0.0};
  double DomainRightEdge[3] = {1.0, 1.0, 1.0};
  double SubDomainLeftEdge[3];
  double SubDomainRightEdge[3];
  double Left[3];
  double Right[3];
  double x0, y0, z0, dx, dy, dz, dpos;
 
  double GridLeft[MAXCPU][3];
  double GridRight[MAXCPU][3];

  double *outbuff = NULL;;
 
  int dbuff_size;
  int TotalParticleCount;
  int NumberOfParticles;
  int ParticleRank;
 
  double Start_Wall_Time, End_Wall_Time, Total_Wall_Time;
 
  MPI_Request req1[100];
  MPI_Request req2[100];
  MPI_Status  stat[100];
 
  char *PPin;
  char *PVin;
  char *PMin;
  char *PTin;
  char *Extension;
 
  int io_log = 1;
  int io_log_d = 0;
 
//
 
  MPI_Arg mpi_argc;
  MPI_Arg mpi_size;
  MPI_Arg mpi_rank;
 
  MPI_Datatype Type;
  MPI_Arg Count;
  MPI_Arg Source;
  MPI_Arg Dest;
  MPI_Arg Tag;

  MPI_Arg lname;

//  MPI_Arg mpi_layout[3] = {0,0,0};
//  MPI_Arg ncpu, jcpu
  
  int mpi_layout[3] = {0,0,0};
  int ncpu, jcpu;

//
 
  mpi_argc = argc;
 
  MPI_Init(&mpi_argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Get_processor_name(proc_name, &lname);

  fprintf(stderr, "Proc %d of %d is %s\n", mpi_rank, mpi_size, proc_name);

 
  ncpu = mpi_size;
  jcpu = mpi_rank;
//  argc = mpi_argc;

/* 
  if ( argc != 5 )
  {
    fprintf(stderr, "USAGE: ring ParticlePositionFile ParticleVelocityFile ParticleMassFile ParticleTypeFile\n");
    MPI_Finalize();
    return FAIL;
  }
*/

/*
  Mode = 1 then Position and Velocity only
  Mode = 2 then Position, Velocity and Mass
  Mode = 3 then Position, Velocity and Type
  Mode = 4 then Position, Velocity, Mass and Type

  If Mode = 1, skip mass and type
  If Mode = 2, skip type
  If Mode = 3, skip mass

  If Mode = 2 or 4, do mass
  If Mode = 3 or 4, do type
*/

/* 
  PPin = argv[1];
  PVin = argv[2];
  PMin = argv[3];
  PTin = argv[4];
*/

  argin = argv[1];
  mode = 0;

  if (argc == 4) {
    strcpy(s1, "pv");
    i = strcmp(argin, s1);
    if ( i == 0 ) {
      mode = 1;
      PPin = argv[2];
      PVin = argv[3];
      if ( mpi_rank == 0 ) {
        printf("Input arg %s should be pv\n", argin);
        printf("PPin = %s\n", PPin);
        printf("PVin = %s\n", PVin);
      }
    }
  }

  if (argc == 5) {

    strcpy(s1, "pvm");
    i = strcmp(argin, s1);
    if ( i == 0 ) {
      mode = 2;
      PPin = argv[2];
      PVin = argv[3];
      PMin = argv[4];
      if ( mpi_rank == 0 ) {
        printf("Input arg %s should be pvm\n", argin);
        printf("PPin = %s\n", PPin);
        printf("PVin = %s\n", PVin);
        printf("PMin = %s\n", PMin);
      }
    }

    strcpy(s1, "pvt");
    i = strcmp(argin, s1);
    if ( i == 0 ) {
      mode = 3;
      PPin = argv[2];
      PVin = argv[3];
      PTin = argv[4];
      if ( mpi_rank == 0 ) {
        printf("Input arg %s should be pvt\n", argin);
        printf("PPin = %s\n", PPin);
        printf("PVin = %s\n", PVin);
        printf("PTin = %s\n", PTin);
      }
    }

  }

  if (argc == 6) {

    strcpy(s1, "pvmt");
    i = strcmp(argin, s1);
    if ( i == 0 ) {
      mode = 4;
      PPin = argv[2];
      PVin = argv[3];
      PMin = argv[4];
      PTin = argv[5];
      if ( mpi_rank == 0 ) {
        printf("Input arg %s should be pvmt\n", argin);
        printf("PPin = %s\n", PPin);
        printf("PVin = %s\n", PVin);
        printf("PMin = %s\n", PMin);
        printf("PTin = %s\n", PTin);
      }
    }

    strcpy(s1, "pvtm");
    i = strcmp(argin, s1);
    if ( i == 0 ) {
      mode = 4;
      PPin = argv[2];
      PVin = argv[3];
      PTin = argv[4];
      PMin = argv[5];
      if ( mpi_rank == 0 ) {
        printf("Input arg %s should be pvtm\n", argin);
        printf("PPin = %s\n", PPin);
        printf("PVin = %s\n", PVin);
        printf("PMin = %s\n", PMin);
        printf("PTin = %s\n", PTin);
      }
    }

  }

  // printf("Mode = %d\n", mode);

  if (mode == 0) {
    if ( mpi_rank == 0 ) {
      printf("Error: no matching keys\n");
    }
  }

  if (mode > 0) {
    if ( mpi_rank == 0 ) {
      printf("Read Position\n");
      printf("Read Velocity\n");
    }
  }

  if ((mode == 2) || (mode == 4)) {
    if ( mpi_rank == 0 ) {
      printf("Read Mass\n");
    }
  }

  if ((mode == 3) || (mode == 4)) {
    if ( mpi_rank == 0 ) {
      printf("Read Type\n");
    }
  }

 
  FILE *log;
  char pid[MAX_GRID_TAG_SIZE];
 
  lex =0;
  fln = strlen(PPin);
  ext = strcspn(PPin, ".");
 
  if ( fln-ext > 0 )
  {
    lex = strlen(strstr(PPin, "."));
    Extension = new char[MAX_LINE_LENGTH];
    strcpy(Extension, strstr(PPin, "."));
    lex = strlen(Extension);
  }
 
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, jcpu);
 
  char *pgname = new char[MAX_LINE_LENGTH];
  strcpy(pgname, "PGlog");
  strcat(pgname, pid);
  if (lex > 0)
    strcat(pgname, Extension);
 
  log = fopen(pgname, "a");
 
  char *PPos = new char[MAX_LINE_LENGTH];
  strcpy(PPos, "PPos");
  strcat(PPos, pid);
  if (lex > 0)
    strcat(PPos, Extension);
 
  char *PVel = new char[MAX_LINE_LENGTH];
  strcpy(PVel, "PVel");
  strcat(PVel, pid);
  if (lex > 0)
    strcat(PVel, Extension);
 
  char *PPro = new char[MAX_LINE_LENGTH];
  strcpy(PPro, "PPro");
  strcat(PPro, pid);
  if (lex > 0)
    strcat(PPro, Extension);
 
  char *PMass = new char[MAX_LINE_LENGTH];
  strcpy(PMass, "PMass");
  strcat(PMass, pid);
  if (lex > 0)
    strcat(PMass, Extension);
 
  char *PType = new char[MAX_LINE_LENGTH];
  strcpy(PType, "PType");
  strcat(PType, pid);
  if (lex > 0)
    strcat(PType, Extension);
 
  mem_type_id = HDF5_R8;
  file_type_id = HDF5_FILE_R8;
 
  file_id = H5Fopen(PPin, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (io_log) fprintf(log, "H5Fopen with Name = %s\n", PPin);
    assert( file_id != h5_error );
 
  dset_id = H5Dopen(file_id, PPin);
    if (io_log) fprintf(log, "H5Dopen with Name = %s\n", PPin);
    assert( dset_id != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "Component_Rank");
    if (io_log) fprintf(log, "H5Aopen with Name = Component_Rank\n");
    assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, &ParticleRank);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "Component_Size");
   if (io_log) fprintf(log, "H5Aopen with Name = Component_Size\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, &TotalParticleCount);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_count = 3;
 
  attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
    assert( attr_dsp_id != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "TopGridStart");
   if (io_log) fprintf(log, "H5Aopen with Name = TopGridStart\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, Starts);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "TopGridEnd");
   if (io_log) fprintf(log, "H5Aopen with Name = TopGridEnd\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, Ends);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  attr_id = H5Aopen_name(dset_id, "TopGridDims");
   if (io_log) fprintf(log, "H5Aopen with Name = TopGridDims\n");
   assert( attr_id != h5_error );
 
  h5_status = H5Aread(attr_id, HDF5_INT, TopGridDims);
    if (io_log) fprintf(log, "H5Aread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(attr_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  // Divide particles among cpus and introduce fake particles
  // to keep FakePC % Ncpu = 0
  // Set the xyz values to -1.0 so they are excluded
 
  int FakeParticleCount;
  int Fakes;
  int ParticlesPerCPU;
  int nfake;
 
  if ( TotalParticleCount % ncpu == 0 )
  {
    ParticlesPerCPU = TotalParticleCount/ncpu;
    FakeParticleCount = ncpu * ParticlesPerCPU;
    Fakes = 0;
  }
  else
  {
    ParticlesPerCPU = TotalParticleCount/ncpu + 1;
    FakeParticleCount = ncpu * ParticlesPerCPU;
    Fakes = FakeParticleCount - TotalParticleCount;
  }
 
  fprintf(stderr, "PROC %"ISYM" : ParticlesPerCPU = %"ISYM", FakeParticleCount = %"ISYM", Fakes = %"ISYM"\n",
          jcpu, ParticlesPerCPU, FakeParticleCount, Fakes);
 
 
 
 
  dbuff_size = ParticlesPerCPU;
 
  double *buff[2];
  buff[0] = new double[dbuff_size];
  buff[1] = new double[dbuff_size];
 
  rank = 3;
  ngrids = 1;
 
  // MPI_Dims_create(ncpu, rank, mpi_layout);

  Enzo_Dims_create(ncpu, rank, mpi_layout);
  
/*
  if (rank == 3 && ncpu == 64)
  {
    for (dim = 0; dim < rank; dim++)
    {
      mpi_layout[dim] = 4;
    }
    printf("NCPU = 64 ==> coerced to 4**3\n");
  }
*/
 
  for (dim = 0; dim < rank; dim++)
  {
    enzo_layout[dim] = mpi_layout[rank-1-dim];
    ngrids *= enzo_layout[dim];
  }

  if ( jcpu == 0 ) {
    fprintf(stderr, "NumberOfGrids = %"ISYM"\n", ngrids);
    fprintf(stderr, "ENZO_layout %"ISYM" x %"ISYM" x %"ISYM"\n", enzo_layout[0], enzo_layout[1], enzo_layout[2]);
  }
 
  fprintf(log, "NumberOfGrids = %"ISYM"\n", ngrids);
  fprintf(log, "ENZO_layout %"ISYM" %"ISYM" %"ISYM"\n", enzo_layout[0], enzo_layout[1], enzo_layout[2]);
 
  fprintf(log, "TopGridDims = %"ISYM" %"ISYM" %"ISYM"\n", TopGridDims[0], TopGridDims[1], TopGridDims[2]);
  fprintf(log, "TopGridStart = %"ISYM" %"ISYM" %"ISYM"\n", Starts[0], Starts[1], Starts[2]);
  fprintf(log, "TopGridEnd = %"ISYM" %"ISYM" %"ISYM"\n", Ends[0], Ends[1], Ends[2]);
 
  for (dim = 0; dim < rank; dim++)
  {
    SubDomainLeftEdge[dim] = Starts[dim] * (DomainRightEdge[dim]-DomainLeftEdge[dim])/((double) TopGridDims[dim]);
    SubDomainRightEdge[dim] = (Ends[dim]+1) * (DomainRightEdge[dim]-DomainLeftEdge[dim])/((double) TopGridDims[dim]);
    //  SubCellWidth[dim] = (SubDomainRightEdge[dim]-SubDomainLeftEdge[dim])/((double) enzo_layout[dim]);
  }
 
  x0 = SubDomainLeftEdge[0];
  y0 = SubDomainLeftEdge[1];
  z0 = SubDomainLeftEdge[2];
 
  dx = (SubDomainRightEdge[0]-SubDomainLeftEdge[0])/((double) enzo_layout[0]);
  dy = (SubDomainRightEdge[1]-SubDomainLeftEdge[1])/((double) enzo_layout[1]);
  dz = (SubDomainRightEdge[2]-SubDomainLeftEdge[2])/((double) enzo_layout[2]);
 
  a=enzo_layout[0];
  b=enzo_layout[1];
  c=enzo_layout[2];
 
  gridcounter = 0;
 
 
  for (kk = 0; kk < enzo_layout[2]; kk++)
    for (jj = 0; jj < enzo_layout[1]; jj++)
      for (ii = 0; ii < enzo_layout[0]; ii++)
      {
        n=gridcounter;
 
        if ( n == jcpu )
        {
      // rank to coordinate
        m = n;
        i = m/(b*c);
        m = m%(b*c);
        j = m/c;
        m = m%c;
        k = m;
      // coordinate to rank check
        m = ((i*b*c) + j*c) + k;
        fprintf(log,"Grid %"ISYM"  {%"ISYM" %"ISYM" %"ISYM"}  %"ISYM"\n",n,i,j,k,m);
 
        Left[0] =  x0 + dx * (double) ii;
        Right[0] = x0 + dx * (double) (ii+1);
        Left[1] =  y0 + dy * (double) jj;
        Right[1] = y0 + dy * (double) (jj+1);
        Left[2] =  z0 + dz * (double) kk;
        Right[2] = z0 + dz * (double) (kk+1);
 
        GridLeft[jcpu][0] = Left[0];
        GridLeft[jcpu][1] = Left[1];
        GridLeft[jcpu][2] = Left[2];
 
        GridRight[jcpu][0] = Right[0];
        GridRight[jcpu][1] = Right[1];
        GridRight[jcpu][2] = Right[2];
 
        for (m = 0; m < rank; m++)
          fprintf(log, "Grid %"ISYM"    Left   %16.8"FSYM"   Right  %16.8"FSYM"\n", gridcounter, Left[m], Right[m]);
 
        }
 
        gridcounter++;
 
      }
 
  // Particle mask array - one bit per particle
 
  Eunsigned_int *BitArray = NULL;
  Eunsigned_int *BitMask = NULL;
  Eunsigned_int BitsPerInt;
  Eunsigned_int BitMaskSize;
  Eunsigned_int BitMaskTrue;
  Eunsigned_int TestBit;
  Eunsigned_int MaskAddr;
  Eunsigned_int WordAddr;
  Eunsigned_int BitAddr;
 
 
  int XMask;
 
 
  BitsPerInt = 8 * sizeof(BitsPerInt);
  BitMaskSize = (FakeParticleCount/BitsPerInt)+1;
 
  BitMask = new Eunsigned_int[BitMaskSize];
  BitArray = new Eunsigned_int[BitsPerInt];
 
  for (i = 0; i < BitMaskSize; i++)
    BitMask[i] = 0;
 
  BitArray[BitsPerInt-1] = 1;
  BitMaskTrue = 1;
 
  for (i=BitsPerInt-1; i>0; i--)
  {
    BitArray[i-1] = 2 * BitArray[i];
    BitMaskTrue = (BitMaskTrue | BitArray[i-1]);
  }
 
  for (i = 0; i < BitMaskSize; i++)
    BitMask[i] = BitMaskTrue;
 
  hits = 0;  // count direct hits on grid boundaries
 
  Start_Wall_Time = MPI_Wtime();
 
  for ( dim = 0; dim < rank; dim++)
  {
 
    dims[0] = 3;
    dims[1] = TotalParticleCount;
 
    mem_type_id = HDF5_R8;

#ifdef MPIO

    //  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);
    //    if (io_log) fprintf(log, "H5Pcreate file_access_template\n");
    //    assert( file_acc_template != h5_error );
 
    //  h5_status = H5Pset_fapl_mpio(file_acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);
    //    if (io_log) fprintf(log, "H5Pset_fapl_mpio: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );

#else
 
    file_acc_template = H5P_DEFAULT;
      if (io_log) fprintf(log, "Default file_access_template\n");

#endif
 
    file_id = H5Fopen(PPin, H5F_ACC_RDONLY, file_acc_template);
      if (io_log) fprintf(log, "H5Fopen with Name = %s\n", PPin);
      assert( file_id != h5_error );
 
    dset_id = H5Dopen(file_id, PPin);
      if (io_log) fprintf(log, "H5Dopen with Name = %s\n", PPin);
      assert( dset_id != h5_error );
 
    mem_offset = 0;
    mem_count = ParticlesPerCPU;
    mem_stride = 1;
 
    if ( jcpu == ncpu - 1 )
    {
      mem_count = ParticlesPerCPU - Fakes;
    }
 
    mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
      assert( mem_dsp_id != h5_error );
 
    h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
      if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    file_dsp_id = H5Screate_simple((Eint32) 2, dims, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
      assert( file_dsp_id != h5_error );
 
    slab_offset[0] = dim;   // x,y,z
    slab_stride[0] = 1;
    slab_count[0] = 1;
 
    slab_offset[1] = jcpu * ParticlesPerCPU;
    slab_stride[1] = 1;
    slab_count[1] = ParticlesPerCPU;
 
    if ( jcpu == ncpu - 1 )
    {
      slab_count[1] = ParticlesPerCPU - Fakes;
    }
 
    h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
      if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

#ifdef MPIO
 
    //  xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
    //    if (io_log) fprintf(log, "H5Pcreate xfer_prop_list\n");
    //    assert( xfer_prop_list != h5_error );
 
    //  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
    //    if (io_log) fprintf(log, "H5Pset_dxpl_mpio: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );

#else
 
    xfer_prop_list = H5P_DEFAULT;
      if (io_log) fprintf(log, "Default xfer_prop_list\n");

#endif
 
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, buff[0]);
      assert( h5_status != h5_error );

#ifdef MPIO
 
    //  h5_status = H5Pclose(xfer_prop_list);
    //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );

#endif
 
    h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(mem_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Fclose(file_id);
      if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

#ifdef MPIO
 
    //  h5_status = H5Pclose(file_acc_template);
    //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
    //    assert( h5_status != h5_error );

#endif
 
    if ( jcpu == ncpu - 1 )
    {
       for ( nfake = 0; nfake < Fakes; nfake++ )
       {
         buff[0][ParticlesPerCPU - Fakes + nfake] = -1.0;
       }
    }
 
 
    if (io_log_d)
    {
          fprintf(log, "Proc %"ISYM" Dim %"ISYM"\n", jcpu, dim);
          for ( i = 0; i < ParticlesPerCPU; i++ )
          {
            fprintf(log, "%8.4"FSYM, buff[0][i]);
            if ( ((i+1) % 16) == 0 )
              fprintf(log, "\n");
          }
          fprintf(log, "\n");
    }
 
 
    MPI_Barrier(MPI_COMM_WORLD);
 
    thisnode = jcpu;
    nextnode = (thisnode+1) % ncpu;
    prevnode = (thisnode-1+ncpu) % ncpu;
 
    fprintf(log, "P:  %"ISYM" %"ISYM" %"ISYM"\n", prevnode, thisnode, nextnode);
 
    stype = 30000+thisnode;
    rtype = 30000+prevnode;
    ltype = 30000+nextnode;
 
    for (k = 0; k < ncpu; k++)
    {
 
    i = (k  ) % 2;
    j = (k+1) % 2;
 
    MPI_Barrier(MPI_COMM_WORLD);
 
  /* Transmit right
    ier = MPI_Irecv( buff[j], ParticlesPerCPU, MPI_DOUBLE, prevnode, rtype, MPI_COMM_WORLD, req1 );
    if (io_log) fprintf(log, "IRECV proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, j, ier);
 
    ier = MPI_Isend( buff[i], ParticlesPerCPU, MPI_DOUBLE, nextnode, stype, MPI_COMM_WORLD, req2 );
    if (io_log) fprintf(log, "SEND proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, i, ier);
  */
 
  /* Transmit left */
 
    Count = ParticlesPerCPU;
    Source = nextnode;
    Tag = ltype;
    Type = MPI_DOUBLE;
 
//    ier = MPI_Irecv( buff[j], ParticlesPerCPU, MPI_DOUBLE, nextnode, ltype, MPI_COMM_WORLD, req1 );
    ier = MPI_Irecv( buff[j], Count, Type, Source, Tag, MPI_COMM_WORLD, req1 );
    if (io_log) fprintf(log, "IRECV proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, j, ier);
 
    Count = ParticlesPerCPU;
    Dest = prevnode;
    Tag = stype;
    Type = MPI_DOUBLE;
 
//    ier = MPI_Isend( buff[i], ParticlesPerCPU, MPI_DOUBLE, prevnode, stype, MPI_COMM_WORLD, req2 );
    ier = MPI_Isend( buff[i], Count, Type, Dest, Tag, MPI_COMM_WORLD, req2 );
    if (io_log) fprintf(log, "SEND proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, i, ier);
 
    // work on buff[i]
 
    if (io_log) fprintf(log, "WORK on buffer[%"ISYM"]\n", i);
 
/*
    if (io_log_d)
    {
          fprintf(log, "Proc %"ISYM" Dim %"ISYM" K %"ISYM" Buffer%"ISYM"\n", jcpu, dim, k, i);
          for (int idim = 0; idim < ParticlesPerCPU; idim++ )
          {
            fprintf(log, "%8.4"FSYM, buff[i][idim]);
            if ( ((idim+1) % 16) == 0 )
              fprintf(log, "\n");
          }
          fprintf(log, "\n");
    }
*/
 
  /* Transmit right
    ipc = (jcpu-k+ncpu) % ncpu;
  */
 
  /* Transmit left */
    ipc = (jcpu+k+ncpu) % ncpu;
 
    ipc = ipc * ParticlesPerCPU;
    jpc = ipc + ParticlesPerCPU - 1;
 
    if (io_log) fprintf(log, "ipc = %"ISYM" to %"ISYM"\n", ipc, jpc);
    if (io_log) fprintf(log, "left %6.4"FSYM", right %6.4"FSYM"\n", Left[dim], Right[dim]);
 
    for (ic = 0; ic < ParticlesPerCPU; ic++)
    {
      dpos = (double) buff[i][ic];
 
      if ( dpos == Left[dim] )
      {
         fprintf(stderr, "Particle on left boundary [%"ISYM"] [%"ISYM"] %16.8"FSYM"\n", i, ic, buff[i][ic]);
         hits = hits + 1;
      }
 
      if ( dpos == Right[dim] )
      {
         fprintf(stderr, "Particle on right boundary [%"ISYM"] [%"ISYM"] %16.8"FSYM"\n", i, ic, buff[i][ic]);
         hits = hits + 1;
      }
 
      if ( buff[i][ic] < Left[dim] || buff[i][ic] > Right[dim] )
      {
    //  Mask[ipc] = 0;
        MaskAddr = ipc;
        WordAddr = MaskAddr/BitsPerInt;
        BitAddr  = MaskAddr%BitsPerInt;
        TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
        if ( TestBit != 0 )
          BitMask[WordAddr] = (BitMask[WordAddr] ^ BitArray[BitAddr]);
      }
      ipc++;
    }
 
/*
    if (io_log_d)
    {
      for (ic = 0; ic < FakeParticleCount; ic++)
      {
        fprintf(log, "%1"ISYM, Mask[ic]);
        if ( ((ic+1) % 128) == 0 )
          fprintf(log, "\n");
      }
      fprintf(log, "\n");
    }
*/

    ier = MPI_Wait( req2, stat );
    if (io_log) fprintf(log, "WAIT2 proc %"ISYM" err %"ISYM"\n", thisnode, ier);
    ier = MPI_Wait( req1, stat );
    if (io_log) fprintf(log, "WAIT1 proc %"ISYM" err %"ISYM"\n", thisnode, ier);
 
    MPI_Barrier(MPI_COMM_WORLD);
 
    } // End loop over ring
 
  } // End loop over {x,y,z}
 
  End_Wall_Time = MPI_Wtime();
  Total_Wall_Time = End_Wall_Time - Start_Wall_Time;
  printf("%"ISYM": %16.8"FSYM" %16.8"FSYM" %16.8"FSYM"\n", jcpu, Start_Wall_Time, End_Wall_Time, Total_Wall_Time);
 
 
 
 
// FAIL here if any particles are on the mesh
// Do not count the fakes
 
  assert( hits == 0 );
 
// CHOP if(0)
// CHOP {
 
  NumberOfParticles = 0;
 
  for (i = 0; i < FakeParticleCount; i++)
  {
    XMask = 1;
    MaskAddr = i;
    WordAddr = MaskAddr/BitsPerInt;
    BitAddr  = MaskAddr%BitsPerInt;
    TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
    if ( TestBit == 0 )
      XMask = 0;
    if (io_log_d)
    {
      fprintf(log, "%1"ISYM, XMask);
      if ( ((i+1) % 128) == 0 )
      {
        fprintf(log, "\n");
      }
    }
    NumberOfParticles += XMask;
  }
 
  if (io_log_d)
    fprintf(log, "\n");
 
  fprintf(log, "Grid %"ISYM":  NumberOfParticles %"ISYM"\n", jcpu, NumberOfParticles);
 
/*
  if (io_log_d)
  {
    for (i = 0; i < FakeParticleCount; i++)
    {
      fprintf(log, "%1"ISYM, Mask[i]);
      if ( ((i+1) % 128) == 0 )
        fprintf(log, "\n");
    }
    fprintf(log, "\n");
  }
*/
 
 
 
 
//  Read particle data under bit mask

  if (mode > 0) {
 
  if (io_log) fprintf(log, "Read particle positions under bit mask\n");
 
//  double *outbuff = NULL;
 
  outbuff = new double[NumberOfParticles];
 
  for (dim = 0; dim < rank; dim++)
  {
 
  ppc = 0; // output particle counter
 
  dims[0] = 3;
  dims[1] = TotalParticleCount;
 
  mem_type_id = HDF5_R8;

#ifdef MPIO
 
  //  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);
  //    if (io_log) fprintf(log, "H5Pcreate file_acc_template\n");
  //    assert( file_acc_template != h5_error );
 
  //  h5_status = H5Pset_fapl_mpio(file_acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);
  //    if (io_log) fprintf(log, "H5Pset_fapl_mpio: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#else
 
  file_acc_template = H5P_DEFAULT;
    if (io_log) fprintf(log, "Default file_acc_template\n");

#endif
 
  file_id = H5Fopen(PPin, H5F_ACC_RDONLY, file_acc_template);
    if (io_log) fprintf(log, "H5Fopen with Name = %s\n", PPin);
    assert( file_id != h5_error );
 
  dset_id = H5Dopen(file_id, PPin);
    if (io_log) fprintf(log, "H5Dopen with name = %s\n", PPin);
    assert( dset_id != h5_error );
 
  mem_offset = 0;
  mem_count = ParticlesPerCPU;
  mem_stride = 1;
 
  if ( jcpu == ncpu - 1 )
  {
    mem_count = ParticlesPerCPU - Fakes;
  }
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );
 
  h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  file_dsp_id = H5Screate_simple((Eint32) 2, dims, NULL);
   if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
   assert( file_dsp_id != h5_error );
 
  slab_offset[0] = dim;   // x,y,z
  slab_stride[0] = 1;
  slab_count[0] = 1;
 
  slab_offset[1] = jcpu * ParticlesPerCPU;
  slab_stride[1] = 1;
  slab_count[1] = ParticlesPerCPU;
 
  if ( jcpu == ncpu - 1 )
  {
    slab_count[1] = ParticlesPerCPU - Fakes;
  }
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
  //    if (io_log) fprintf(log, "H5Pcreate: %"ISYM"\n", xfer_prop_list);
  //    assert( xfer_prop_list != h5_error );
 
  //  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
  //    if (io_log) fprintf(log, "H5Pset_dxpl_mpio: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#else
 
  xfer_prop_list = H5P_DEFAULT;
    if (io_log) fprintf(log, "Default xfer_prop_list\n");

#endif
 
  h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, buff[0]);
    if (io_log) fprintf(log, "H5Dread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  h5_status = H5Pclose(xfer_prop_list);
  //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#endif
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  h5_status = H5Pclose(file_acc_template);
  //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#endif
 
/*
  if ( jcpu == ncpu - 1 )
  {
     for ( nfake = 0; nfake < Fakes; nfake++ )
     {
       buff[0][ParticlesPerCPU - Fakes + nfake] = -1.0;
     }
  }
*/
 
  MPI_Barrier(MPI_COMM_WORLD);
 
  thisnode = jcpu;
  nextnode = (thisnode+1) % ncpu;
  prevnode = (thisnode-1+ncpu) % ncpu;
 
  stype = 30000+thisnode;
  rtype = 30000+prevnode;
  ltype = 30000+nextnode;
 
  iblk = - ( (ncpu - jcpu) % ncpu );
 
  for (k = 0; k < 2*ncpu; k++) // two passes
  {
 
  i = (k  ) % 2;
  j = (k+1) % 2;
 
  MPI_Barrier(MPI_COMM_WORLD);
 
  // Transmit left
 
  Count = ParticlesPerCPU;
  Source = nextnode;
  Tag = ltype;
  Type = MPI_DOUBLE;
 
//  ier = MPI_Irecv( buff[j], ParticlesPerCPU, MPI_DOUBLE, nextnode, ltype, MPI_COMM_WORLD, req1 );
  ier = MPI_Irecv( buff[j], Count, Type, Source, Tag, MPI_COMM_WORLD, req1 );
  if (io_log) fprintf(log, "IRECV proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, j, ier);
 
  Count = ParticlesPerCPU;
  Dest = prevnode;
  Tag = stype;
  Type = MPI_DOUBLE;
 
//  ier = MPI_Isend( buff[i], ParticlesPerCPU, MPI_DOUBLE, prevnode, stype, MPI_COMM_WORLD, req2 );
  ier = MPI_Isend( buff[i], Count, Type, Dest, Tag, MPI_COMM_WORLD, req2 );
  if (io_log) fprintf(log, "SEND proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, i, ier);
 
  // work on buff[i]
 
  if (io_log) fprintf(log, "K = %"ISYM", IBLK = %"ISYM"\n", k, iblk);
 
  if ( iblk > -1 && iblk < ncpu )
  {
 
    ipc = (jcpu+k+ncpu) % ncpu;
 
    ipc = ipc * ParticlesPerCPU;
    jpc = ipc + ParticlesPerCPU - 1;
 
    if (io_log) fprintf(log, "  ipc = %"ISYM" to %"ISYM"\n", ipc, jpc);
 
    for (ic = 0; ic < ParticlesPerCPU; ic++)
    {
      MaskAddr = ipc;
      WordAddr = MaskAddr/BitsPerInt;
      BitAddr  = MaskAddr%BitsPerInt;
      TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
      if ( TestBit == 0 )
        XMask = 0;
      if ( TestBit != 0 )
      {
        XMask = 1;
        outbuff[ppc] = buff[i][ic];
        ppc++;
      }
      ipc++;
    }
 
  }
 
  iblk++;

  ier = MPI_Wait( req2, stat ); 
  ier = MPI_Wait( req1, stat );
 
  if (io_log) fprintf(log, "WAIT proc %"ISYM" err %"ISYM"\n", thisnode, ier);
 
  } // end loop over rings
 
  fprintf(log, "L %8.4"FSYM"  R %8.4"FSYM"\n", GridLeft[jcpu][dim], GridRight[jcpu][dim]);
 
  if (io_log_d)
  {
    for (i = 0; i < NumberOfParticles; i++)
    {
      fprintf(log, "%8.4"FSYM, outbuff[i]);
      if ( ((i+1) % 16) == 0 )
        fprintf(log, "\n");
    }
    fprintf(log, "\n");
  }
 
  Slab_Rank = 2;
  Slab_Dims[0] = 3;
  Slab_Dims[1] = NumberOfParticles;
 
  if ( NumberOfParticles == 0 )
  {
    Slab_Dims[1] = 1;
  }
 
/*
  component_rank_attr = 3;
  component_size_attr = NumberOfParticles;
  field_rank_attr = 1;
  field_dims_attr = NumberOfParticles;
*/
 
  // Data in memory is considered 1D, stride 1, with zero offset
 
  mem_stride = 1;                   // contiguous elements
  mem_count = NumberOfParticles;    // number of elements in field
  mem_offset = 0;                   // zero offset in buffer
 
  if ( NumberOfParticles == 0 )
  {
    mem_count = 1;
  }
 
  // 1D memory model
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );
 
  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  // Data in the file is (1+Rank)D with Npart components per grid point.
  // Offset[0] is the component Part of Npart components.  Data for each
  // Part are contiguous in the file, so stride = 1.
 
  slab_stride[0] = 1;      // contiguous elements
  slab_count[0] = 1;       // one component per call
  slab_offset[0] = dim;    // component Part of Npart
 
  slab_stride[1] = 1;                   // contiguous elements
  slab_count[1] = NumberOfParticles;    // field dimensions
  slab_offset[1] = 0;                   // complete field, no offset
 
  if ( NumberOfParticles == 0 )
  {
    slab_count[1] = 1;
  }
 
  file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
    assert( file_dsp_id != h5_error );
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  if ( dim == 0 )
  {
    file_id = H5Fcreate(PPos, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fcreate with name = %s\n", PPos);
      assert( file_id != h5_error );
 
    dset_id =  H5Dcreate(file_id, PPos, file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Dcreate with Name = %s\n", PPos);
      assert( dset_id != h5_error );
 
    attr_count = 1;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "NumberOfParticles", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = NumberOfParticles\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &NumberOfParticles);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "TotalParticleCount", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = TotalParticleCount\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &TotalParticleCount);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_count = 3;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridLeft", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridLeft\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Left);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridRight", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridRight\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Right);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
  }
  else
  {
    file_id = H5Fopen(PPos, H5F_ACC_RDWR, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fopen with Name = %s\n", PPos);
      assert( file_id != h5_error );
 
    dset_id =  H5Dopen(file_id, PPos);
      if (io_log) fprintf(log, "H5Dopen with Name = %s\n", PPos);
      assert( dset_id != h5_error );
  }
 
  if ( NumberOfParticles > 0 )
  {
    h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, outbuff);
      if (io_log) fprintf(log, "H5Dwrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
  }
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
 
  } // end of loop over dims
 
  } 
 
 
//  Read particle velocities under bit mask

  if (mode > 0) {
 
  if (io_log) fprintf(log, "Read particle velocities under bit mask\n");
 
  for (dim = 0; dim < rank; dim++)
  {
 
  ppc = 0; // output particle counter
 
  dims[0] = 3;
  dims[1] = TotalParticleCount;
 
  mem_type_id = HDF5_R8;

#ifdef MPIO
 
  //  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);
  //    if (io_log) fprintf(log, "H5Pcreate file_access_template\n");
  //    assert( file_acc_template != h5_error );
 
  //  h5_status = H5Pset_fapl_mpio(file_acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);
  //    if (io_log) fprintf(log, "H5Pset_fapl_mpio: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#else
 
  file_acc_template = H5P_DEFAULT;
    if (io_log) fprintf(log, "Default file_access_template\n");

#endif
 
  file_id = H5Fopen(PVin, H5F_ACC_RDONLY, file_acc_template);
    if (io_log) fprintf(log, "H5Fopen with Name = %s\n", PVin);
    assert( file_id != h5_error );
 
  dset_id = H5Dopen(file_id, PVin);
    if (io_log) fprintf(log, "H5Dopen with Name = %s\n", PVin);
    assert( dset_id != h5_error );
 
  mem_offset = 0;
  mem_count = ParticlesPerCPU;
  mem_stride = 1;
 
  if ( jcpu == ncpu - 1 )
  {
    mem_count = ParticlesPerCPU - Fakes;
  }
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );
 
  h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  file_dsp_id = H5Screate_simple((Eint32) 2, dims, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
    assert( file_dsp_id != h5_error );
 
  slab_offset[0] = dim;   // x,y,z
  slab_stride[0] = 1;
  slab_count[0] = 1;
 
  slab_offset[1] = jcpu * ParticlesPerCPU;
  slab_stride[1] = 1;
  slab_count[1] = ParticlesPerCPU;
 
  if ( jcpu == ncpu - 1 )
  {
    slab_count[1] = ParticlesPerCPU - Fakes;
  }
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
  //    if (io_log) fprintf(log, "H5Pcreate xfer_prop_list\n");
  //    assert( xfer_prop_list != h5_error );
 
  //  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
  //    if (io_log) fprintf(log, "H5Pset_dxpl_mpio: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#else
 
  xfer_prop_list = H5P_DEFAULT;
   if (io_log) fprintf(log, "Default xfer_prop_list\n");

#endif
 
  h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, buff[0]);
    if (io_log) fprintf(log, "H5Dread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  h5_status = H5Pclose(xfer_prop_list);
  //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#endif
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  h5_status = H5Pclose(file_acc_template);
  //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#endif
 
/*
  if ( jcpu == ncpu - 1 )
  {
     for ( nfake = 0; nfake < Fakes; nfake++ )
     {
       buff[0][ParticlesPerCPU - Fakes + nfake] = 0.0;
     }
  }
*/
 
  MPI_Barrier(MPI_COMM_WORLD);
 
  thisnode = jcpu;
  nextnode = (thisnode+1) % ncpu;
  prevnode = (thisnode-1+ncpu) % ncpu;
 
  stype = 30000+thisnode;
  rtype = 30000+prevnode;
  ltype = 30000+nextnode;
 
  iblk = - ( (ncpu - jcpu) % ncpu );
 
  for (k = 0; k < 2*ncpu; k++) // two passes
  {
 
  i = (k  ) % 2;
  j = (k+1) % 2;
 
  MPI_Barrier(MPI_COMM_WORLD);
 
  // Transmit left
 
 
  Count = ParticlesPerCPU;
  Source = nextnode;
  Tag = ltype;
  Type = MPI_DOUBLE;
 
//  ier = MPI_Irecv( buff[j], ParticlesPerCPU, MPI_DOUBLE, nextnode, ltype, MPI_COMM_WORLD, req1 );
  ier = MPI_Irecv( buff[j], Count, Type, Source, Tag, MPI_COMM_WORLD, req1 );
  if (io_log) fprintf(log, "IRECV proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, j, ier);
 
  Count = ParticlesPerCPU;
  Dest = prevnode;
  Tag = stype;
  Type = MPI_DOUBLE;
 
  ier = MPI_Isend( buff[i], ParticlesPerCPU, MPI_DOUBLE, prevnode, stype, MPI_COMM_WORLD, req2 );
  if (io_log) fprintf(log, "SEND proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, i, ier);
 
  // work on buff[i]
 
  if (io_log) fprintf(log, "K = %"ISYM", IBLK = %"ISYM"\n", k, iblk);
 
  if ( iblk > -1 && iblk < ncpu )
  {
 
    ipc = (jcpu+k+ncpu) % ncpu;
 
    ipc = ipc * ParticlesPerCPU;
    jpc = ipc + ParticlesPerCPU - 1;
 
    if (io_log) fprintf(log, "  ipc = %"ISYM" to %"ISYM"\n", ipc, jpc);
 
    for (ic = 0; ic < ParticlesPerCPU; ic++)
    {
      MaskAddr = ipc;
      WordAddr = MaskAddr/BitsPerInt;
      BitAddr  = MaskAddr%BitsPerInt;
      TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
      if ( TestBit == 0 )
        XMask = 0;
      if ( TestBit != 0 )
      {
        XMask = 1;
        outbuff[ppc] = buff[i][ic];
        ppc++;
      }
      ipc++;
    }
 
  }
 
  iblk++;

  ier = MPI_Wait( req2, stat ); 
  ier = MPI_Wait( req1, stat );
 
  if (io_log) fprintf(log, "WAIT proc %"ISYM" err %"ISYM"\n", thisnode, ier);
 
  } // end loop over rings
 
  fprintf(log, "L %8.4"FSYM"  R %8.4"FSYM"\n", GridLeft[jcpu][dim], GridRight[jcpu][dim]);
 
  if (io_log_d)
  {
    for (i = 0; i < NumberOfParticles; i++)
    {
      fprintf(log, "%8.4"FSYM, outbuff[i]);
      if ( ((i+1) % 16) == 0 )
        fprintf(log, "\n");
    }
    fprintf(log, "\n");
  }
 
  Slab_Rank = 2;
  Slab_Dims[0] = 3;
  Slab_Dims[1] = NumberOfParticles;
 
  if ( NumberOfParticles == 0 )
  {
    Slab_Dims[1] = 1;
  }
 
/*
  component_rank_attr = 3;
  component_size_attr = NumberOfParticles;
  field_rank_attr = 1;
  field_dims_attr = NumberOfParticles;
*/
 
  // Data in memory is considered 1D, stride 1, with zero offset
 
  mem_stride = 1;                   // contiguous elements
  mem_count = NumberOfParticles;    // number of elements in field
  mem_offset = 0;                   // zero offset in buffer
 
  if ( NumberOfParticles == 0 )
  {
    mem_count = 1;
  }
 
  // 1D memory model
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );
 
  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  // Data in the file is (1+Rank)D with Npart components per grid point.
  // Offset[0] is the component Part of Npart components.  Data for each
  // Part are contiguous in the file, so stride = 1.
 
  slab_stride[0] = 1;      // contiguous elements
  slab_count[0] = 1;       // one component per call
  slab_offset[0] = dim;    // component Part of Npart
 
  slab_stride[1] = 1;                   // contiguous elements
  slab_count[1] = NumberOfParticles;    // field dimensions
  slab_offset[1] = 0;                   // complete field, no offset
 
  if ( NumberOfParticles == 0 )
  {
    slab_count[1] = 1;
  }
 
  file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
   if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
   assert( file_dsp_id != h5_error );
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  if ( dim == 0 )
  {
    file_id = H5Fcreate(PVel, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fcreate with name = %s\n", PVel);
      assert( file_id != h5_error );
 
    dset_id =  H5Dcreate(file_id, PVel, file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Dcreate with Name = %s\n", PVel);
      assert( dset_id != h5_error );
 
    attr_count = 1;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "NumberOfParticles", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = NumberOfParticles\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &NumberOfParticles);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "TotalParticleCount", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = TotalParticleCount\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &TotalParticleCount);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_count = 3;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridLeft", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridLeft\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Left);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridRight", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridRight\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Right);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
  }
  else
  {
    file_id = H5Fopen(PVel, H5F_ACC_RDWR, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fopen with Name = %s\n", PVel);
      assert( file_id != h5_error );
 
    dset_id =  H5Dopen(file_id, PVel);
      if (io_log) fprintf(log, "H5Dopen with Name = %s\n", PPos);
      assert( dset_id != h5_error );
 
  }
 
  if ( NumberOfParticles > 0 )
  {
    h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, outbuff);
      if (io_log) fprintf(log, "H5Dwrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
  }
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
 
  } // end of loop over dims

  }
 
// CHOP  }
 
 
//  Read particle masses under bit mask

  if ((mode == 2) || (mode == 4)) {
 
  if (io_log) fprintf(log, "Read particle masses under bit mask\n");
 
  dim = 0;
 
  ppc = 0; // output particle counter
 
  dims[0] = 1;
  dims[1] = TotalParticleCount;
 
  mem_type_id = HDF5_R8;
  file_type_id = HDF5_FILE_R8;

#ifdef MPIO
 
  //  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);
  //    if (io_log) fprintf(log, "H5Pcreate file_access_template\n");
  //    assert( file_acc_template != h5_error );
 
  //  h5_status = H5Pset_fapl_mpio(file_acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);
  //    if (io_log) fprintf(log, "H5Pset_fapl_mpio: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#else
 
  file_acc_template = H5P_DEFAULT;
    if (io_log) fprintf(log, "Default file_access_template\n");

#endif
 
  file_id = H5Fopen(PMin, H5F_ACC_RDONLY, file_acc_template);
    if (io_log) fprintf(log, "H5Fopen with Name = %s\n", PMin);
    assert( file_id != h5_error );
 
  dset_id = H5Dopen(file_id, PMin);
    if (io_log) fprintf(log, "H5Dopen with Name = %s\n", PMin);
    assert( dset_id != h5_error );
 
  mem_offset = 0;
  mem_count = ParticlesPerCPU;
  mem_stride = 1;
 
  if ( jcpu == ncpu - 1 )
  {
    mem_count = ParticlesPerCPU - Fakes;
  }
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );
 
  h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  file_dsp_id = H5Screate_simple((Eint32) 2, dims, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
    assert( file_dsp_id != h5_error );
 
  slab_offset[0] = dim;   // x,y,z
  slab_stride[0] = 1;
  slab_count[0] = 1;
 
  slab_offset[1] = jcpu * ParticlesPerCPU;
  slab_stride[1] = 1;
  slab_count[1] = ParticlesPerCPU;
 
  if ( jcpu == ncpu - 1 )
  {
    slab_count[1] = ParticlesPerCPU - Fakes;
  }
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
  //    if (io_log) fprintf(log, "H5Pcreate xfer_prop_list\n");
  //    assert( xfer_prop_list != h5_error );
 
  //  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
  //    if (io_log) fprintf(log, "H5Pset_dxpl_mpio: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#else
 
  xfer_prop_list = H5P_DEFAULT;
   if (io_log) fprintf(log, "Default xfer_prop_list\n");

#endif
 
  h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, buff[0]);
    if (io_log) fprintf(log, "H5Dread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  h5_status = H5Pclose(xfer_prop_list);
  //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#endif
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  h5_status = H5Pclose(file_acc_template);
  //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#endif
 
/*
  if ( jcpu == ncpu - 1 )
  {
     for ( nfake = 0; nfake < Fakes; nfake++ )
     {
       buff[0][ParticlesPerCPU - Fakes + nfake] = 0.0;
     }
  }
*/
 
  MPI_Barrier(MPI_COMM_WORLD);
 
  thisnode = jcpu;
  nextnode = (thisnode+1) % ncpu;
  prevnode = (thisnode-1+ncpu) % ncpu;
 
  stype = 30000+thisnode;
  rtype = 30000+prevnode;
  ltype = 30000+nextnode;
 
  iblk = - ( (ncpu - jcpu) % ncpu );
 
  for (k = 0; k < 2*ncpu; k++) // two passes
  {
 
  i = (k  ) % 2;
  j = (k+1) % 2;
 
  MPI_Barrier(MPI_COMM_WORLD);
 
  // Transmit left
 
  Count = ParticlesPerCPU;
  Source = nextnode;
  Tag = ltype;
  Type = MPI_DOUBLE;
 
//  ier = MPI_Irecv( buff[j], ParticlesPerCPU, MPI_DOUBLE, nextnode, ltype, MPI_COMM_WORLD, req1 );
  ier = MPI_Irecv( buff[j], Count, Type, Source, Tag, MPI_COMM_WORLD, req1 );
  if (io_log) fprintf(log, "IRECV proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, j, ier);
 
  Count = ParticlesPerCPU;
  Dest = prevnode;
  Tag = stype;
  Type = MPI_DOUBLE;
 
//  ier = MPI_Isend( buff[i], ParticlesPerCPU, MPI_DOUBLE, prevnode, stype, MPI_COMM_WORLD, req2 );
  ier = MPI_Isend( buff[i], Count, Type, Dest, Tag, MPI_COMM_WORLD, req2 );
  if (io_log) fprintf(log, "SEND proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, i, ier);
 
  // work on buff[i]
 
  if (io_log) fprintf(log, "K = %"ISYM", IBLK = %"ISYM"\n", k, iblk);
 
  if ( iblk > -1 && iblk < ncpu )
  {
 
    ipc = (jcpu+k+ncpu) % ncpu;
 
    ipc = ipc * ParticlesPerCPU;
    jpc = ipc + ParticlesPerCPU - 1;
 
    if (io_log) fprintf(log, "  ipc = %"ISYM" to %"ISYM"\n", ipc, jpc);
 
    for (ic = 0; ic < ParticlesPerCPU; ic++)
    {
      MaskAddr = ipc;
      WordAddr = MaskAddr/BitsPerInt;
      BitAddr  = MaskAddr%BitsPerInt;
      TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
      if ( TestBit == 0 )
        XMask = 0;
      if ( TestBit != 0 )
      {
        XMask = 1;
        outbuff[ppc] = buff[i][ic];
        ppc++;
      }
      ipc++;
    }
 
  }
 
  iblk++;

  ier = MPI_Wait( req2, stat ); 
  ier = MPI_Wait( req1, stat );
 
  if (io_log) fprintf(log, "WAIT proc %"ISYM" err %"ISYM"\n", thisnode, ier);
 
  } // end loop over rings
 
  fprintf(log, "L %8.4"FSYM"  R %8.4"FSYM"\n", GridLeft[jcpu][dim], GridRight[jcpu][dim]);
 
  if (io_log_d)
  {
    for (i = 0; i < NumberOfParticles; i++)
    {
      fprintf(log, "%8.4"FSYM, outbuff[i]);
      if ( ((i+1) % 16) == 0 )
        fprintf(log, "\n");
    }
    fprintf(log, "\n");
  }
 
  Slab_Rank = 2;
  Slab_Dims[0] = 1;
  Slab_Dims[1] = NumberOfParticles;
 
  if ( NumberOfParticles == 0 )
  {
    Slab_Dims[1] = 1;
  }
 
  // Data in memory is considered 1D, stride 1, with zero offset
 
  mem_stride = 1;                   // contiguous elements
  mem_count = NumberOfParticles;    // number of elements in field
  mem_offset = 0;                   // zero offset in buffer
 
  if ( NumberOfParticles == 0 )
  {
    mem_count = 1;
  }
 
  // 1D memory model
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );
 
  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  // Data in the file is (1+Rank)D with Npart components per grid point.
  // Offset[0] is the component Part of Npart components.  Data for each
  // Part are contiguous in the file, so stride = 1.
 
  slab_stride[0] = 1;      // contiguous elements
  slab_count[0] = 1;       // one component per call
  slab_offset[0] = dim;    // component Part of Npart
 
  slab_stride[1] = 1;                   // contiguous elements
  slab_count[1] = NumberOfParticles;    // field dimensions
  slab_offset[1] = 0;                   // complete field, no offset
 
  if ( NumberOfParticles == 0 )
  {
    slab_count[1] = 1;
  }
 
  file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
   if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
   assert( file_dsp_id != h5_error );
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  if ( dim == 0 )
  {
    file_id = H5Fcreate(PMass, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fcreate with name = %s\n", PMass);
      assert( file_id != h5_error );
 
    dset_id =  H5Dcreate(file_id, PMass, file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Dcreate with Name = %s\n", PMass);
      assert( dset_id != h5_error );
 
    attr_count = 1;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "NumberOfParticles", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = NumberOfParticles\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &NumberOfParticles);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "TotalParticleCount", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = TotalParticleCount\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &TotalParticleCount);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_count = 3;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridLeft", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridLeft\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Left);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridRight", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridRight\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Right);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
  }
  else
  {
    file_id = H5Fopen(PMass, H5F_ACC_RDWR, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fopen with Name = %s\n", PMass);
      assert( file_id != h5_error );
 
    dset_id =  H5Dopen(file_id, PMass);
      if (io_log) fprintf(log, "H5Dopen with Name = %s\n", PMass);
      assert( dset_id != h5_error );
 
  }
 
  if ( NumberOfParticles > 0 )
  {
    h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, outbuff);
      if (io_log) fprintf(log, "H5Dwrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
  }
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  } 
 
 
//  Read particle types under bit mask

  if ((mode == 3) || (mode == 4)) {
 
  if (io_log) fprintf(log, "Read particle types under bit mask\n");
 
  int *ibuff[2];
  ibuff[0] = new int[dbuff_size];
  ibuff[1] = new int[dbuff_size];
 
  int *outibuff = NULL;
  outibuff = new int[NumberOfParticles];
 
  dim = 0;
 
  ppc = 0; // output particle counter
 
  dims[0] = 1;
  dims[1] = TotalParticleCount;
 
  int_mem_type_id = HDF5_INT;
  int_file_type_id = HDF5_FILE_INT;

#ifdef MPIO
 
  //  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);
  //    if (io_log) fprintf(log, "H5Pcreate file_access_template\n");
  //    assert( file_acc_template != h5_error );
 
  //  h5_status = H5Pset_fapl_mpio(file_acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);
  //    if (io_log) fprintf(log, "H5Pset_fapl_mpio: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#else
 
  file_acc_template = H5P_DEFAULT;
    if (io_log) fprintf(log, "Default file_access_template\n");

#endif
 
  file_id = H5Fopen(PTin, H5F_ACC_RDONLY, file_acc_template);
    if (io_log) fprintf(log, "H5Fopen with Name = %s\n", PTin);
    assert( file_id != h5_error );
 
  dset_id = H5Dopen(file_id, PTin);
    if (io_log) fprintf(log, "H5Dopen with Name = %s\n", PTin);
    assert( dset_id != h5_error );
 
  mem_offset = 0;
  mem_count = ParticlesPerCPU;
  mem_stride = 1;
 
  if ( jcpu == ncpu - 1 )
  {
    mem_count = ParticlesPerCPU - Fakes;
  }
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );
 
  h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  file_dsp_id = H5Screate_simple((Eint32) 2, dims, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
    assert( file_dsp_id != h5_error );
 
  slab_offset[0] = dim;   // x,y,z
  slab_stride[0] = 1;
  slab_count[0] = 1;
 
  slab_offset[1] = jcpu * ParticlesPerCPU;
  slab_stride[1] = 1;
  slab_count[1] = ParticlesPerCPU;
 
  if ( jcpu == ncpu - 1 )
  {
    slab_count[1] = ParticlesPerCPU - Fakes;
  }
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
  //    if (io_log) fprintf(log, "H5Pcreate xfer_prop_list\n");
  //    assert( xfer_prop_list != h5_error );
 
  //  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
  //    if (io_log) fprintf(log, "H5Pset_dxpl_mpio: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#else
 
  xfer_prop_list = H5P_DEFAULT;
   if (io_log) fprintf(log, "Default xfer_prop_list\n");

#endif
 
  h5_status = H5Dread(dset_id, int_mem_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, ibuff[0]);
    if (io_log) fprintf(log, "H5Dread: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  h5_status = H5Pclose(xfer_prop_list);
  //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#endif
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

#ifdef MPIO
 
  //  h5_status = H5Pclose(file_acc_template);
  //    if (io_log) fprintf(log, "H5Pclose: %"ISYM"\n", h5_status);
  //    assert( h5_status != h5_error );

#endif
 
/*
  if ( jcpu == ncpu - 1 )
  {
     for ( nfake = 0; nfake < Fakes; nfake++ )
     {
       ibuff[0][ParticlesPerCPU - Fakes + nfake] = 0;
     }
  }
*/
 
  MPI_Barrier(MPI_COMM_WORLD);
 
  thisnode = jcpu;
  nextnode = (thisnode+1) % ncpu;
  prevnode = (thisnode-1+ncpu) % ncpu;
 
  stype = 30000+thisnode;
  rtype = 30000+prevnode;
  ltype = 30000+nextnode;
 
  iblk = - ( (ncpu - jcpu) % ncpu );
 
  for (k = 0; k < 2*ncpu; k++) // two passes
  {
 
  i = (k  ) % 2;
  j = (k+1) % 2;
 
  MPI_Barrier(MPI_COMM_WORLD);
 
  // Transmit left
 
  Count = ParticlesPerCPU;
  Source = nextnode;
  Tag = ltype;
  Type = IntDataType;
 
//  ier = MPI_Irecv( ibuff[j], ParticlesPerCPU, MPI_INT, nextnode, ltype, MPI_COMM_WORLD, req1 );
  ier = MPI_Irecv( ibuff[j], Count, Type, Source, Tag, MPI_COMM_WORLD, req1 );
  if (io_log) fprintf(log, "IRECV proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, j, ier);
 
  Count = ParticlesPerCPU;
  Dest = prevnode;
  Tag = stype;
  Type = IntDataType;
 
//  ier = MPI_Isend( ibuff[i], ParticlesPerCPU, MPI_INT, prevnode, stype, MPI_COMM_WORLD, req2 );
  ier = MPI_Isend( ibuff[i], Count, Type, Dest, Tag, MPI_COMM_WORLD, req2 );
  if (io_log) fprintf(log, "SEND proc %"ISYM" buffer[%"ISYM"] err %"ISYM"\n", thisnode, i, ier);
 
  // work on ibuff[i]
 
  if (io_log) fprintf(log, "K = %"ISYM", IBLK = %"ISYM"\n", k, iblk);
 
  if ( iblk > -1 && iblk < ncpu )
  {
 
    ipc = (jcpu+k+ncpu) % ncpu;
 
    ipc = ipc * ParticlesPerCPU;
    jpc = ipc + ParticlesPerCPU - 1;
 
    if (io_log) fprintf(log, "  ipc = %"ISYM" to %"ISYM"\n", ipc, jpc);
 
    for (ic = 0; ic < ParticlesPerCPU; ic++)
    {
      MaskAddr = ipc;
      WordAddr = MaskAddr/BitsPerInt;
      BitAddr  = MaskAddr%BitsPerInt;
      TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
      if ( TestBit == 0 )
        XMask = 0;
      if ( TestBit != 0 )
      {
        XMask = 1;
        outibuff[ppc] = ibuff[i][ic];
        ppc++;
      }
      ipc++;
    }
 
  }
 
  iblk++;

  ier = MPI_Wait( req2, stat ); 
  ier = MPI_Wait( req1, stat );
 
  if (io_log) fprintf(log, "WAIT proc %"ISYM" err %"ISYM"\n", thisnode, ier);
 
  } // end loop over rings
 
  fprintf(log, "L %8.4"FSYM"  R %8.4"FSYM"\n", GridLeft[jcpu][dim], GridRight[jcpu][dim]);
 
  if (io_log_d)
  {
    for (i = 0; i < NumberOfParticles; i++)
    {
      fprintf(log, "%8"ISYM, outibuff[i]);
      if ( ((i+1) % 16) == 0 )
        fprintf(log, "\n");
    }
    fprintf(log, "\n");
  }
 
  Slab_Rank = 2;
  Slab_Dims[0] = 1;
  Slab_Dims[1] = NumberOfParticles;
 
  if ( NumberOfParticles == 0 )
  {
    Slab_Dims[1] = 1;
  }
 
/*
  component_rank_attr = 3;
  component_size_attr = NumberOfParticles;
  field_rank_attr = 1;
  field_dims_attr = NumberOfParticles;
*/
 
  // Data in memory is considered 1D, stride 1, with zero offset
 
  mem_stride = 1;                   // contiguous elements
  mem_count = NumberOfParticles;    // number of elements in field
  mem_offset = 0;                   // zero offset in buffer
 
  if ( NumberOfParticles == 0 )
  {
    mem_count = 1;
  }
 
  // 1D memory model
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );
 
  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  // Data in the file is (1+Rank)D with Npart components per grid point.
  // Offset[0] is the component Part of Npart components.  Data for each
  // Part are contiguous in the file, so stride = 1.
 
  slab_stride[0] = 1;      // contiguous elements
  slab_count[0] = 1;       // one component per call
  slab_offset[0] = dim;    // component Part of Npart
 
  slab_stride[1] = 1;                   // contiguous elements
  slab_count[1] = NumberOfParticles;    // field dimensions
  slab_offset[1] = 0;                   // complete field, no offset
 
  if ( NumberOfParticles == 0 )
  {
    slab_count[1] = 1;
  }
 
  file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
   if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
   assert( file_dsp_id != h5_error );
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
    if (io_log) fprintf(log, "H5Sselect_hyperslab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  if ( dim == 0 )
  {
    file_id = H5Fcreate(PType, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fcreate with name = %s\n", PType);
      assert( file_id != h5_error );
 
    dset_id =  H5Dcreate(file_id, PType, int_file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Dcreate with Name = %s\n", PType);
      assert( dset_id != h5_error );
 
    attr_count = 1;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "NumberOfParticles", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = NumberOfParticles\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &NumberOfParticles);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "TotalParticleCount", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = TotalParticleCount\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_INT, &TotalParticleCount);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_count = 3;
 
    attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridLeft", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridLeft\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Left);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_id = H5Acreate(dset_id, "GridRight", HDF5_FILE_R8, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate with Name = GridRight\n");
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, HDF5_R8, Right);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
  }
  else
  {
    file_id = H5Fopen(PType, H5F_ACC_RDWR, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fopen with Name = %s\n", PType);
      assert( file_id != h5_error );
 
    dset_id =  H5Dopen(file_id, PType);
      if (io_log) fprintf(log, "H5Dopen with Name = %s\n", PType);
      assert( dset_id != h5_error );
 
  }
 
  if ( NumberOfParticles > 0 )
  {
    h5_status = H5Dwrite(dset_id, int_mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, outibuff);
      if (io_log) fprintf(log, "H5Dwrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
  }
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  }
 
 
  fclose(log);

  MPI_Barrier(MPI_COMM_WORLD);

  if ( mpi_rank == 0 )
  {
    fprintf(stdout, "Sort completed\n");
  }
 
  MPI_Finalize();
 
}
