#define _XOPEN_SOURCE 500
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
//#include <mpi.h>
//#define OFFSET_MAX 274877906944 // 4*4096^3
#define OFFSET_MAX 2199023255552 // 4*8192^3
int parallel_write(char*,long,long,void*);
int parallel_read(char*,long,long,void*);
void f77_parallel_read_(char*,int*,long*,long*,void*);
void f77_parallel_write_(char*,int*,long*,long*,void*);


// This routine writes size bytes of buffer at position offset
// in the file filename.
// Can be used in multiple concurrent access

int parallel_write(char *filename, long size, long offset, void *buffer) {
  
 
  int fd; //file descriptor
  int stat; // I/O status
  if (size + offset > OFFSET_MAX) {
    fprintf(stderr,"You are trying to access a file location\n");
    fprintf(stderr,"which is bigger than %ld\n",(long)OFFSET_MAX);
    fprintf(stderr,"Verify your code and/or change OFFSET_MAX\n");
    return(1);
  }
  fd = open(filename,O_RDWR|O_CREAT,S_IRUSR|S_IWUSR);
  stat=pwrite(fd,buffer,size,offset);
  close(fd);
  if (stat != size) return(1);
  else return(0);
}

// This routine reads size bytes of buffer at position offset
// in the file filename.
// Can be used in multiple concurrent access

int parallel_read(char *filename, long size, long offset, void *buffer) {

  int fd; //file descriptor
  int stat; // I/O status
  if (size + offset > OFFSET_MAX) {
    fprintf(stderr,"You are trying to access a file location\n");
    fprintf(stderr,"which is bigger than %ld\n",(long)OFFSET_MAX);
    fprintf(stderr,"Verify your code and/or change OFFSET_MAX\n");
    return(1);
  }
  fd = open(filename,O_RDONLY);
  stat = pread(fd,buffer,size,offset);
  close(fd);
  if (stat != size) return(1);
  else return(0);
}

// Fortran wrappers

void f77_parallel_read_(char *filename, int *fnamelen, long *size,
			long *offset, void *buffer) {

  int stat;
  filename[*fnamelen]='\0';
  stat = parallel_read(filename,*size,*offset,buffer);
  if (stat !=0) fprintf(stderr,"Failure to read with parallel_read\n");
}
void f77_parallel_write_(char *filename, int *fnamelen, long *size,
			 long *offset, void *buffer) {

  int stat;
  filename[*fnamelen]='\0';
  stat = parallel_write(filename,*size,*offset,buffer);
  if (stat !=0) fprintf(stderr,"Failure to write with parallel_read\n");
}



///
/*
TMP_main(int argc, char **argv) {

  int fd;
  char filename[128];
  long nsize, noffset;
  float *buffer;
  int i, stat;
  int myid, nproc;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if (argc != 3 && myid==0) {
    fprintf(stderr,"io filename nsize\n");
    MPI_Finalize();
    exit(1);
  }
    

  strcpy(filename,argv[1]);
  nsize=atoi(argv[2]); // number of floats to read/write
  // noffset=atoi(argv[3]); // offset in sizeof(float)
  noffset=myid*nsize;
  
  
  fd = open(filename,O_RDWR|O_CREAT,S_IRUSR|S_IWUSR);
  buffer = (float*)malloc(nsize*sizeof(float));
  for (i=0;i<nsize;i++) buffer[i] = (float)(i+nsize*myid);
  stat=pwrite(fd,(void*)buffer,nsize*sizeof(float),noffset*sizeof(float));
  close(fd);
  fprintf(stderr,"coucou from proc # %d: stat = %d\n",myid,stat);
  MPI_Barrier(MPI_COMM_WORLD);
  return;
  

  
  fd = open(filename,O_RDONLY);
  buffer = (float*)malloc(nsize*sizeof(float));
  stat=pread(fd,(void*)buffer,nsize*sizeof(float),noffset*sizeof(float));
  for (i=0;i<nsize;i++) fprintf(stderr,"%d %f\n",i+nsize*myid,buffer[i]);
  close(fd);
  return;
  
  MPI_Finalize();

}
*/
