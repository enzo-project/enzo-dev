#ifndef ENZO_TIMING__
#define ENZO_TIMING__
#include <stdio.h>
#include <math.h>
#include <string>

#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))

double ReturnWallTime(void);
void Reduce_Times(double time, double *time_array);

namespace enzo_timing{

  class section_performance
  {
  public:
    
    section_performance(string myname){
      name = myname;
      next = NULL;
      total_time = 0.0;
      current_time = 0.0;
    }
    
    string myname;

    void start_timer(void){
      t0 = ReturnWallTime();
    }

    void stop_and_add_timer(void){
      t1 = ReturnWallTime();
      current_time += t1-t0;
      total_time += t1-t0;
    }

    double get_total_time(void){
      return total_time;
    }
    double get_current_time(void){
      return current_time;
    }
    void reset_current_time(void){
      current_time = 0.0;  
    }

    section_performance *next;
  
  private:
    double t0;
    double t1; 
    double total_time;
    double current_time;
    int ngrids;
    long int ncells;
    
  };

 

  class level_performance
  {
  public:
    
    level_performance(int this_level){
      level = this_level;
      next = NULL;
      total_time = 0.0;
      current_time = 0.0;
    }
    
    int level;

    void start_timer(void){
      t0 = ReturnWallTime();
    }

    void stop_and_add_timer(void){
      t1 = ReturnWallTime();
      current_time += t1-t0;
      total_time += t1-t0;
    }

    template <typename T>
    void set_ngrids(T grids){
      ngrids = (int) (grids);
    }
    template <typename T>
    void set_ncells(T cells){
      ngrids = (long int)(cells);
    }
    
    int get_ngrids(void){
      return (int)ngrids;
    }
    long int get_ncells(void){
      return ncells;
    }

    double get_total_time(void){
      return total_time;
    }
    double get_current_time(void){
      return current_time;
    }
    void reset_current_time(void){
      current_time = 0.0;  
    }

    level_performance *next;
  
  private:
    double t0;
    double t1; 
    double total_time;
    double current_time;
    int ngrids;
    long int ncells;
    
  };

  /* --------------------------------------------------------- */

  class enzo_timer
  {
  public:

    enzo_timer(void){
      total_time = 0.0;
      current_time = 0.0;
      level_zero_performance = new level_performance(0);
      filename = (char *)("chronos.out");
      set_mpi_environment();
    }

    enzo_timer(char *performance_name){
      total_time = 0.0;
      current_time = 0.0;
      level_zero_performance = new level_performance(0);
      filename = performance_name;
      set_mpi_environment();
    }
    
    void set_mpi_environment(void){
#ifdef USE_MPI
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
      nprocs = 1;
      my_rank = 0;
#endif
    }

    level_performance *level_zero_performance;
    level_performance *temp_performance;

    double get_total_time(void){
      temp_performance = level_zero_performance;
      while(temp_performance != NULL){
	total_time += temp_performance->get_total_time();
	temp_performance = temp_performance->next;
      }
      return total_time;
    }

    int get_total_grids(void){
      temp_performance = level_zero_performance;
      int total_grids = 0;
      while(temp_performance != NULL){
	total_grids += temp_performance->get_ngrids();
	temp_performance = temp_performance->next;
      }
      
      return total_grids;
    }

    double get_current_time(void){
      current_time = 0.0;
      temp_performance = level_zero_performance;
      while(temp_performance != NULL){
	current_time += temp_performance->get_current_time();
	temp_performance = temp_performance->next;
      }
      return current_time;
    }

    level_performance* get_level_performance(int level){
      temp_performance = level_zero_performance;
      while(temp_performance != NULL){
	if (temp_performance->level == level){
	  break;
	}
	if (temp_performance->next == NULL){
	  temp_performance->next = new level_performance(level);
	  fprintf(stderr,"Creating Performance Level %d\n",level);
	}
	temp_performance = temp_performance->next;
      }
      if (temp_performance == NULL){
	fprintf(stderr,"Couldn't find the performance level!\n");
      }
      return temp_performance;
    }

    void write_out(int step, bool verbose=false){
      if (my_rank == 0){
	chronos_file = fopen(filename,"a");
      }
      double *time_array;
      if (my_rank == 0){
	time_array = new double[nprocs];
      }

      double thistime = this->get_current_time();
      Reduce_Times(thistime, time_array);

      double mean_time = 0.0;
      double min_time, max_time, stddev_time;
      min_time = max_time = time_array[0];
      for (int i=0; i<nprocs; i++){
	mean_time += time_array[i];
	min_time = min(min_time, time_array[i]);
	max_time = max(max_time, time_array[i]);
      }
      mean_time /= nprocs;

      stddev_time = 0.0;
      for (int i=0; i<nprocs; i++){
	stddev_time += pow((double)(time_array[i]-mean_time),(double)(2.0));
      }
      stddev_time = sqrt(stddev_time/nprocs);
      
      if (my_rank == 0){
	fprintf(chronos_file, "Step %d %d %f %f %e %e",
		step, this->get_total_grids(), mean_time, stddev_time, min_time, max_time);
	if(verbose){
	  for (int i=0; i<nprocs; i++){
	    fprintf(chronos_file, " %e", time_array[i]);
	  }
	}
	fprintf(chronos_file, "\n");
      }

      temp_performance = level_zero_performance;
      while(temp_performance != NULL){
	current_time = temp_performance->get_current_time();
	temp_performance->reset_current_time();
	Reduce_Times(current_time, time_array);
	
	if (my_rank == 0){
	  mean_time = 0.0;
	  min_time = max_time = time_array[0];
	  for (int i=0; i<nprocs; i++){
	    mean_time += time_array[i];
	    min_time = min(min_time, time_array[i]);
	    max_time = max(max_time, time_array[i]);
	  }
	  mean_time /= nprocs;
	  
	  stddev_time = 0.0;
	  for (int i=0; i<nprocs; i++){
	    stddev_time += pow(time_array[i]-mean_time,2.0);
	  }
	  stddev_time = sqrt(stddev_time/nprocs);
	  fprintf(chronos_file, "Level %d %d %e %e %e %e",
		  temp_performance->level, temp_performance->get_ngrids(), mean_time, stddev_time, min_time, max_time);
	  if(verbose){
	    for (int i=0; i<nprocs; i++){
	      fprintf(chronos_file, " %e", time_array[i]);
	    }
	  }	  
	  fprintf(chronos_file, "\n");
	}
	temp_performance = temp_performance->next;
      }

      if (my_rank == 0){
	fprintf(chronos_file, "\n");
	fclose(chronos_file);	  
	delete [] time_array;
      }
    }
	      
  private:
    FILE * chronos_file;
    double total_time;
    double current_time;
    char * filename;
    int my_rank;
    int nprocs;
  };
}

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif
EXTERN enzo_timing::enzo_timer *my_enzo_timer;                                                         
EXTERN enzo_timing::enzo_timer *hydro_timer;  
EXTERN enzo_timing::enzo_timer *section_timer;

#endif //ENZO_TIMING
