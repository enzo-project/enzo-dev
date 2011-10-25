#ifndef ENZO_TIMING__
#define ENZO_TIMING__
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <map>

#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))

double ReturnWallTime(void);
void Reduce_Times(double time, double *time_array);

namespace enzo_timing{
  
  class section_performance
  {
  public:
    
    section_performance(char* myname){
      name = myname;
      next = NULL;
      total_time = 0.0;
      current_time = 0.0;
    }
    
    std::string name;

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


  /* --------------------------------------------------------- */

  class enzo_timer
  {
  public:

    enzo_timer(void){
      total_time = 0.0;
      current_time = 0.0;
      filename = (char *)("chronos.out");
      set_mpi_environment();
    }

    enzo_timer(char *performance_name){
      total_time = 0.0;
      current_time = 0.0;
      filename = performance_name;
      set_mpi_environment();
    }

    ~enzo_timer(void){
      for( iter=timers.begin(); iter!=timers.end(); ++iter)	{
	delete [] iter->second;
	timers.erase(iter);
      }
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

    std::map<std::string, section_performance *> timers;
    std::map<std::string, section_performance *>::iterator iter;

    section_performance * get_or_add_new(char *name){
      iter = timers.find(name);
      if (iter == timers.end()){
	fprintf(stderr, "%s Being Created\n", name);
	timers[name] = new section_performance(name);
      }
      return timers[name];
    }

    section_performance * get(char *name){
      return timers[name];
    }
    
    double get_total_time(void){
      total_time = 0.0;
      for( iter=timers.begin(); iter!=timers.end(); ++iter){
	total_time += iter->second->get_total_time();
      }
      return total_time;
    }

    double get_total_levels_time(void){
      total_time = 0.0;
      std::string keyname;
      for( iter=timers.begin(); iter!=timers.end(); ++iter){
	keyname = iter->first;
	if (std::strncmp(keyname.c_str(), "Level", 5) == 0){
	    total_time += iter->second->get_total_time();
	  }
	}
      return total_time;
    }

    int get_total_grids(void){
      return 0;
      /* temp_performance = level_zero_performance; */
      /* int total_grids = 0; */
      /* while(temp_performance != NULL){ */
      /* 	total_grids += temp_performance->get_ngrids(); */
      /* 	temp_performance = temp_performance->next; */
      /* } */
      
      /* return total_grids; */
    }

    double get_current_time(void){
      current_time = 0.0;
      for( iter=timers.begin(); iter!=timers.end(); ++iter){
	current_time += iter->second->get_current_time();
      }
      return current_time;
    }

    double get_current_levels_time(void){
      total_time = 0.0;
      std::string keyname;
      for( iter=timers.begin(); iter!=timers.end(); ++iter){
	keyname = iter->first;
	if (std::strncmp(keyname.c_str(), "Level", 5) == 0){
	    total_time += iter->second->get_current_time();
	  }
	}
      return total_time;
    }

    void write_out(int step, bool verbose=false){
      if (my_rank == 0){
	chronos_file = fopen(filename,"a");
      }
      double *time_array;
      if (my_rank == 0){
	time_array = new double[nprocs];
      }

      double thistime = this->get_current_levels_time();
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
	fprintf(chronos_file, "CycleNumber %d %d %f %f %e %e",
		step, this->get_total_grids(), mean_time, stddev_time, min_time, max_time);
	if(verbose){
	  for (int i=0; i<nprocs; i++){
	    fprintf(chronos_file, " %e", time_array[i]);
	  }
	}
	fprintf(chronos_file, "\n");
      }

      for( iter=timers.begin(); iter!=timers.end(); ++iter){
      	current_time = iter->second->get_current_time();
      	iter->second->reset_current_time();
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
      	  fprintf(chronos_file, "%s %e %e %e %e",
      		  iter->first.c_str(), mean_time, stddev_time, min_time, max_time);
      	  if(verbose){
      	    for (int i=0; i<nprocs; i++){
      	      fprintf(chronos_file, " %e", time_array[i]);
      	    }
      	  }
      	  fprintf(chronos_file, "\n");
      	}
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
/* EXTERN enzo_timing::enzo_timer *hydro_timer;   */
/* EXTERN enzo_timing::enzo_timer *section_timer; */

#endif //ENZO_TIMING
