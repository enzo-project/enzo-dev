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
      ngrids = 0;
      ncells = 0;
      total_time = 0.0;
      current_time = 0.0;
    }
    
    std::string name;

    void start(void){
      t0 = ReturnWallTime();
    }

    void stop(void){
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

    long int get_cells(void){
      return ncells;
    }
    long int get_grids(void){
      return ngrids;
    }

    void add_stats(long int my_grids, long int my_cells){
      ngrids = my_grids;
      ncells = my_cells;
    }

    section_performance *next;
  
  private:
    double t0;
    double t1;
    double total_time;
    double current_time;
    long int ngrids;
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
      timers.clear();
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

    section_performance * get(char *name){
      iter = timers.find(name);
      if (iter == timers.end()){
	fprintf(stderr, "%s Being Created\n", name);
	timers[name] = new section_performance(name);
      }
      return timers[name];
    }

    section_performance * get_level(int level){
      char level_name[256];
      sprintf(level_name, "Level_%d", level);
      return get(level_name);
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
	if (step == 1){
	  fprintf(chronos_file, "# This file contains timing information\n");
	  fprintf(chronos_file, "# For instructions on how to decipher this information,\n");
	  fprintf(chronos_file, "# see http://enzo-project.org/docs/somewhere\n");
	  fprintf(chronos_file, "\n");
	}
      }
      


      double *time_array;
      if (my_rank == 0){
	time_array = new double[nprocs];
      }
      std::string keyname;

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
	fprintf(chronos_file, "Cycle_Number %d\n",step);
      }
      
      long int total_cells = 0;
      
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
	  stddev_time = sqrt(stddev_time/nprocs);

      	  fprintf(chronos_file, "%s %e %e %e %e",
      		  iter->first.c_str(), mean_time, stddev_time, min_time, max_time);
	  keyname = iter->first;
	  if (std::strncmp(keyname.c_str(), "Total", 5) == 0){
	    total_time = mean_time;
	  }
	  if (std::strncmp(keyname.c_str(), "Level", 5) == 0){
	    total_cells += iter->second->get_cells();
	    fprintf(chronos_file, "% ld %ld %e", 
		    iter->second->get_cells(),
		    iter->second->get_grids(),
		    (double)(iter->second->get_cells()/mean_time/nprocs));
	  }
      	  if(verbose){
      	    for (int i=0; i<nprocs; i++){
      	      fprintf(chronos_file, " %e", time_array[i]);
      	    }
      	  }
      	  fprintf(chronos_file, "\n");
      	}
      }

      if (my_rank == 0){
	fprintf(chronos_file, "Total_Cells/Sec/Processor: %e\n", (double)(total_cells/total_time/nprocs));
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

EXTERN enzo_timing::enzo_timer *enzo_timer;                                                         

#define TIMER_START(section_name) enzo_timer->get(section_name)->start()
#define TIMER_STOP(section_name) enzo_timer->get(section_name)->stop()

#endif //ENZO_TIMING
