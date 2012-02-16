/***********************************************************************
/
/  EnzoTiming.h  
/
/  written by: Samuel Skillman
/  date:       February, 2012
/  modified1:
/
/  PURPOSE: The framework for lightweight timing functions for Enzo 
/   Routines. Sets up the enzo_timing namespace, which contains the 
/   following classes:
/       section_performance: Contains timing information for a single 
/           code section
/       enzo_timer: Contains general information and section_performance 
/           objects.
/
************************************************************************/

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

  // Section Performance Class
  // Contains timers, stop/start functions, and accessors. 
  // The main things kept track of are current and total times.
  // Current is meant to be reset to 0.0 every time the timer 
  // is written out, whereas total is kept until the simulation 
  // ends.
  class section_performance
  {
  public:
    // Constructor 
    section_performance(char* myname){
      name = myname;
      next = NULL;
      ngrids = 0;
      ncells = 0;
      total_time = 0.0;
      current_time = 0.0;
    }
   
    // Start Timer 
    void start(void){
      t0 = ReturnWallTime();
    }
  
    // Stop Timer, add to current/total times
    void stop(void){
      t1 = ReturnWallTime();
      current_time += t1-t0;
      total_time += t1-t0;
    }

    // Access total_time
    double get_total_time(void){
      return total_time;
    }
    
    // Access current_time
    double get_current_time(void){
      return current_time;
    }

    // Reset the current_time timer.
    void reset_current_time(void){
      current_time = 0.0;
    }

    // Access the ncells counter
    long int get_cells(void){
      return ncells;
    }

    // Access the ngrids counter
    long int get_grids(void){
      return ngrids;
    }

    // Add the number of grids, and cells to this timer.
    void add_stats(long int my_grids, long int my_cells){
      ngrids = my_grids;
      ncells = my_cells;
    }

    std::string name;           // Name of the timer
    section_performance *next;  // Pointer to the next timer 
  
  private:
    double t0;            // Start Time
    double t1;            // End Time
    double total_time;    // Total time during the simulation
    double current_time;  // Time spent in this timer since last write-out
    long int ngrids;      // Number of Grids (For Level Timers)
    long int ncells;      // Number of Cells (For Level Timers)
    
  };


  /* --------------------------------------------------------- */

  // enzo_timer class definition
  // Contains linked list of timers, helper functions to 
  // start/stop timers, write out the data, and manage the 
  // timers. 
  class enzo_timer
  {
  public:

    // Constructor, uses default filename.
    enzo_timer(void){
      total_time = 0.0;
      current_time = 0.0;
      filename = (char *)("chronos.out");
      set_mpi_environment();
    }

    // Constructor, accepts non-standard filename
    enzo_timer(char *performance_name){
      total_time = 0.0;
      current_time = 0.0;
      filename = performance_name;
      set_mpi_environment();
    }

    // Destructor, erases each timer.
    ~enzo_timer(void){
      for( iter=timers.begin(); iter!=timers.end(); ++iter){
        delete [] iter->second;
        timers.erase(iter);
      }
      timers.clear();
    }

    // Sets up nprocs/my_rank properly.
    void set_mpi_environment(void){
#ifdef USE_MPI
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
      nprocs = 1;
      my_rank = 0;
#endif
    }

    // Use std::map objects and iterators to store timers
    // in a linked list.
    std::map<std::string, section_performance *> timers;
    std::map<std::string, section_performance *>::iterator iter;

    // Accessor for section performance object
    section_performance * get(char *name){
      iter = timers.find(name);
      if (iter == timers.end()){
        fprintf(stderr, "%s Being Created\n", name);
        timers[name] = new section_performance(name);
      }
      return timers[name];
    }

    // Start a timer by name
    void start(char *name){
      iter = timers.find(name);
      if (iter == timers.end()){
        fprintf(stderr, "%s Being Created\n", name);
        timers[name] = new section_performance(name);
      }
      timers[name]->start();
    }

    // Stop a timer by name
    void stop(char *name){
      iter = timers.find(name);
      if (iter == timers.end()){
        fprintf(stderr, "%s Being Created\n", name);
        timers[name] = new section_performance(name);
      }
      timers[name]->stop();
    }

    // Get a level section_performance by level
    section_performance * get_level(int level){
      char level_name[256];
      sprintf(level_name, "Level_%d", level);
      return get(level_name);
    }

    // Access total_time by summing from timers
    double get_total_time(void){
      total_time = 0.0;
      for( iter=timers.begin(); iter!=timers.end(); ++iter){
        total_time += iter->second->get_total_time();
      }
      return total_time;
    }

    // Get the total time for the level timers.
    double get_total_levels_time(void){
      total_time = 0.0;
      std::string keyname;
      for( iter=timers.begin(); iter!=timers.end(); ++iter){
        keyname = iter->first;
        if (strncmp(keyname.c_str(), "Level", 5) == 0){
          total_time += iter->second->get_total_time();
        }
      }
      return total_time;
    }

    // Not implemented.
    int get_total_grids(void){
      return 0;
      /* temp_performance = level_zero_performance; */
      /* int total_grids = 0; */
      /* while(temp_performance != NULL){ */
      /*     total_grids += temp_performance->get_ngrids(); */
      /*     temp_performance = temp_performance->next; */
      /* } */
      
      /* return total_grids; */
    }

    // Get sum of current_time's
    double get_current_time(void){
      current_time = 0.0;
      for( iter=timers.begin(); iter!=timers.end(); ++iter){
        current_time += iter->second->get_current_time();
      }
      return current_time;
    }

    // Get sum of level current_time's
    double get_current_levels_time(void){
      total_time = 0.0;
      std::string keyname;
      for( iter=timers.begin(); iter!=timers.end(); ++iter){
        keyname = iter->first;
        if (strncmp(keyname.c_str(), "Level", 5) == 0){
          total_time += iter->second->get_current_time();
        }
      }
      return total_time;
    }

    // Write out performance measures to a file, optionally specifying
    // verbose to get all timers from all processors.
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

      // Collect times from all processors.
      double thistime = this->get_current_levels_time();
      Reduce_Times(thistime, time_array);

      // Write out Cycle Number and calculate stats on total time
      // for this cycle.
      double mean_time = 0.0;
      double min_time, max_time, stddev_time;
      if (my_rank == 0){
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
      
        fprintf(chronos_file, "Cycle_Number %d\n",step);
      }
      
      long int total_cells = 0;
      
      // Print out info for each timer.
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
          if (strncmp(keyname.c_str(), "Total", 5) == 0){
            total_time = mean_time;
          }
          if (strncmp(keyname.c_str(), "Level", 5) == 0){
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

      // Write out total cells divided by processor-seconds.
      if (my_rank == 0){
        fprintf(chronos_file, "Total_Cells/Sec/Processor %e\n", (double)(total_cells/total_time/nprocs));
        fprintf(chronos_file, "\n");
        fclose(chronos_file);      
        delete [] time_array;
      }
    }
          
  private:
    FILE * chronos_file;    // File to write stats to
    double total_time;      // Total time
    double current_time;    // Current Time since last write_out
    char * filename;        // Filename
    int my_rank;            // MPI Rank
    int nprocs;             // MPI Size
  };
}

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

EXTERN enzo_timing::enzo_timer *enzo_timer;   // Add global timer                                                         

// Macros to start, stop, and write timers.
#define TIMING_ON
#ifdef TIMING_ON
#define TIMER_START(section_name) enzo_timer->start(section_name)
#define TIMER_STOP(section_name) enzo_timer->stop(section_name)
#define TIMER_WRITE(cycle_number) enzo_timer->write_out(cycle_number)
#else
#define TIMER_START(section_name)
#define TIMER_STOP(section_name)
#define TIMER_WRITE(cycle_number)
#endif

#endif //ENZO_TIMING
