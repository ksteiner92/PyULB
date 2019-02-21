/**
 * This file is part of edge.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * logger.h
 *
 *  Created on: 11.06.2017
 *      Author: Klaus Steiner
 */

#ifndef LOGGER_H_
#define LOGGER_H_

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <thread>
#include <iomanip>
#include <stdexcept>

#ifdef _OPENMP

#include <omp.h>

#else
///Macros used to disguise the fact that we do not have multithreading enabled.
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#endif

/* FOREGROUND */
#define RST  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST

#define BOLD(x)   "\x1B[1m" x RST
#define UNDL(x)   "\x1B[4m" x RST

enum LogLevel
{
   DEBUG, INFO, WARNING, ERROR, MSG
};

enum LogFlags
{
   ENDL, TIME
};

class Logger
{
private:
   static Logger *singleton;

   Logger();

public:
   static Logger *get();

   void log(const std::string &message);

   Logger &operator<<(const std::string &message);

   Logger &operator<<(long message);

   Logger &operator<<(const char *message);

   Logger &operator<<(LogLevel level);

   Logger &operator<<(LogFlags flag);

   static inline std::string logLevelColorize(LogLevel level, const char *msg);

   static inline std::string logLevelColorize(LogLevel level, const std::string &msg);

private:
   std::vector<std::stringstream> buffer;
   bool preinfos;
   LogLevel level;
   bool logtime;
};

/*
 * Progressbar class with multi-thread support taken from
 * https://stackoverflow.com/questions/28050669/can-i-report-progress-for-openmp-tasks
 */
///@brief Used to time how intervals in code.
///
///Such as how long it takes a given function to run, or how long I/O has taken.
class Timer
{
private:
   typedef std::chrono::high_resolution_clock clock;
   typedef std::chrono::duration<double, std::ratio<1> > second;

   std::chrono::time_point<clock> start_time; ///< Last time the timer was started
   double accumulated_time;                   ///< Accumulated running time since creation
   bool running;                              ///< True when the timer is running

public:
   Timer()
   {
      accumulated_time = 0;
      running = false;
   }

   ///Start the timer. Throws an exception if timer was already running.
   void start()
   {
      if (running)
         throw std::runtime_error("Timer was already started!");
      running = true;
      start_time = clock::now();
   }

   ///Stop the timer. Throws an exception if timer was already stopped.
   ///Calling this adds to the timer's accumulated time.
   ///@return The accumulated time in seconds.
   double stop()
   {
      if (!running)
         throw std::runtime_error("Timer was already stopped!");

      accumulated_time += lap();
      running = false;

      return accumulated_time;
   }

   ///Returns the timer's accumulated time. Throws an exception if the timer is
   ///running.
   double accumulated()
   {
      if (running)
         throw std::runtime_error("Timer is still running!");
      return accumulated_time;
   }

   ///Returns the time between when the timer was started and the current
   ///moment. Throws an exception if the timer is not running.
   double lap()
   {
      if (!running)
         throw std::runtime_error("Timer was not started!");
      return std::chrono::duration_cast<second>(clock::now() - start_time).count();
   }

   ///Stops the timer and resets its accumulated time. No exceptions are thrown
   ///ever.
   void reset()
   {
      accumulated_time = 0;
      running = false;
   }
};


///@brief Manages a console-based progress bar to keep the user entertained.
///
///Defining the global `NOPROGRESS` will
///disable all progress operations, potentially speeding up a program. The look
///of the progress bar is shown in ProgressBar.hpp.
class ProgressBar
{
private:
   uint32_t total_work;    ///< Total work to be accomplished
   uint32_t next_update;   ///< Next point to update the visible progress bar
   uint32_t call_diff;     ///< Interval between updates in work units
   uint32_t work_done;
   uint16_t old_percent;   ///< Old percentage value (aka: should we update the progress bar) TODO: Maybe that we do not need this
   Timer timer;         ///< Used for generating ETA

   ///Clear current line on console so a new progress bar can be written
   void clearConsoleLine() const
   {
      std::cerr << "\r\033[2K" << std::flush;
   }

public:
   ///@brief Start/reset the progress bar.
   ///@param total_work  The amount of work to be completed, usually specified in cells.
   void start(uint32_t total_work)
   {
      timer = Timer();
      timer.start();
      this->total_work = total_work;
      next_update = 0;
      call_diff = total_work / 200;
      old_percent = 0;
      work_done = 0;
      clearConsoleLine();
   }

   ///@brief Update the visible progress bar, but only if enough work has been done.
   ///
   ///Define the global `NOPROGRESS` flag to prevent this from having an
   ///effect. Doing so may speed up the program's execution.
   void update(uint32_t work_done0)
   {
      //Provide simple way of optimizing out progress updates
#ifdef NOPROGRESS
      return;
#endif

      //Quick return if this isn't the main thread
      if (omp_get_thread_num() != 0)
         return;

      //Update the amount of work done
      work_done = work_done0;

      //Quick return if insufficient progress has occurred
      if (work_done < next_update)
         return;

      //Update the next time at which we'll do the expensive update stuff
      next_update += call_diff;

      //Use a uint16_t because using a uint8_t will cause the result to print as a
      //character instead of a number
      uint16_t percent = (uint8_t) (work_done * omp_get_num_threads() * 100 / total_work);

      //Handle overflows
      if (percent > 100)
         percent = 100;

      //In the case that there has been no update (which should never be the case,
      //actually), skip the expensive screen print
      if (percent == old_percent)
         return;

      //Update old_percent accordingly
      old_percent = percent;

      //Print an update string which looks like this:
      //  [================================================  ] (96% - 1.0s - 4 threads)
      std::cerr << "\r\033[2K["
                << std::string(percent / 2, '|') << std::string(50 - percent / 2, ' ')
                << "] ("
                << percent << "% - "
                << std::fixed << std::setprecision(1) << timer.lap() / percent * (100 - percent)
                << "s - "
                << omp_get_num_threads() << " threads)" << std::flush;
   }

   ///Increment by one the work done and update the progress bar
   ProgressBar &operator++()
   {
      //Quick return if this isn't the main thread
      if (omp_get_thread_num() != 0)
         return *this;

      work_done++;
      update(work_done);
      return *this;
   }

   ///Stop the progress bar. Throws an exception if it wasn't started.
   ///@return The number of seconds the progress bar was running.
   double stop()
   {
      clearConsoleLine();
      timer.stop();
      return timer.accumulated();
   }

   ///@return Return the time the progress bar ran for.
   double time_it_took()
   {
      return timer.accumulated();
   }

   uint32_t cellsProcessed() const
   {
      return work_done;
   }
};

#define _LOG (*Logger::get())

#define LOG(Level) _LOG << LogLevel::Level

#define LOG_T(Level) _LOG << LogLevel::Level << LogFlags::TIME

#endif /* LOGGER_H_ */
