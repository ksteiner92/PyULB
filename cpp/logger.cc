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
 * logger.cc
 *
 *  Created on: 11.06.2017
 *      Author: Klaus Steiner
 */

#ifdef _OPENMP

#include <omp.h>

#endif

#include <ctime>
#include "logger.h"

using namespace std;

Logger::Logger() : level(INFO), logtime(false), preinfos(false)
{
   size_t num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
   {
#pragma omp single
      num_threads = omp_get_num_threads();
   }
#endif
   buffer.resize(num_threads);
}

Logger *Logger::get()
{
   if (singleton == nullptr)
      singleton = new Logger();
   return singleton;
}

inline string Logger::logLevelColorize(LogLevel level, const char *msg)
{
   stringstream ss;
   switch (level) {
      case DEBUG:
         ss << KBLU << msg << RST;
         break;
      case INFO:
         ss << KGRN << msg << RST;
         break;
      case WARNING:
         ss << KYEL << msg << RST;
         break;
      case ERROR:
         ss << KRED << msg << RST;
         break;
   }
   return ss.str();
}

inline string Logger::logLevelColorize(LogLevel level, const string &msg)
{
   return logLevelColorize(level, msg.c_str());
}

void Logger::log(const string &message)
{
   size_t thread = 0;
#ifdef _OPENMP
   thread = omp_get_thread_num();
#endif
   if (!preinfos && level != MSG) {
      switch (level) {
         case DEBUG:
            buffer[thread] << FBLU("[Debug");
            break;
         case INFO:
            buffer[thread] << FGRN("[Info");
            break;
         case WARNING:
            buffer[thread] << FYEL("[Warning");
            break;
         case ERROR:
            buffer[thread] << FRED("[Error");
            break;
      }
      if (logtime) {
         time_t t = time(nullptr);
         tm tm = *localtime(&t);
         stringstream ss;
         time_t rawtime;
         struct tm *timeinfo;
         char tbuff[80];
         time(&rawtime);
         timeinfo = localtime(&rawtime);
         strftime(tbuff, 80, "%T", timeinfo);
         ss << ", " << tbuff;
         buffer[thread] << logLevelColorize(level, ss.str());
      }
      buffer[thread] << logLevelColorize(level, "]: ");
      preinfos = true;
   }
   buffer[thread] << message;
}

Logger &Logger::operator<<(LogLevel level)
{
   this->level = level;
   return *this;
}

Logger &Logger::operator<<(const string &message)
{
   log(message);
   return *this;
}

Logger &Logger::operator<<(long message)
{
   log(to_string(message));
   return *this;
}

Logger &Logger::operator<<(LogFlags flag)
{
   size_t thread = 0;
#ifdef _OPENMP
   thread = omp_get_thread_num();
#endif
   switch (flag) {
      case ENDL:
         cerr << buffer[thread].str() << endl;
         preinfos = false;
         logtime = false;
         buffer[thread].str("");
         break;
      case TIME:
         logtime = true;
         break;
   }
   return *this;
}

Logger &Logger::operator<<(const char *message)
{
   log(string(message));
   return *this;
}

Logger *Logger::singleton = nullptr;
