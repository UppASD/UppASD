// Stopwatch class
// Niklas Fejes 2012-2013

#pragma once

#include <list>
#include <string>

#include "c_headers.hpp"
#include "real_type.h"
#include "stopwatch.hpp"

class StopwatchPool {
private:
   class StopwatchNode {
   public:
      std::string name;
      Stopwatch watch;

      StopwatchNode(const std::string &n) : name(n) {
      }
   };

   std::list<StopwatchNode> watchlist;

public:
   // Keep the list even if no timing is done
   Stopwatch &get(const char *name) {
      return get(std::string(name));
   }

   Stopwatch &get(const char *name, usd_int len) {
      return get(std::string(name, len));
   }

   Stopwatch &get(const std::string &name) {
      // Check if the name is already in the list
      std::list<StopwatchNode>::iterator it;
      for(it = watchlist.begin(); it != watchlist.end(); it++) {
         if(name.compare(it->name) == 0) {
            return it->watch;
         }
      }

      // Otherwise add new
      StopwatchNode node(name);
      watchlist.push_back(node);
      return watchlist.back().watch;
   }

#ifdef DUMMY_STOPWATCH
   void printAll() {
   }

   bool empty() {
      return true;
   }
#else
   // Prints the result
   void printAll() {
      std::list<StopwatchNode>::iterator it;

      // Minimum length
      int minlen = 5;

      // Find the length of the longest name
      for(it = watchlist.begin(); it != watchlist.end(); it++) {
         int len = it->watch.minNameLen();
         if(minlen < len) {
            minlen = len;
         }
      }

      // Display all
      for(it = watchlist.begin(); it != watchlist.end(); it++) {
         std::printf("%s timing:\n", it->name.c_str());
         it->watch.print(2, minlen);
      }
   }

   bool empty() {
      std::list<StopwatchNode>::iterator it;
      for(it = watchlist.begin(); it != watchlist.end(); it++) {
         if(!it->watch.empty()) {
            return false;
         }
      }
      return true;
   }
#endif  // #ifdef STOPWATCH_DUMMY
};

class GlobalStopwatchPool {
private:
   class __StopwatchPool : public StopwatchPool {
   public:
      ~__StopwatchPool() {
         if(!empty()) {
            std::printf("\n\n");
            std::printf("========================================\n");
            std::printf("============ Stopwatch pool ============\n");
            std::printf("========================================\n");
            printAll();
            std::printf("========================================\n");
            std::printf("============ Stopwatch pool ============\n");
            std::printf("========================================\n");
         }
      }
   };

   static __StopwatchPool pool;

public:
   static Stopwatch &get(const char *name) {
      return pool.get(name);
   }

   static Stopwatch &get(const char *name, usd_int len) {
      return pool.get(name, len);
   }

   static Stopwatch &get(const std::string &name) {
      return pool.get(name);
   }

   // NOTE Will be printed automatically when the program finishes
   static void printAll() {
      pool.printAll();
   }
};

