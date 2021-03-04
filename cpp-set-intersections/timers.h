
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <unistd.h>
#include <errno.h>

#include <pthread.h>

#include <sys/time.h>
#include <time.h>



#ifndef timersub
void
timersub(struct timeval *a, struct timeval *b, struct timeval *res)
{
    res->tv_usec = a->tv_usec - b->tv_usec;
    res->tv_sec = a->tv_sec - b->tv_sec;
    if (res->tv_usec > a->tv_usec) {
        res->tv_sec--;
        res->tv_usec += 1000;
    }
}
#endif

#define CLOCKTIME



#ifndef CLOCKTIME

class Timer {

public:

    Timer() { if(1) start(); }

    struct timeval s;

        void start() {
                gettimeofday(&s,0);
        }


    double duration()
    {
        struct timeval end, dt;

        gettimeofday(&end,0);

#ifdef timersub
        timersub(&end,&s, &dt);
#endif


        return dt.tv_sec * 1000.0 + dt.tv_usec/1000.0;
    }





};

#endif 


class Timer {

public:

  timespec st,end;

  Timer() { start(); }

  void start() { 
      clock_gettime( CLOCK_MONOTONIC, &st ); 
  }

  double duration() {
          clock_gettime( CLOCK_MONOTONIC, &end );
          double diff;
          #define BILLION 1000000000L

          //printf("%-15s: %10jd.%03ld (", "CLOCK1", (intmax_t) end.tv_sec, end.tv_nsec / 1000000);

          //printf("Sec = %llu %llu", (int64_t)st.tv_sec,(int64_t) end.tv_sec);

          diff = BILLION * (end.tv_sec - st.tv_sec)  + (end.tv_nsec - st.tv_nsec);
          printf("elapsed process CPU time = %lf nanoseconds\n", (double) diff);
          return (double)(diff/1e6);
  }


};


//#define Timer ClockTime


