#pragma once
#ifndef __STOPWATCH_H__
#define __STOPWATCH_H__

#include <chrono>
#include <ctime>

namespace mns 
{
typedef std::chrono::high_resolution_clock clock;
typedef std::chrono::microseconds microseconds;
typedef std::chrono::milliseconds milliseconds;

class StopWatch final
{
	public:
		StopWatch() : start_( clock::now()  ) {}

		clock::time_point Restart();
		
		double StopWatch::Elapsed() const;

		milliseconds ElapsedMs() const;

		microseconds ElapsedUs() const;

	    clock::time_point Now() const;

		static std::time_t Convert(clock::time_point time_point);
	private:
		clock::time_point start_;

    	StopWatch(const StopWatch&);
		StopWatch& operator =(const StopWatch&);
		StopWatch& operator =(StopWatch&&);
};

} // end of MNS namespace

#endif // __STOPWATCH_H__
