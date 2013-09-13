#include "stopwatch.h"

namespace mns 
{
	double StopWatch::Elapsed() const
	{ 
		return 0.001 * std::chrono::duration_cast<milliseconds>(clock::now() - start_).count(); 
	}


	microseconds StopWatch::ElapsedUs() const
	{ 
		return std::chrono::duration_cast<microseconds>(clock::now() - start_);
	}
		
	milliseconds StopWatch::ElapsedMs() const
	{
		return std::chrono::duration_cast<milliseconds>(clock::now() - start_); 
	}

	clock::time_point StopWatch::Now() const
	{
		return clock::now(); 
	}

	clock::time_point StopWatch::Restart()
	{ 
		start_ = clock::now();
		return start_; 
	}

	std::time_t StopWatch::Convert(clock::time_point time_point)
	{
		return std::chrono::system_clock::to_time_t(time_point);
	}	
} 
 

 