/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/
#pragma once
#ifndef __DEFS_H__
#define __DEFS_H__

#include <array>
#include <vector>

namespace mns 
{
	template <typename T, int N=0>
	class Defs
	{
	public:
		typedef typename std::vector<T> VectorT;
		typedef typename VectorT SpdMatrixT;
		typedef typename std::array<T, N> ArrayT;
	};

	typedef std::size_t Size_T;

/* MNS Status Codes */
	enum class Status : unsigned int 
	{
		Success = 0x0,
		BadParameter = 0x1,
		IllConditionedMatrix = 0x3,
		IterationLimit = 0x7,
		OutOfMempory = 0xB,
		Failure = 0xE
	};

} // end of MNS namespace

#endif /* __DEFS_H__ */

