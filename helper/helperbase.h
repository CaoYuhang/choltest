/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/
#pragma once
#ifndef __HELPERBASE_H__
#define __HELPERBASE_H__

#include "../common/defs.h"

namespace mns 
{
	inline Size_T GetIndex(int i, int j) { return (i >= j) ? j + ((Size_T)i) * (i + 1) / 2 : i + ((Size_T)j) * (j + 1) / 2; };

	template <typename T>
	class HelperBase
	{
	// Defines interface for common vector/matrix operations
	public:
		typedef typename Defs<T>::VectorT VectorT;
		typedef typename Defs<T>::SpdMatrixT SpdMatrixT;

		VectorT GetResidual(int n, const SpdMatrixT& a, const VectorT& x, const VectorT& b) const { return GetResidualImpl(n, a, x, b); };
		T		GetVectorNorm2(Size_T n, const VectorT& v) const { return GetVectorNorm2Impl(n, v); }
		virtual ~HelperBase() {};
	protected:
		virtual VectorT GetResidualImpl(int n, const SpdMatrixT& a, const VectorT& x, const VectorT& b) const = 0; 
		virtual T GetVectorNorm2Impl(Size_T n, const VectorT& v) const = 0;

		HelperBase() {};
	private:
    	HelperBase(const HelperBase&);
		HelperBase& operator =(const HelperBase&);
		HelperBase& operator =(HelperBase&&);
	};

} // end of MNS namespace

#endif // __HELPERBASE_H__
