/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/
#pragma once
#ifndef __HELPER_H__
#define __HELPER_H__

#include "helperbase.h"
#include <cmath>
#include <functional>
#include <numeric>

namespace mns 
{
	template <typename T>
	class Helper1 final : public HelperBase<T>
	{
	// Implements common vector/matrix operations
	public:
		Helper1() {};
	private:
		virtual VectorT GetResidualImpl(int n, const SpdMatrixT& a, const VectorT& x, const VectorT& b) const override;
		virtual T GetVectorNorm2Impl(Size_T n, const VectorT& v) const override;

		Helper1(const Helper1&);
		Helper1& operator =(const Helper1&);
		Helper1& operator =(Helper1&&);
	};

	template<typename T> 
	T Helper1<T>::GetVectorNorm2Impl(Size_T n, const VectorT& v) const
	{
		T s = T();
		s = inner_product(v.begin(), v.end(), v.begin(), s);
		return std::sqrt(s);
	}

	template<typename T> 
	typename Helper1<T>::VectorT Helper1<T>::GetResidualImpl(int n, const SpdMatrixT& a, const VectorT& x, const VectorT& b) const
	{
		VectorT r(n);
		for( int i = 0; i < n; ++i )
		{
			r[i] = T();
			for( int j = 0; j < i; ++j )
			{
				r[i] += a[j + ((Size_T)i) * (i + 1) / 2] * x[j];
			}
			for( int j = i; j < n; ++j )
			{
				r[i] += a[i + ((Size_T)j) * (j + 1) / 2] * x[j];
			}
			r[i] = b[i] - r[i];
		}
		return r;
	}


} // end of MNS namespace

#endif // __HELPER_H__
