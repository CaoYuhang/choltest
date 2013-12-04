/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/
#pragma once
#ifndef __SPDBASE_H__
#define __SPDBASE_H__

#include "../common/defs.h"

namespace mns 
{
	template <typename T>
	class SpdBase
	{
	// Defines interface for solving the system of linear equations with symmetric positive-definite matrix
	public:
		typedef typename Defs<T>::VectorT VectorT;
		typedef typename Defs<T>::SpdMatrixT SpdMatrixT;

		Status Factorize() { return FactorizeImpl(); };
		Status UpdateAdd(VectorT& a) { return UpdateAddImpl(a); };
		Status UpdateDel(int ix) { return UpdateDelImpl(ix); };
		Status Solve(VectorT& b) const { return SolveImpl(b); };
		T      GetRCond() const { return GetRCondImpl(); };

		int    GetMatrixDim() const { return GetMatrixDimImpl(); };
		bool   IsFactorized() const { return IsFactorizedImp(); } //{ return GetRCondImpl() > T(0.0); };
		explicit operator bool() const { return IsFactorized(); }
		
		virtual ~SpdBase() = default;

		SpdBase(const SpdBase&) = delete;
		SpdBase& operator =(const SpdBase&) = delete;
		SpdBase& operator =(SpdBase&&) = delete;
	private:
		virtual	Status FactorizeImpl() = 0;
		virtual Status SolveImpl(VectorT& b) const  = 0;
		virtual Status UpdateAddImpl(VectorT& a) = 0;
		virtual Status UpdateDelImpl(int ix) = 0;
		virtual int    GetMatrixDimImpl() const = 0;
		virtual T	   GetRCondImpl() const  = 0;
		virtual bool   IsFactorizedImp() const  = 0;
	};

} // end of MNS namespace

#endif // __SPDBASE_H__
