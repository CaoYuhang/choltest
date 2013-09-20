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

		int    GetMatrixDim() const { return n_; };
		bool   IsFactorized() const { return isFactorized_; };

		// explicit operator bool() const { return isFactorized_; }

		virtual ~SpdBase() {};
	protected:
		virtual	Status FactorizeImpl() = 0;
		virtual Status SolveImpl(VectorT& b) const = 0;
		virtual Status UpdateAddImpl(VectorT& a) { return Status::Failure; };
		virtual Status UpdateDelImpl(int ix) { return Status::Failure; };
		virtual T	   GetRCondImpl() const { return T(); };

		SpdBase() {};
		SpdBase(int n, bool isFactorized, T cond) : n_(n), isFactorized_(isFactorized), cond_(cond) {};

		int n_;
		bool isFactorized_;
		T cond_;
	private:
    	SpdBase(const SpdBase&);
		SpdBase& operator =(const SpdBase&);
		SpdBase& operator =(SpdBase&&);
	};

} // end of MNS namespace

#endif // __SPDBASE_H__
