/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/
#pragma once
#ifndef __SPDCHOL_H__
#define __SPDCHOL_H__

#include <cmath>
#include <limits>
#include <numeric>
#include "spdbase.h"

namespace mns 
{
	template <typename T>
	class SpdChol : public SpdBase<T> 
	{
	// Calculates Cholesky decomposition and solves the system of linear equations with symmetric positive-definite matrix
	// Cholesky factor is being updated by means of Givens rotations
	public:
		SpdChol(SpdMatrixT&& spdMatrixT, int n) : SpdBase(n, false, T(0.0)), m_(std::move(spdMatrixT)) {};
		const SpdMatrixT& GetMatrix();
		~SpdChol() {};
	private:
		// Interface implementation
		virtual Status FactorizeImpl() override;
		virtual Status SolveImpl(VectorT& b) const override;
		virtual T	   GetRCondImpl() const override;
		virtual Status UpdateAddImpl(VectorT& a) override final;
		virtual Status UpdateDelImpl(int ix) override final;
		void   Compress(int ix);
		void   GetGivensRotation(T x, T y, T& c, T& s) const;

		SpdMatrixT m_;
	};

	template<typename T> 
	const typename SpdChol<T>::SpdMatrixT& SpdChol<T>::GetMatrix()
	{
		int n = GetMatrixDim(); 
		m_.resize(((Size_T)n) * (n + 1) / 2); 
		return m_; 
	};

	template<typename T> 
	Status SpdChol<T>::FactorizeImpl()
	{
	// Computes Cholessky factor of the symmetric positive-definite matrix
		if( IsFactorized() )
		{
			return Status::Success;
		}

		int n = GetMatrixDim();

		for( int i = 0; i < n; ++i ) 
		{
			for( int k = 0; k <= i; ++k ) 
			{
				T s = T(0.0);
				for( int j = 0; j < k; ++j )
				{
					s += m_[j + ((Size_T)i) * (i + 1) / 2] * m_[j + ((Size_T)k) * (k + 1) / 2]; // m[i][j] * m[k][j]
				}

				if ( i == k )
				{
					Size_T ii = i + ((Size_T)i) * (i + 1) / 2;
					T d = m_[ii]  - s;
					if( d <= std::numeric_limits<T>::epsilon() )
					{
					// Matrix is not positive definite one
						this->isFactorized_ = false;
						return Status::IllConditionedMatrix;
					}
					m_[ii] = std::sqrt(d);
				}
				else
				{
					Size_T ik = k + ((Size_T)i) * (i + 1) / 2;
					m_[ik] = ( T(1.0) / m_[k + ((Size_T)k) * (k + 1) / 2] * (m_[ik] - s) ); // m[k][k] * (m[i][k] - s))
				}
			}
		}
		this->isFactorized_ = true;
		return Status::Success;
	}

	template<typename T> 
	Status SpdChol<T>::SolveImpl(VectorT& b) const
	{
	// Solves the linear equations system using the Cholesky decomposition
		if( !IsFactorized() )
		{
			return Status::Failure;
		}

		int n = GetMatrixDim();

		if( b.size() < n )
		{
			return Status::BadParameter;
		}

		T  s;
		for( int i = 0; i < n; ++i )  
		{
			s = b[i];
			for( int j = 0; j < i; ++j ) 
			{
				s -= m_[j + ((Size_T)i) * (i + 1) / 2] * b[j];
			}
			b[i] = s / m_[i + ((Size_T)i) * (i + 1) / 2];
		}

		for( int i = n - 1; i >= 0; --i ) 
		{
    		b[i] /= m_[i + ((Size_T)i) * (i + 1) / 2]; 
			for( int j = 0; j < i; ++j )
			{
				b[j] -= m_[j + ((Size_T)i) * (i + 1) / 2] * b[i];
			}
		}

		return Status::Success;
	}

	template <typename T>
	T SpdChol<T>::GetRCondImpl() const
	{ 
	//  Alg. 5.1 from
	//  Nicholas J. Higham "A Survey of Condition Number Estimation for Triangular Matrices" // SIAM Review Vol.29, No.4, 1987
	//  http://eprints.ma.man.ac.uk/695/01/covered/MIMS_ep2007_10.pdf

		T cond = T();
		if( !this->isFactorized_ )
		{
			return cond;
		}
		else
		{
			int n = GetMatrixDim();
			T norm1 = T(0.0);
			for( int j = 0; j < n; ++j )  
			{
				T s = T();
				for( int i = 0; i < j; ++i ) 
				{
					s += fabs(m_[i + ((Size_T)j) * (j + 1) / 2]);
				}
				for( int i = j; i < n; ++i ) 
				{
					s += fabs(m_[j + ((Size_T)i) * (i + 1) / 2]);
				}
				if ( s > norm1 )
				{
					norm1 = s;
				}
			}

			T renorm1;
			VectorT x(n, T(1.0)/n);
			VectorT e(n);
			int ix;
			for( int k = 0; k < 4; ++k )
			{
				SolveImpl(x);

				for( int i = 0; i < n; ++i )
				{
					e[i] = ( x[i] >= T(0.0) ) ? T(1.0) : T(0.0);
				}

				SolveImpl(e);

				T r;
				if( k == 0 )
				{
					r = std::accumulate(std::begin(e), std::end(e), 0.0) / n; 
				}
				else
				{
					r = e[ix]; 
				}

				T maxAbsEl = std::fabs(e[0]);
				ix = 0;
				for( int i = 0; i < n; ++i ) 
				{
					T w = std::fabs(e[i]);
					if( w > maxAbsEl )
					{
						ix = i;
						maxAbsEl = w;  
					}
				}

				if ( maxAbsEl <= r )
				{
					renorm1 = T(0.0);
					for( int i = 0; i < n; ++i ) 
					{
						renorm1 += fabs(x[i]);
					}
					break;
				}

				x.assign(n, T(0.0));
				x[ix] = T(1.0);
			}
			cond = norm1 * renorm1;
		}
		return T(1.0) / cond;
	}

	template<typename T> 
	Status SpdChol<T>::UpdateAddImpl(VectorT& d)
	{
		// Updates the Cholesky factor after a symmetric column/row	addition
		// d -  new matrix column

		if( !IsFactorized() )
		{
			return Status::Failure;
		}

		int n = GetMatrixDim();
		Size_T msize = n * (n + 1) / 2;
		if( d.size() < n + 1 )
		{
			return Status::BadParameter;
		}

		if( m_.size() < msize + n + 1 )
		{
			m_.resize(msize + n + 1);
		}

		// Calculate a new row of the matrix decomposition
		// Solve L * y = d 
		T s;
		int i, j, k;
		for( j = 0; j < n; ++j ) 
		{
			s = d[j];
			for( k = 0; k < j; ++k ) 
			{
				s -= m_[k + ((Size_T)j) * (j + 1) / 2] * d[k];
			}
			d[j] = s / m_[j + ((Size_T)j) * (j + 1) / 2];
		}

		s = T(0.0);
		for( i = 0; i < n; ++i ) 
		{
			s += d[i] * d[i];
		}

		s = d[n] - s;
		if( s <= std::numeric_limits<T>::epsilon() ) 
		{
			return Status::IllConditionedMatrix;
		}

		d[n] = std::sqrt(s);
		
		for( i = 0; i <= n; ++i ) 
		{
			m_[msize + i] = d[i];
		}

		++this->n_;
		return Status::Success;
	}

	template<typename T> 
	Status SpdChol<T>::UpdateDelImpl(int ix)
	{
		int n = GetMatrixDim();
		// Calculates a new Cholesky factor for a matrix with deleted row and column 
		if ( ix < 0 || ix > n - 1 )
		{
			return Status::BadParameter;
		}

		if( ix < n - 1 )
		{
			Size_T ii1, ii2;
			T m1, m2, c, s;
			for ( int i = ix; i < n - 1; ++i )
			{
				int ip1 = i + 1;
				ii1 = i   + ((Size_T)ip1) * (ip1 + 1) / 2; 
				ii2 = ip1 + ((Size_T)ip1) * (ip1 + 1) / 2; 
				m1 = m_[ii1];
				m2 = m_[ii2];
				GetGivensRotation(m1, m2, c, s);
				m_[ii1] =  c * m1 + s * m2;
				m_[ii2] = -s * m1 + c * m2;
				if ( i < n_ - 2 )
				{
					for ( int k = i + 2; k < n_; ++k )
					{
						ii1 = i   + ((Size_T)k) * (k + 1) / 2; 
						ii2 = ip1 + ((Size_T)k) * (k + 1) / 2;
						m1 = m_[ii1];
						m2 = m_[ii2];
						m_[ii1] =  c * m1 + s * m2;
						m_[ii2] = -s * m1 + c * m2;
					}
				}
			}

			Compress(ix);
		}           // if ( ix == n - 1 ) then we do nothing besides
		--this->n_; // reducing the matrix dimension
		return Status::Success;
	}

	template<typename T> 
	void  SpdChol<T>::Compress(int ix)
	{
		if( ix < n_ - 1 )
		{
			Size_T ij = ((Size_T)ix) * (ix + 1) / 2;
			for ( int i = ix + 1; i < n_; ++i )
			{
				for ( int j = 0; j < i; ++j )
				{
					m_[ij++] = m_[j + ((Size_T)i) * (i + 1) / 2];
				}
			}
		} 
	}

	template<typename T> 
	void SpdChol<T>::GetGivensRotation(T x, T y, T& c, T& s) const
	{
	// Computes a Givens plane rotation for values x and y

		T ax, ay, t, u, w, r;

		ax = std::fabs(x);
		ay = std::fabs(y);
		t = std::max(ax, ay);
		u = std::min(ax, ay);

		if( t != T(0.0) ) 
		{
			w = u / t;
			r = t * std::sqrt(T(1.0) + w * w);
			c = x / r;
			s = y / r;
		}
		else 
		{
			c = T(1.0);
			s = T(0.0);
		};
	}


} // end of MNS namespace

#endif // __SPDCHOL_H__
