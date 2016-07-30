
#include "mixture_mvn.hpp"

#include <cmath>
#include <cstdlib>
#include <numeric>
#include <cassert>
#include <stdexcept>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


//#define MYDEBUG
#include <iostream>

using namespace std;


namespace subpavings { 
	
	namespace kde {
		
		MixtureMVN::MixtureMVN(const vector < vector < double > >& m,
					const vector < vector < vector < double > > >& s,
					const vector < double >& mx,
					unsigned long int rsd) 
					: mixn(m.size()), mixes(mx), r(NULL)
		{
			if (!mixn)
				throw std::invalid_argument("MixtureMVN(...) : m empty");
			if (mixn != s.size())
				throw std::invalid_argument("MixtureMVN(...) : m and s different sizes");
			if (mixn != mx.size())
				throw std::invalid_argument("MixtureMVN(...) : m and mx different sizes");
			if (m.front().empty())	
				throw std::invalid_argument("MixtureMVN(...) : m first element empty");
				
			dim = m.front().size();
			
			try {
				
				double mixTotal = std::accumulate(mixes.begin(), mixes.end(), 0.0);
				
				for (size_t i = 0; i < mixn; ++i) {
					
					vector < double > thisMean = m[i];
					vector < vector < double > > thisScale = s[i];
					
					
					
					#ifdef MYDEBUG
						std::cout << "\nindex is i " << std::endl;
						std::cout << "mean is " << std::endl;
						for (size_t k = 0; k < thisMean.size(); ++k) cout << thisMean[k] << "\t";
						std::cout << std::endl;
						
						std::cout << "scale  is " << std::endl;
						for (size_t k = 0; k < thisScale.size(); ++k) {
							for (size_t q = 0; q < thisScale[k].size(); ++q) cout << thisScale[k][q] << "\t";
							std::cout << std::endl;
						}
						std::cout << std::endl;
					
					#endif
					
					
					
					if (thisMean.size() != dim) 
						throw std::invalid_argument("MixtureMVN(...) : m inconsistent sizes");
					if (thisScale.size() != dim) 
						throw std::invalid_argument("MixtureMVN(...) : s inconsistent sizes");
					means.push_back(gsl_vector_alloc (dim));
					scales.push_back(gsl_matrix_alloc (dim, dim));
					
					for (size_t j = 0; j < dim; ++j) {
						gsl_vector_set(means[i], j, thisMean[j]);
						
						if (thisScale[j].size() != dim) 
						throw std::invalid_argument("MixtureMVN(...) : s element not square");
					
						for (size_t k = 0; k < dim; ++k) {
							gsl_matrix_set(scales[i], j, k, thisScale[j][k]);
						}
					}
					
					// makes sure mixes all add to 1.0 or close ...
					mixes[i] /= mixTotal;
					
				}// end i loop
				
				resetPRNG(rsd);
				
			}
			catch  (...) {
				clean();
			}
		}
		
		
		MixtureMVN::~MixtureMVN()
		{
			clean();
			
		}
		
		void MixtureMVN::resetPRNG(unsigned long int rsd)
		{
			if (r != NULL) gsl_rng_free(r);
			r = gsl_rng_alloc (gsl_rng_mt19937); 
			
			seed = rsd;
		
			gsl_rng_set (r, seed);
			
		}
		
		void MixtureMVN::resetPRNG()
		{
			if (r != NULL) gsl_rng_free(r);
			r = gsl_rng_alloc (gsl_rng_mt19937); 
			
			gsl_rng_set (r, seed); 
			
		}
		
		std::vector < std::vector < double > >&
			MixtureMVN::prn(std::vector < std::vector < double > >& rvs,
							size_t n) const
		{
			std::vector < std::vector < double > >tmp;
			tmp.reserve(n);
			
			for (size_t i = 0; i < n; ++i) {
				tmp.push_back(vector < double >());
				prn(tmp.back());
			}
			
			rvs.swap(tmp);
			return rvs;
			
		}
		
		std::vector < std::vector < real > >&
			MixtureMVN::prn(std::vector < std::vector < real > >& rvs,
							size_t n) const
		{
			std::vector < std::vector < real > >tmp;
			tmp.reserve(n);
			
			for (size_t i = 0; i < n; ++i) {
				tmp.push_back(vector < real >());
				prn(tmp.back());
			}
			
			rvs.swap(tmp);
			return rvs;
			
		}
		
		std::vector < cxsc::rvector >&
			MixtureMVN::prn(std::vector < cxsc::rvector >& rvs,
							size_t n) const
		{
			std::vector < cxsc::rvector >tmp;
			tmp.reserve(n);
			
			for (size_t i = 0; i < n; ++i) {
				tmp.push_back(cxsc::rvector());
				prn(tmp.back());
			}
			
			rvs.swap(tmp);
			return rvs;
			
		}
		
		std::vector < double >& MixtureMVN::prn(std::vector < double >& rv) const
		{
			size_t index = getDistIndex();
			
			gsl_vector* x = gsl_vector_alloc (dim);
			
			prn_gsl_vec(index, 	x);
			std::vector < double >result(dim);
			for (size_t i = 0; i < dim; ++i) result[i] = gsl_vector_get(x, i);
			rv.swap(result);
			
			return rv;
		}
		
		std::vector < real >& MixtureMVN::prn(std::vector < real >& rv) const
		{
			size_t index = getDistIndex();
			
			gsl_vector* x = gsl_vector_alloc (dim);
			
			prn_gsl_vec(index, 	x);
			std::vector < real >result(dim);
			for (size_t i = 0; i < dim; ++i) result[i] = gsl_vector_get(x, i);
			rv.swap(result);
			
			return rv;
		}
		
		cxsc::rvector& MixtureMVN::prn(cxsc::rvector& rv) const
		{
			size_t index = getDistIndex();
			
			gsl_vector* x = gsl_vector_alloc (dim);
			
			prn_gsl_vec(index, 	x);
			cxsc::Resize(rv, dim);
			int lb = Lb(rv);
			for (size_t i = 0; i < dim; ++i) rv[i+lb] = gsl_vector_get(x, i);
			
			return rv;
		}
		
		real MixtureMVN::f(const std::vector < double >& dv) const
		{
			
			gsl_vector* x = gsl_vector_alloc (dim);
			for (size_t i = 0; i < dim; ++i) gsl_vector_set(x, i, dv[i]);
			
			real result(0.0);
			for (size_t i = 0; i < mixn; ++i) {
				result += (mixes[i] * d_gsl_vec(i, x));
			}
				
			return result;
		}
		
		real MixtureMVN::f(const std::vector < real >& rv) const
		{
			
			gsl_vector* x = gsl_vector_alloc (dim);
			for (size_t i = 0; i < dim; ++i) gsl_vector_set(x, i, _double(rv[i]));
			
			real result(0.0);
			for (size_t i = 0; i < mixn; ++i) {
				result += (mixes[i] * d_gsl_vec(i, x));
			}
				
			return result;
		}
		
		real MixtureMVN::f(const cxsc::rvector& rv) const
		{
			
			gsl_vector* x = gsl_vector_alloc (dim);
			int lb = Lb(rv);
			for (size_t i = 0; i < dim; ++i) gsl_vector_set(x, i, _double(rv[i+lb]));
			
			real result(0.0);
			for (size_t i = 0; i < mixn; ++i) {
				result += (mixes[i] * d_gsl_vec(i, x));
			}
				
			return result;
		}
		
		size_t MixtureMVN::getDistIndex() const 
		{
			/* generate a random number u in [0, 1) and use to pick the distn*/
			double u = gsl_rng_uniform(r); 
			
			size_t index = 0;
			double cmix = 0.0;
			for (; index < mixn-1 ; ++index) {
				cmix += mixes[index];
				if (u < cmix) break;
			}
			
			assert (((u < cmix) && index < (mixn-1)) 
					|| (!(u < cmix) && index == (mixn-1)));
			
			return index;
		}
		
		
		void MixtureMVN::clean()
		{
			try {
				for (size_t i = 0; i < means.size(); ++i) {
					if (NULL != means[i]) gsl_vector_free(means[i]);
					means[i] = NULL;
				}
				
			}
			catch(...) {} // catch and swallow
			try {
				for (size_t i = 0; i < scales.size(); ++i) {
					
					if (NULL != scales[i]) gsl_matrix_free(scales[i]);
					scales[i] = NULL;
				}
				
			}
			catch(...) {} // catch and swallow
			try {
				
				if (NULL != r) gsl_rng_free(r);
				r = NULL;
			}
			catch(...) {} // catch and swallow
			
		}
		
		
		
		void MixtureMVN::prn_gsl_vec(size_t i, 	gsl_vector *result) const
		{
		
			gsl_matrix *work = NULL;
			try {
				work = gsl_matrix_alloc(dim,dim);

				gsl_matrix_memcpy(work,scales[i]);
				
				gsl_linalg_cholesky_decomp(work);

				for(int k = 0; k < dim; ++k)
					gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

				gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
				gsl_vector_add(result,means[i]);

				gsl_matrix_free(work);
				work = NULL;
			}
			catch(...) {
				try {
					if (work != NULL) gsl_matrix_free(work);
				}
				catch(...) {}
				try {
					int n = result->size;
					gsl_vector_free(result);
					result = gsl_vector_calloc(n);
				}
				catch(...) {}
				throw;
			}

		}

		double MixtureMVN::d_gsl_vec(size_t i, const gsl_vector *x) const
		{
		
			gsl_vector *ym = NULL;
			gsl_vector *xm = NULL;
			gsl_matrix *work = NULL;
			gsl_matrix *winv = NULL;
			gsl_permutation *p = NULL;

			try {
				double ax,ay;
				int s;
				work = gsl_matrix_alloc(dim,dim), 
				winv = gsl_matrix_alloc(dim,dim);
				p = gsl_permutation_alloc(dim);

				gsl_matrix_memcpy( work, scales[i] );
				gsl_linalg_LU_decomp( work, p, &s );
				gsl_linalg_LU_invert( work, p, winv );
				ax = gsl_linalg_LU_det( work, s );
				
				gsl_matrix_free( work );
				work = NULL;
				gsl_permutation_free( p );
				p = NULL;

				xm = gsl_vector_alloc(dim);
				gsl_vector_memcpy( xm, x);
				gsl_vector_sub( xm, means[i] );
				ym = gsl_vector_alloc(dim);
				gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
				
				gsl_matrix_free( winv );
				winv = NULL;
				
				gsl_blas_ddot( xm, ym, &ay);
				
				gsl_vector_free(xm);
				xm = NULL;
				gsl_vector_free(ym);
				ym = NULL;
				
				ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),dim)*ax );

				return ay;
			}
			catch (...) {
				try { if (xm != NULL) gsl_vector_free(xm); }
				catch (...) {}
				try { if (ym != NULL) gsl_vector_free(ym); }
				catch (...) {}
				try { if (work != NULL) gsl_matrix_free(work); }
				catch (...) {}
				try { if (winv != NULL) gsl_matrix_free(winv); }
				catch (...) {}
				try { if (p != NULL) gsl_permutation_free(p); }
				catch (...) {}
				
				throw;
			}
		}

	
	}
}
