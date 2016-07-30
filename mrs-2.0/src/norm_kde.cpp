/*
* Copyright (C) 2012 Jennifer Harlow
*
* This file is part of mrs, a C++ class library for statistical set processing.
*
* mrs is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/*! \file
\brief Definitions for a class for making a kernel density 
* estimate using normal reference rule 
* (Scott, 'Multivariate Density Estimation: Theory, Practice, and Visualization',
* 1992, pge 152).

*/

#include "norm_kde.hpp"


#include <cmath>
#include <cstdlib> // srand

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <cassert>

/*18/11/2003: multivariate density estimation: product kernel*/

/*16/11/2003*/
/*negative log posterior: to sample the bandwidths*/

//#define DEBUG
		
namespace subpavings { 
	
	namespace kde {
	
		
		NormKDE::NormKDE(const std::vector < std::vector < real > >& dx)
			: data_x(dx),
			productBandwidths(0.0)
		
		{
			data_num = data_x.size();
			if (!data_num) throw std::invalid_argument("NormKDE(...): data_x empty");
			dim = data_x.front().size();
			if (!dim) throw std::invalid_argument("NormKDE(...): data dim = 0");
			cont = cxsc::exp(-0.5*dim*cxsc::ln(cxsc::Pi2_real));
		
		}	
		
		NormKDE::NormKDE(const std::vector < std::vector < double > >& dx)
			: productBandwidths(0.0)
		
		{
			data_num = dx.size();
			
			if (!data_num) throw std::invalid_argument("NormKDE(...): data_x empty");
			dim = dx.front().size();
			if (!dim) throw std::invalid_argument("NormKDE(...): data dim = 0");
			cont = cxsc::exp(-0.5*dim*cxsc::ln(cxsc::Pi2_real));
		
			data_x.reserve(data_num);
			for (size_t i = 0; i < data_num; ++i) {
				data_x.push_back(std::vector < real >(dx[i].begin(), dx[i].end()));
			}
		
		}			


		cxsc::real NormKDE::kde(const std::vector < cxsc::real >& x) const
		{
			
			if (bandwidths.empty())
				throw std::runtime_error("NormKDE::kde(...) : no bandwidths set");
			
			// rescale
			std::vector < cxsc::real > rx(x);
			std::transform( rx.begin(), rx.end(), bandwidths.begin(), 
			rx.begin(), std::divides <real>() );
			
			return _kde(rx);
			
		}
		
		cxsc::real NormKDE::kde(const std::vector < double >& x) const
		{
			// rescale
			std::vector < cxsc::real > rx(x.begin(), x.end());
			std::transform( rx.begin(), rx.end(), bandwidths.begin(), 
			rx.begin(), std::divides <real>() );
			
			return _kde(rx);
			
		}
		
		cxsc::real NormKDE::kde(const cxsc::rvector& x) const
		{
			if (bandwidths.empty())
				throw std::runtime_error("NormKDE::kde(...) : no bandwidths set");
			
			int dd = VecLen(x);
			// rescale
			std::vector < cxsc::real > rx(dd);
			for (int d = 1; d < dd; ++d) rx[d-1] = x[d]/bandwidths[d-1];
			
			return _kde(rx);
			
		}
			
		std::vector < cxsc::real > NormKDE::df(	
									const std::string& logFilename) const 
		{
			std::vector < cxsc::real >().swap(bandwidths);
			std::vector < real > h(dim, 0.0);
			
			std::ofstream ofs;
			ofs.open(logFilename.c_str());         // overwrite
			
			if (ofs.is_open()) {
				
				ofs << "sample size = " << data_num << std::endl;
				ofs << std::endl; 
				ofs.close();
			}
			else {
				std::cerr << "Error: could not open file named "
				<< logFilename << std::endl << std::endl;
			} 
			
			/*normal reference rule, or rule of thumb*/
			for(size_t j = 0; j < dim; ++j) 
			{ 
				cxsc::real temp1 = 0.0; 
				cxsc::real temp2 = 0.0; 
				for(size_t i = 0; i < data_num; ++i) 
				{ 
					temp1 += data_x[i][j]; 
					temp2 += data_x[i][j]*data_x[i][j]; 
					
				}
				cxsc::real rdnum(1.0*data_num); 
				cxsc::real sigma=cxsc::sqrt(temp2/rdnum-(temp1/rdnum)*(temp1/rdnum)); 
				cxsc::real temp = cxsc::exp(1.0/(dim+4.0)*std::log(4.0/(dim+2.0))); 
				h[j] = temp*sigma*cxsc::exp(-1.0/(dim+4)*log(1.0*data_num)); 
			} 
			
			ofs.open(logFilename.c_str(), std::ios::app);         // append
			
			if (ofs.is_open()) {
			
				ofs << "normal reference rule:" << std::endl;
				ofs << std::fixed << std::setprecision(8); 
				for(size_t j = 0; j < dim; ++j) ofs << _double(h[j]) << std::endl; 
				ofs << std::endl; 
				ofs.close();
			}
			else {
				std::cerr << "Error: could not open file named "
				<< logFilename << std::endl << std::endl;
			} 
			
			productBandwidths = h[0];
			for(size_t i = 1; i < dim; ++i) productBandwidths *= h[i];
			

			bandwidths.swap(h);
			
			return bandwidths; 

		}  

		//assumes rescaled data
		cxsc::real NormKDE::_kde(const std::vector < cxsc::real >& rx) const
		{
			if (bandwidths.empty())
				throw std::runtime_error("NormKDE::kde(...) : Bandwidths not set");
			
			if (rescaled_x.empty()) rescaleData();
						
			/*for each (rescaled) sample pt, calculate f( ( x_d - x_i,d) /h_d ) 
			 * and acculate produce over d*/
			
						
			cxsc::real sum(0.0);  
			for(size_t i = 0; i < data_num; ++i) { 
				cxsc::real temp = 0.0;
				for(size_t j = 0; j < dim; ++j) {  
					cxsc::real xa = rx[j] - rescaled_x[i][j]; 
					temp += xa*xa; 
				} 
				sum += cont*cxsc::exp(-0.5*temp); 
			}
			cxsc::real hatf = sum/(productBandwidths * (1.0*data_num)); 
			
			return hatf;
			
		}
		
		
		void NormKDE::rescaleData() const
		{
			rescaled_x = data_x;
			/*rescale by dividing by h_j, j = 1... dim */
			for (size_t i = 0; i < data_num; ++i) {
				
				std::transform( rescaled_x[i].begin(), rescaled_x[i].end(),
					bandwidths.begin(),
					rescaled_x[i].begin(), std::divides<real>() );
			}
			
		}
	}
}
