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
* estimate using a data-adaptive method with a MCMC to get
* the bandwidth estimates.  See
* (Zhang, X., King, Maxwell, Hyndman, Rob. J. (2006), 'A Bayesian approach
* to bandwidth selection for multivariate kernel density estimation',
* Computational Statistics and Data Analysis, vol. 50, pp. 3009--3031).

*/

#include "mcmc_kde.hpp"


#include <cmath>
#include <cstdlib> // srand

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <ctime>
#include <cassert>

/*18/11/2003: multivariate density estimation: product kernel*/

/*16/11/2003*/
/*negative log posterior: to sample the bandwidths*/

//#define DEBUG
		
namespace subpavings { 
	
	namespace kde {
	
		const cxsc::real MCMCKDE::mutsizp = 0.15;

		MCMCKDE::MCMCKDE(const std::vector < std::vector < real > >& dx)
			: accept_h(0), total_h(1), iset(0), gset(0.0), data_x(dx),
			productBandwidths(0.0)
		
		{
			data_num = data_x.size();
			if (!data_num) throw std::invalid_argument("MCMCKDE(...): data_x empty");
			dim = data_x.front().size();
			if (!dim) throw std::invalid_argument("MCMCKDE(...): data dim = 0");
			cont = cxsc::exp(-0.5*dim*cxsc::ln(cxsc::Pi2_real));
		
		}	
		
		MCMCKDE::MCMCKDE(const std::vector < std::vector < double > >& dx)
			: accept_h(0), total_h(1), iset(0), gset(0.0),
			productBandwidths(0.0)
		
		{
			data_num = dx.size();
			
			if (!data_num) throw std::invalid_argument("MCMCKDE(...): data_x empty");
			dim = dx.front().size();
			if (!dim) throw std::invalid_argument("MCMCKDE(...): data dim = 0");
			cont = cxsc::exp(-0.5*dim*cxsc::ln(cxsc::Pi2_real));
		
			data_x.reserve(data_num);
			for (size_t i = 0; i < data_num; ++i) {
				data_x.push_back(std::vector < real >(dx[i].begin(), dx[i].end()));
			}
		
		}			


		cxsc::real MCMCKDE::kde(const std::vector < cxsc::real >& x) const
		{
			
			// rescale
			std::vector < cxsc::real > rx(x);
			std::transform( rx.begin(), rx.end(), bandwidths.begin(), 
			rx.begin(), std::divides <real>() );
			
			return _kde(rx);
			
		}
		
		cxsc::real MCMCKDE::kde(const std::vector < double >& x) const
		{
			// rescale
			std::vector < cxsc::real > rx(x.begin(), x.end());
			std::transform( rx.begin(), rx.end(), bandwidths.begin(), 
			rx.begin(), std::divides <real>() );
			
			return _kde(rx);
			
		}
		
		cxsc::real MCMCKDE::kde(const cxsc::rvector& x) const
		{
			int dd = VecLen(x);
			// rescale
			std::vector < cxsc::real > rx(dd);
			for (int d = 1; d < dd; ++d) rx[d-1] = x[d]/bandwidths[d-1];
			
			return _kde(rx);
			
		}
			

		std::vector < cxsc::real > MCMCKDE::df(	size_t size_batch, 
									size_t num_batch, 
									size_t warm, 
									size_t step,
									const std::string& logFilename,
									unsigned int seed) const 
		{
			std::string resultsFilename("");
			return df(	size_batch, 
						num_batch, 
						warm, 
						step,
						logFilename,
						resultsFilename,
						seed);
			
		}

		std::vector < cxsc::real > MCMCKDE::df(	size_t size_batch, 
									size_t num_batch, 
									size_t warm, 
									size_t step,
									const std::string& logFilename,
									const std::string& resultsFilename,
									unsigned int seed) const 
		{
			std::vector < cxsc::real >().swap(bandwidths);
			size_t M = size_batch*num_batch;
			
			std::vector < cxsc::real > x(dim, 0.0);  // initialise to 0.0
			std::vector < real > h(dim, 0.0);
			
			// for logging results
			std::vector < std::vector < real > >logh;
			std::vector < real > costs;
			bool logging = !resultsFilename.empty();
			if (logging) {
				logh.reserve(ceil((M + warm)/(1.0*step))+1);
				costs.reserve(ceil((M + warm)/(1.0*step))+1);
			}
			
			//std::vector < std::vector < cxsc::real > > 
						//cov(dim, std::vector < cxsc::real >(dim));
			
			std::vector < std::vector < cxsc::real > > batch_h;
			batch_h.reserve(num_batch);
			
			std::vector < cxsc::real > sum_h(dim, 0.0);
			std::vector < cxsc::real > var_h(dim, 0.0);
			
			// seed rand
			srand(seed); 

			std::ofstream ofs;
			ofs.open(logFilename.c_str());         // overwrite
			
			if (ofs.is_open()) {
				
				ofs << "random seed = " << seed << std::endl; 
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

			

			/*Initial values*/
			cxsc::real xCost = cost(x);

			std::cout << "Initial cost = " << _double(xCost) << std::endl;
			
			if(logging) {
				costs.push_back(xCost);
				/* initial xs are all 0 so initial h all 1.0*/
				logh.push_back(std::vector < cxsc::real >(dim, 1.0));
			}
			
			clock_t startWarm = clock();
			
			for(size_t k = 0; k < warm; ++k) {
				xCost = kn_gibbs(xCost, x);
				
				if(logging && (k%step == 0)) {
					costs.push_back(xCost);
					logh.push_back(std::vector < cxsc::real >(dim, 0.0));
					for(size_t j = 0; j < dim; ++j) logh.back()[j] = cxsc::exp(x[j]); 
				}
			}
			
			clock_t endWarm = clock();
			double timeWarm = (static_cast<double>(endWarm-startWarm)/CLOCKS_PER_SEC);	
			std::cout << "Warm finished, cost = " << _double(xCost);
			std::cout << ", time = " << timeWarm << " sec" << std::endl;
			
			//reset totals
			total_h = 1; // we give it one state to start with
			accept_h = 0; 

			clock_t startEst = clock();
			
			for(size_t k = 0; k < M; ++k) { 
				xCost = kn_gibbs(xCost, x);
				
				bool logThisOne = (logging && (k%step == 0));
				if(logThisOne) {
					logh.push_back(std::vector < cxsc::real >(dim, 0.0));
					costs.push_back(xCost);
				}
				
				if (k%size_batch == 0) 
					batch_h.push_back(std::vector < cxsc::real >(dim, 0.0));				
				
				for(size_t j = 0; j < dim; ++j) {
					real this_h = cxsc::exp(x[j]);
					sum_h[j] += this_h;
					batch_h.back()[j] += this_h; 
					if(logThisOne)
						logh.back()[j] = this_h; 
				} 

			} 
			
			clock_t endEst = clock();
			double timeEst = (static_cast<double>(endEst-startEst)/CLOCKS_PER_SEC);	
			std::cout << "Bandwidth estimation finished: ";
			std::cout << " estimation time = " << timeEst << " sec" << std::endl;
			std::cout << " (total time = " << (timeEst+timeWarm) << " sec)" << std::endl;
						
			// log before assert
			if (logging) logResults(resultsFilename, costs, logh, warm, step);
		
			assert(batch_h.size() == num_batch);
			

			for(size_t i = 0; i < dim; ++i) { 
				sum_h[i]=sum_h[i]/(1.0*M); 
				for(size_t j = 0; j < num_batch; ++j) {
					batch_h[j][i] /= (1.0*size_batch); 
					cxsc::real temp = batch_h[j][i] - sum_h[i]; 
					var_h[i] += temp*temp; 
				} 
				var_h[i] = cxsc::sqrt(var_h[i]/(1.0*num_batch*num_batch-num_batch)); 
			} 

			ofs.open(logFilename.c_str(), std::ios::app);         // append
			
			if (ofs.is_open()) {
			
				ofs << "Bayesian estimates of bandwidths:" << std::endl;
				ofs << "accept rate = " << ((1.0*accept_h)/total_h) << std::endl;
				ofs << std::fixed << std::setprecision(8); 
				for(size_t i = 0; i < dim; ++i) ofs << _double(sum_h[i]) 
						<< "\t" << _double(var_h[i]) << std::endl; 
				ofs << "\n\nbatch means:" <<std::endl;
				
				/* output batch means as well */
				for(size_t j = 0; j < num_batch; ++j) { 
					ofs << (j+1);
					for(size_t i = 0; i < dim; ++i) ofs << "\t" << _double(batch_h[j][i]);
					ofs << std::endl;
				} 
				 
				ofs.close();
			}
			else {
				std::cerr << "Error: could not open file named "
				<< logFilename << std::endl << std::endl;
			} 
			
			productBandwidths = sum_h[0];
			for(size_t i = 1; i < dim; ++i) productBandwidths *= sum_h[i];
			
			bandwidths.swap(sum_h);
			
			return bandwidths; 
		}  

		//assumes rescaled data
		cxsc::real MCMCKDE::_kde(const std::vector < cxsc::real >& rx) const
		{
			if (bandwidths.empty())
				throw std::runtime_error("MCMCKDE::kde(...) : Bandwidths not set");
			
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
		

		cxsc::real MCMCKDE::cost(const std::vector < cxsc::real >& h) const
		{
			std::vector < cxsc::real > he(dim);

			cxsc::real hprod = 1.0;
			for(size_t k = 0; k < dim; ++k) 
			{ 
				he[k]=cxsc::exp(h[k]); 
				hprod *= he[k]; 
			}
			 
			/*	
			if(he[1]<=0.001) return 1.0*exp(20.0); 
			if(he[2]<=0.001) return 1.0*exp(20.0); 
			*/ 
			
			
			#ifdef DEBUG
				std::cout << "hprod = " << _double(hprod) 
					<< "\tcont = " << _double(cont) << std::endl;
			#endif 
			
			cxsc::real sum=0.0;  
			for(size_t j=1; j < data_num; ++j) 
			{ 
				cxsc::real temp = 0.0;
				for(size_t k = 0; k < dim; ++k)  
				{  
					cxsc::real xa=(data_x[0][k] - data_x[j][k])/he[k]; 
					temp += xa*xa; 
				} 
				//sum += cont*cxsc::exp(-0.5*temp)/hprod; 
				sum += cxsc::exp(-0.5*temp); 
			} 
			
			cxsc::real hatf = (cont*sum/hprod)/(data_num - 1.0); 
			cxsc::real cv = cxsc::ln(hatf); 
			
			#ifdef DEBUG
				std::cout << "sum = " << _double(sum) << "\thatf = " << _double(hatf) 
						<< "\tcv = " << _double(cv) << std::endl;
			#endif
			
			for(size_t i= 1; i < data_num-1; ++i) 
			{ 
				cxsc::real thisSum = 0.0; 
				for(size_t j = 0; j < i; ++j) 
				{ 
					cxsc::real temp = 0.0;
					for(size_t k = 0; k < dim; ++k) 
					{ 
						cxsc::real xa=(data_x[i][k] - data_x[j][k])/he[k]; 
						temp += xa*xa; 
					} 
					thisSum += cxsc::exp(-0.5*temp); 
					//thisSum += cont*cxsc::exp(-0.5*temp)/hprod; 
				} 
				for(size_t j = i+1; j < data_num; ++j) 
				{ 
					cxsc::real temp = 0.0;
					for(size_t k = 0; k < dim; ++k) 
					{ 
						cxsc::real xa=(data_x[i][k]-data_x[j][k])/he[k]; 
						temp += xa*xa; 
					} 
					thisSum += cxsc::exp(-0.5*temp); 
					//thisSum += cont*cxsc::exp(-0.5*temp)/hprod; 
				} 
				//cxsc::real thisHatf=thisSum/(data_num - 1.0); 
				cxsc::real thisHatf=(cont*thisSum/hprod)/(data_num - 1.0); 

				cv += cxsc::ln(thisHatf); 
			}
			
			#ifdef DEBUG
				std::cout << "cv = " << cv << std::endl;
			#endif
			
			cxsc::real finalSum=0.0; 
			for(size_t j = 0; j < data_num-1; ++j) 
			{ 
				cxsc::real temp = 0.0;
				for(size_t k = 0; k < dim; ++k) 
				{ 
					cxsc::real xa=(data_x[data_num-1][k]-data_x[j][k])/he[k]; 
					temp += xa*xa; 
				} 
				finalSum += cxsc::exp(-0.5*temp); 
				//finalSum += cont*cxsc::exp(-0.5*temp)/hprod; 
			} 
			//cxsc::real finalHatf = finalSum/(data_num - 1.0); 
			cxsc::real finalHatf = (cont*finalSum/hprod)/(data_num - 1.0); 
			cv += cxsc::ln(finalHatf);
				
			for(size_t i = 0;i < dim; ++i) {
				/*log Jacobi*/ 
				cv += h[i]; 
				/*log priors*/
				cv += -1.0*cxsc::ln(1.0+he[i]*he[i]); 
			}
			#ifdef DEBUG
				std::cout << "finalSum = " << _double(finalSum) 
					<< "\tfinalHatf = " << _double(finalHatf) 
					<< "\tcv = " << _double(cv) << std::endl;
			#endif
			return -1.0*cv; 
		} 


		cxsc::real MCMCKDE::kn_gibbs(
				cxsc::real xCost,
				std::vector < cxsc::real >& x) const
		{
			
			std::vector < cxsc::real > dv(dim);
			std::vector < cxsc::real > rn(dim);
			
			cxsc::real sum = 0.0;
			for(size_t i = 0; i < dim; ++i) { 
				rn[i] = gasdev(); 
				sum += rn[i]*rn[i]; 
			} 
			
			std::vector < cxsc::real > temp_x(x);
			
			
			for(size_t i = 0; i < dim; ++i) { 
				dv[i] = rn[i]/cxsc::sqrt(sum)*gasdev()*mutsizp; 
				temp_x[i] += dv[i]; 
			}
			
			 
			cxsc::real yCost = cost(temp_x);
			
			cxsc::real r = xCost - yCost; 
			
			bool accept = false;
			if(r > 0.0) accept = true; 
			else 
			{ 
				cxsc::real un = 0.0; 
				while( !(un > 0.0)) un = rand()*1.0/RAND_MAX; 
				if(un < exp(r)) accept = true; 
				
			} 
			if(accept) 
			{ 
				++accept_h;  
				xCost = yCost;
				x = temp_x;
			} 
			// else x unchanged
			 
			++total_h;

			return xCost; 
		} 




		real MCMCKDE::gasdev() const
		{ 
			if(!iset) 
			{ 
				real r(0.0);
				real v1(0.0); 
				real v2(0.0); 
					
				do 
				{ 
					v1 = rand()*2.0/RAND_MAX - 1.0; 
					v2 = rand()*2.0/RAND_MAX - 1.0; 
					r = v1*v1 + v2*v2; 
				} 
				while(!(r < 1.0)); 
				
				real fac = cxsc::sqrt(-2.0*cxsc::ln(r)/r); 
				gset = v1*fac; 
				iset = 1; 
				
				return v2*fac; 
			} 
			else 
			{ 
				iset = 0; 
				
				return gset; 
			} 
		} 
		
		void MCMCKDE::rescaleData() const
		{
			rescaled_x = data_x;
			/*rescale by dividing by h_j, j = 1... dim */
			for (size_t i = 0; i < data_num; ++i) {
				
				std::transform( rescaled_x[i].begin(), rescaled_x[i].end(),
					bandwidths.begin(),
					rescaled_x[i].begin(), std::divides<real>() );
			}
			
		}
		
		
		void MCMCKDE::logResults(const std::string& resultsFilename,
				const std::vector < cxsc::real >& costs,
				const std::vector <std::vector < cxsc::real > >& results,
				size_t warm, size_t step) const
		{
			std::ofstream rofs;
			rofs.open(resultsFilename.c_str());         // overwrite
			
			if (rofs.is_open()) {
				
				size_t n = resultsFilename.size();
				if (costs.size() < n) n = costs.size();
						
				rofs << std::fixed; 
				size_t ind = 0;
				for(size_t i = 0; i < n; ++i) { 
					if (i == warm/step+1) ind = warm+1;
					else if (i) ind += step;
					
					rofs << ind << "\t";
					rofs << std::setprecision(2) << _double(costs[i]); 
						rofs << std::setprecision(8); 
						for(size_t j = 0; j < dim; ++j)  rofs << "\t" << _double(results[i][j]); 
						rofs << std::endl;
				} 
				 
				rofs.close(); 
				 
			}
			else {
				std::cerr << "Error: could not open file named "
				<< resultsFilename << std::endl << std::endl;
			} 
		}

				
	} // end namespace subpavings::kde
} // end namespace subpavings


