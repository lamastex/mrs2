/*   
     Computes the integrals ...
     

     Author:  Tomas Johnson <johnson@math.uu.se>
     Latest edit: jan 27 13:43:40 CET 2010

     By: Tomas Johnson <johnson@math.uu.se>
*/


#include "Int.h"

class subDomain{
public:
  ivector sDomain;
  real sTol;
};

//Performance variables:
int split(0);

interval integrate(taylor::dim2taylor (*integrand)(taylor::dim2taylor_vector), 
		   ivector &domain, int order) {
 
  order=2*(order/2);
  ivector midDomain((Inf(domain)+Sup(domain))/2.0);
  taylor::dim2taylor_vector mp;
  mp=taylor::init_var(order,midDomain);  
  taylor::dim2taylor midInt=integrand(mp);
  interval result(0.0);
  interval rad1=interval(Sup(domain[1])-midDomain[1]);
  interval rad2=interval(Sup(domain[2])-midDomain[2]);
  for(int i=0; i<=order; i+=2){
    for(int j=0; j<=order-i; j+=2){
      //      result+=interval(4.0)*interval(binomial[i+j][j])*midInt[i][j]*power(rad1,i+1)*power(rad2,j+1)
      //	/((i+1.0)*(j+1.0));
      result+=interval(4.0)*midInt[i][j]*power(rad1,i+1)*power(rad2,j+1)
      	/((i+1.0)*(j+1.0));
    }
  }
  taylor::dim2taylor_vector rem=taylor::init_var(order+2,domain);
  taylor::dim2taylor remInt=integrand(rem);
  interval remTerm(0.0);
  for(int i=0; i<=order+2; i++) {
    remTerm+=remInt[i][order+2-i]
      *domain[1]*domain[2]*power(-rad1|rad1,i)*power(-rad2|rad2,order+2-i);
  }

  return result+remTerm;

}


interval integrateWithSplitting(taylor::dim2taylor (*integrand)(taylor::dim2taylor_vector), 
				const ivector &domain, int order, real tol) {
	     
  subDomain workDomain;
  workDomain.sDomain=domain;
  workDomain.sTol=tol;
  stack<subDomain> workStack;
  workStack.push(workDomain);
  interval result(0.0);
  while(!workStack.empty()) {
    workDomain=workStack.top();
    workStack.pop();
  
      interval sResult=integrate(integrand,workDomain.sDomain,order);
      if(diam(sResult)<workDomain.sTol)
	//if(diam(sResult)/Sup(sResult)<workDomain.sTol)
	result+=sResult;
      else {
	subDomain newDomain1,newDomain2;
	if(diam(workDomain.sDomain[1])<=diam(workDomain.sDomain[2])) {
	  newDomain1.sDomain=workDomain.sDomain;
	  newDomain1.sDomain[2]=interval(Inf(workDomain.sDomain[2]),mid(workDomain.sDomain[2]));
	  newDomain2.sDomain=workDomain.sDomain;
	  newDomain2.sDomain[2]=interval(mid(workDomain.sDomain[2]),Sup(workDomain.sDomain[2]));
	}
	else {
	  newDomain1.sDomain=workDomain.sDomain;
	  newDomain1.sDomain[1]=interval(Inf(workDomain.sDomain[1]),mid(workDomain.sDomain[1]));
	  newDomain2.sDomain=workDomain.sDomain;
	  newDomain2.sDomain[1]=interval(mid(workDomain.sDomain[1]),Sup(workDomain.sDomain[1]));
	}
	newDomain1.sTol=workDomain.sTol/2.0;
	newDomain2.sTol=workDomain.sTol/2.0; //2
	workStack.push(newDomain1);
	workStack.push(newDomain2);
	++split;
	//	cout<<workStack.size()<<endl;
      }
    
  }
  cout<<"The Integral is enclosed by : "<<result<<endl
      <<"The Domain was split "<<split<<" times"<<endl;

  return result;    

}

//gloria's addition for integration for absolute error calculations
interval integrate(taylor::dim2taylor (*integrand)(taylor::dim2taylor_vector, interval), 
		   interval fhat, ivector &domain, int order) {
 
  order=2*(order/2);
  ivector midDomain((Inf(domain)+Sup(domain))/2.0);
  taylor::dim2taylor_vector mp;
  mp=taylor::init_var(order,midDomain);  
  taylor::dim2taylor midInt=integrand(mp, fhat);
  interval result(0.0);
  interval rad1=interval(Sup(domain[1])-midDomain[1]);
  interval rad2=interval(Sup(domain[2])-midDomain[2]);
  for(int i=0; i<=order; i+=2){
    for(int j=0; j<=order-i; j+=2){
      //      result+=interval(4.0)*interval(binomial[i+j][j])*midInt[i][j]*power(rad1,i+1)*power(rad2,j+1)
      //	/((i+1.0)*(j+1.0));
      result+=interval(4.0)*midInt[i][j]*power(rad1,i+1)*power(rad2,j+1)
      	/((i+1.0)*(j+1.0));
    }
  }
  taylor::dim2taylor_vector rem=taylor::init_var(order+2,domain);
  taylor::dim2taylor remInt=integrand(rem, fhat);
  interval remTerm(0.0);
  for(int i=0; i<=order+2; i++) {
    remTerm+=remInt[i][order+2-i]
      *domain[1]*domain[2]*power(-rad1|rad1,i)*power(-rad2|rad2,order+2-i);
  }

  return result+remTerm;

}

interval integrateWithSplitting(taylor::dim2taylor (*integrand)(taylor::dim2taylor_vector, interval), 
				interval fhat, const ivector &domain, int order, real tol) {
	     
  subDomain workDomain;
  workDomain.sDomain=domain;
  workDomain.sTol=tol;
  stack<subDomain> workStack;
  workStack.push(workDomain);
  interval result(0.0);
  while(!workStack.empty()) {
    workDomain=workStack.top();
    workStack.pop();
  
      interval sResult=integrate(integrand,fhat,workDomain.sDomain,order);
      if(diam(sResult)<workDomain.sTol)
	//if(diam(sResult)/Sup(sResult)<workDomain.sTol)
	result+=sResult;
      else {
	subDomain newDomain1,newDomain2;
	if(diam(workDomain.sDomain[1])<=diam(workDomain.sDomain[2])) {
	  newDomain1.sDomain=workDomain.sDomain;
	  newDomain1.sDomain[2]=interval(Inf(workDomain.sDomain[2]),mid(workDomain.sDomain[2]));
	  newDomain2.sDomain=workDomain.sDomain;
	  newDomain2.sDomain[2]=interval(mid(workDomain.sDomain[2]),Sup(workDomain.sDomain[2]));
	}
	else {
	  newDomain1.sDomain=workDomain.sDomain;
	  newDomain1.sDomain[1]=interval(Inf(workDomain.sDomain[1]),mid(workDomain.sDomain[1]));
	  newDomain2.sDomain=workDomain.sDomain;
	  newDomain2.sDomain[1]=interval(mid(workDomain.sDomain[1]),Sup(workDomain.sDomain[1]));
	}
	newDomain1.sTol=workDomain.sTol/2.0;
	newDomain2.sTol=workDomain.sTol/2.0; //2
	workStack.push(newDomain1);
	workStack.push(newDomain2);
	++split;
	//	cout<<workStack.size()<<endl;
      }
    
  }
  cout<<"The Integral is enclosed by : "<<result<<endl
      <<"The Domain was split "<<split<<" times"<<endl;

  return result;    

}
////////////////////////////////////////////////////////////////////////
