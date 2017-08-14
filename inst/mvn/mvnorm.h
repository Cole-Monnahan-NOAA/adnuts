//  --------------------------------------------------------------------------
// Copyright (c) 2008,2009,2010,2011,2012, Anders Nielsen <an@aqua.dtu.dk> 
// and Casper Berg <cbe@aqua.dtu.dk>. All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//   * Neither the name of the assessment tool SAM nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// ARE DISCLAIMED. IN NO EVENT SHALL ANDERS NIELSEN OR CASPER BERG BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY 
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
// DAMAGE.
//  --------------------------------------------------------------------------

#ifndef __mvnorm_h__
#define __mvnorm_h__

#include <df1b2fun.h>
#include <math.h>

const double log2pi = log(2.0*M_PI);

dvar_vector bksolve(const dmatrix& L, const dvar_vector& b)
{
  RETURN_ARRAYS_INCREMENT();
  int i, k, r, R, c, C;
  r=b.indexmin();
  R=b.indexmax();
  dvariable sum;
  dvar_vector x(r,R);
  for(i=r; i<=R; ++i){
    sum=b(i);
    for(k=i-1; k>=r; k--){
      sum-=L(i,k)*x(k);
    }
    x(i)=sum/L(i,i);
  }
  RETURN_ARRAYS_DECREMENT();
  return x;
}

class MVNORM_t{
  int r, R, N;
  dmatrix L;      
  double logdet; 
  dmatrix Sigma;  
public:
  MVNORM_t(){}
  MVNORM_t(dmatrix Sigma_){
    setSigma(Sigma_);
  } 
  dmatrix cov(){return Sigma;}
  void setSigma(dmatrix Sigma_){
    r=Sigma_.rowmin();
    R=Sigma_.rowmax();
    N=R-r+1;
    Sigma = Sigma_;
    double logdet=0;
    L=choleski_decomp(Sigma);
    for(int i=r; i<=R; ++i){logdet+=log(L(i,i));}
    logdet*=2.0;
  }
  dvariable operator()(dvar_vector x){
    dvar_vector tmp=bksolve(L,x);
    dvariable ret=0.5*(log2pi*N+logdet+sum(square(tmp)));
    return ret;
  }
};
#endif
