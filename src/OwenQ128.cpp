// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <cmath>

namespace mp = boost::multiprecision;
namespace sf = boost::math;

mp::float128 pnorm_mp(mp::float128 x){
  const mp::float128 sqrt2 = 1.4142135623730950488016887242096980785696718753769Q;
  return sf::erfc(-x / sqrt2)/2;
}
mp::float128 pnorm_mp(mp::float128 x);

// [[Rcpp::export]]
double ptowen(double q, int nu, double delta){
  if(nu==1){
    double C = R::pnorm(-delta*sqrt(nu/(nu+q*q)), 0.0, 1.0, 1, 0);
    C += 2*sf::owens_t(delta*sqrt(nu/(nu+q*q)), R::sign(q)*sqrt(q*q/nu));
    return C;
  }
  const mp::float128 q_mp(q);
  const mp::float128 a = R::sign(q) * mp::sqrt(q_mp*q_mp/nu);
  const mp::float128 b = nu/(nu+q_mp*q_mp);
  const mp::float128 dsB = delta * mp::sqrt(b);
  const mp::float128 delta2(delta*delta);
  const mp::float128 twopi = 6.2831853071795864769252867665590057683943387987502Q;
  std::vector<mp::float128> M (nu-1);
  M[0] = a * mp::sqrt(b) * mp::exp(-delta2*b/2)/mp::sqrt(twopi) * pnorm_mp(a*dsB);
  if(nu>2){
    M[1] = b * (delta * a * M[0] + a * mp::exp(-delta2/2) / twopi);
    if(nu>3){
      std::vector<mp::float128> A(nu-3);
      A[0] = 1;
      int k;
      if(nu>4){
        for(k=1; k<nu-3; k++){
          A[k] = 1.0/k/A[k-1];
        }
      }
      for(k=2; k<nu-1; k++){
        M[k] = (k-1) * b * (A[k-2] * delta * a * M[k-1] + M[k-2L]) / k;
      }
    }
  }
  int i;
  if(nu%2==1){
    double C = R::pnorm(-delta*sqrt(nu/(nu+q*q)), 0.0, 1.0, 1, 0);
    C += 2*sf::owens_t(delta*sqrt(nu/(nu+q*q)), R::sign(q)*sqrt(q*q/nu));
    mp::float128 sum;
    for(i=1; i<nu-1; i+=2){
      sum += M[i];
    }
    mp::float128 res = C + 2*sum;
    return res.convert_to<double>();
  }
  mp::float128 sum;
  for(i=0; i<nu-1; i+=2){
    sum += M[i];
  }
  mp::float128 res = pnorm_mp(-delta) + mp::sqrt(twopi)*sum;
  return res.convert_to<double>();
}

// [[Rcpp::export]]
double OwenQ1(int nu, double t, double delta, double R){
  double C;
  if(nu % 2 == 1){
    double a = R::sign(t)*sqrt(t*t/nu);
    double ab = 0;
    double b = nu/(nu+t*t);
    double sB = sqrt(b);
    if(fabs(t)<DBL_MAX){
      ab = a*b;
    }
    C = R::pnorm(R, 0.0, 1.0, 1, 0) - (delta >=0) + 
      2*(sf::owens_t(delta*sB, a) - sf::owens_t(R, a-delta/R) - 
        sf::owens_t(delta*sB, (ab-R/delta)/b));
    if(nu==1){
      return C;
    }
  }
  const mp::float128 a(R::sign(t)*sqrt(t*t/nu));
  const mp::float128 b(nu/(nu+t*t));
  const mp::float128 sB = mp::sqrt(b);
  const mp::float128 R2(R*R);
  mp::float128 ab;
  mp::float128 asB;
  if(fabs(t)>DBL_MAX){
    ab = 0;
    asB = R::sign(t);
  }else{
    ab = a*b;
    asB = R::sign(t)*sqrt(t*t/(nu+t*t));
  }
  const mp::float128 twopi = 6.2831853071795864769252867665590057683943387987502Q;
  const mp::float128 sqrt2pi = mp::sqrt(twopi);
  const int n = nu-1;
  std::vector<mp::float128> H(n);
  std::vector<mp::float128> M(n);
  H[0] = - mp::exp(-R2/2)/sqrt2pi * pnorm_mp(a*R-delta);
  M[0] = asB*mp::exp(-delta*delta*b/2)/sqrt2pi * 
    (pnorm_mp(delta*asB)-pnorm_mp((delta*ab-R)/sB));
  if(nu >= 3){
    H[1] = R * H[0];
    M[1] = delta*ab*M[0] + ab*mp::exp(-delta*delta*b/2) * 
      (mp::exp(-delta*delta*a*a*b/2)-mp::exp(-(delta*ab-R)*(delta*ab-R)/b/2))/twopi;
    if(nu >= 4){
      std::vector<mp::float128> A(n);
      std::vector<mp::float128> L(n-2);
      A[0] = 1;
      A[1] = 1;
      L[0] = ab * R * mp::exp(-R2/2) * mp::exp(-(a*R-delta)*(a*R-delta)/2) / 2 / twopi;
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          L[k] = A[k+2] * R * L[k-1];
        }
      }
      for(k=2; k<n; k++){
        H[k] = A[k] * R * H[k-1];
        M[k] = (k-1.0)/k * (A[k-2] * delta * ab * M[k-1] + b*M[k-2]) - L[k-2];
      }
    }
  }
  if(nu % 2 == 0){
    mp::float128 sum;
    int i;
    for(i=0; i<nu-1; i+=2){
      sum += M[i]+H[i];
    }
    mp::float128 res = pnorm_mp(-delta) + sqrt2pi * sum;
    return res.convert_to<double>();
  }else{
    mp::float128 sum;
    int i;
    for(i=1; i<nu-1; i+=2){
      sum += M[i]+H[i];
    }
    return C + 2*sum.convert_to<double>();
  }
}
