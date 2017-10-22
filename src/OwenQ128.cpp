// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <cmath>
#include <boost/math/distributions/non_central_t.hpp>

namespace mp = boost::multiprecision;
namespace sf = boost::math;
using namespace Rcpp;

// [[Rcpp::export]]
double pt_boost(double q, double nu, double delta){
  return boost::math::cdf(boost::math::non_central_t(nu, delta), q);
}

// [[Rcpp::export]]
NumericVector vecstudent(double q, int nu, NumericVector delta){
  const double a = R::sign(q)*sqrt(q*q/nu);
  const double b = nu/(nu+q*q);
  const double sB = sqrt(b);
  const int J = delta.size();
  if(nu==1){
    NumericVector C = pnorm(-delta*sB);
    int i;
    for(i=0; i<J; i++){
      C[i] += 2*sf::owens_t(delta[i] * sB, a);
    }
    return C;
  }
  std::vector<double> dsB(J); 
  int j;
  for(j=0; j<J; j++){
    dsB[j] = delta[j] * sB;
  }
  std::vector< std::vector<double> > M(nu-1, std::vector<double>(J));
  for(j=0; j<J ; j++){
    M[0][j] = a * sB * R::dnorm(dsB[j], 0.0, 1.0, 0) * R::pnorm(a*dsB[j], 0.0, 1.0, 1, 0);
  }
  const double sqrt2pi = 2.506628274631000502415765284811;
  if(nu>2){
    for(j=0; j<J; j++){
      M[1][j] = b * (delta[j] * a * M[0][j] + a * R::dnorm(delta[j], 0.0, 1.0, 0) / sqrt2pi);
    }
    if(nu>3){
      std::vector<double> A(nu-3);
      A[0] = 1.0;
      int k;
      if(nu>4){
        for(k=1; k<nu-3; k++){
          A[k] = 1.0/k/A[k-1];
        }
      }
      for(k=2; k<nu-1; k++){
        for(j=0; j<J; j++){
          M[k][j] = (k-1) * b * (A[k-2] * delta[j] * a * M[k-1][j] + M[k-2][j]) / k;
        }
      }
    }
  }
  int i;
  if(nu%2==1){
    NumericVector C = pnorm(-delta*sB);
    int i;
    for(i=0; i<J; i++){
      C[i] += 2*sf::owens_t(dsB[i], a);
    }
    std::vector<double> sum(J);
    for(i=1; i<nu-1; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j];
      }
    }
    NumericVector out(J);
    for(j=0; j<J; j++){
      out[j] = C[j] + 2*sum[j];
    }
    return out;
  }
  std::vector<double> sum(J);
  for(i=0; i<nu-1; i+=2){
    for(j=0; j<J; j++){
      sum[j] += M[i][j];
    }
  }
  NumericVector out(J);
  for(j=0; j<J; j++){
    out[j] = R::pnorm(-delta[j], 0.0, 1.0, 1, 0) + sqrt2pi*sum[j];
  }
  return out;
}

NumericVector isPositive(NumericVector x){
  int n = x.size();
  NumericVector out(n);
  int i;
  for(i=0; i<n; i++){
    out[i] = x[i] >= 0;
  }
  return out;
}
NumericVector isPositive(NumericVector x);

// [[Rcpp::export]]
NumericVector vecOwenQ1(int nu, double t, NumericVector delta, NumericVector R){
  const double a = R::sign(t)*sqrt(t*t/nu);
  const double b = nu/(nu+t*t);
  const double sB = sqrt(b);
  double ab;
  double asB;
  if(fabs(t)>DBL_MAX){
    ab = 0;
    asB = R::sign(t);
  }else{
    ab = a*b;
    asB = R::sign(t)*sqrt(t*t/(nu+t*t));
  }
  const int J = delta.size();
  if(nu==1){
    NumericVector C = pnorm(R) - isPositive(delta);
    int i;
    for(i=0; i<J; i++){
      double C1 =
        sf::owens_t(delta[i]*sB, a);
      double C2 =
        sf::owens_t(R[i], a-delta[i]/R[i]);
      double C3 =
        sf::owens_t(delta[i]*sB, (ab-R[i]/delta[i])/b);
      C[i] += 2*(C1 - C2 - C3);
    }
    return C;
  }
  std::vector<double> dsB(J); 
  std::vector<double> dnormdsB(J); 
  std::vector<double> dabminusRoversB(J); 
  std::vector<double> dnormR(J); 
  int j;
  for(j=0; j<J; j++){
    dsB[j] = delta[j] * sB;
    dnormdsB[j] = R::dnorm(dsB[j], 0.0, 1.0, 0);
    dabminusRoversB[j] = (delta[j]*ab - R[j])/sB; 
    dnormR[j] = R::dnorm(R[j], 0.0, 1.0, 0);
  }
  const int n = nu-1;
  std::vector< std::vector<double> > M(n, std::vector<double>(J));
  std::vector< std::vector<double> > H(n, std::vector<double>(J));
  for(j=0; j<J; j++){
    H[0][j] = -dnormR[j] * R::pnorm(a*R[j]-delta[j], 0.0, 1.0, 1, 0);
    M[0][j] = asB * dnormdsB[j] * (R::pnorm(dsB[j]*a, 0.0, 1.0, 1, 0) - R::pnorm(dabminusRoversB[j], 0.0, 1.0, 1, 0));
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = R[j] * H[0][j];
      M[1][j] = delta[j]*ab*M[0][j] + ab * dnormdsB[j] * (R::dnorm(dsB[j]*a, 0.0, 1.0, 0) - R::dnorm(dabminusRoversB[j], 0.0, 1.0, 0));
    }
    if(nu >= 4){
      std::vector<double> A(n);
      std::vector< std::vector<double> > L(n-2, std::vector<double>(J));
      A[0] = 1;
      A[1] = 1;
      for(j=0; j<J; j++){
        L[0][j] = ab * R[j] * dnormR[j] * R::dnorm(a*R[j]-delta[j], 0.0, 1.0, 0)/2;
      }
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          for(j=0; j<J; j++){
            L[k][j] = A[k+2] * R[j] * L[k-1][j];
          }
        }
      }
      for(k=2; k<n; k++){
        for(j=0; j<J; j++){
          H[k][j] = A[k] * R[j] * H[k-1][j];
          M[k][j] = (k-1.0)/k * (A[k-2] * delta[j] * ab * M[k-1][j] + b*M[k-2][j]) - L[k-2][j];
        }
      }
    }
  }
  if(nu % 2 == 0){
    const double sqrt2pi = 2.506628274631000502415765284811;
    int i;
    std::vector<double> sum(J);
    for(i=0; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j]+H[i][j];
      }
    }
    NumericVector out(J);
    for(j=0; j<J; j++){
      out[j] = R::pnorm(-delta[j], 0.0, 1.0, 1, 0) + sqrt2pi*sum[j];
    }
    return out;
  }else{
    std::vector<double> sum(J);
    int i;
    for(i=1; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j]+H[i][j];
      }
    }
    NumericVector out(J);
    NumericVector C = pnorm(R) - isPositive(delta);
    for(j=0; j<J; j++){
      double C1 =
        sf::owens_t(dsB[j], a);
      double C2 =
        sf::owens_t(R[j], a-delta[j]/R[j]);
      double C3 =
        sf::owens_t(dsB[j], (ab-R[j]/delta[j])/b);
      C[j] += 2*(C1 - C2 - C3);
      out[j] = C[j] + 2*sum[j];
    }
    return out;
  }
}

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

mp::float128 dnorm_mp(mp::float128 x){
  return mp::exp(-x*x/2)/2.5066282746310005024157652848110452530069867406099Q;
}
mp::float128 dnorm_mp(mp::float128 x); 
// [[Rcpp::export]]
double powen4(int nu, double t1, double t2, double delta1, double delta2){
  double C;
  if(nu % 2 == 1){
    const double a1 = R::sign(t1)*sqrt(t1*t1/nu);
    const double b1 = nu/(nu+t1*t1);
    const double sB1 = sqrt(b1);
    const double ab1 = a1*b1;
    const double asB1 = R::sign(t1)*sqrt(t1*t1/(nu+t1*t1));
    const double a2 = R::sign(t2)*sqrt(t2*t2/nu);
    const double b2 = nu/(nu+t2*t2);
    const double sB2 = sqrt(b2);
    const double ab2 = a2*b2;
    const double asB2 = R::sign(t2)*sqrt(t2*t2/(nu+t2*t2));
    const double R = sqrt(nu)*(delta1 - delta2)/(t1-t2);
    C = (delta1>=0) - (delta2>=0) + 
      2*(sf::owens_t(delta2 *sB2, a2) -
          sf::owens_t(delta1 *sB1, a1) -
          sf::owens_t(R , (a2*R -delta2 )/R ) +
          sf::owens_t(R , (a1*R -delta1 )/R ) -
          sf::owens_t(delta2 *sB2, (delta2 *ab2-R )/b2/delta2) +
          sf::owens_t(delta1 *sB1, (delta1 *ab1-R )/b1/delta1) );
    if(nu == 1){
      return C;
    }
  }
  const mp::float128 a1(R::sign(t1)*sqrt(t1*t1/nu));
  const mp::float128 b1(nu/(nu+t1*t1));
  const mp::float128 sB1(sqrt(b1));
  const mp::float128 ab1(a1*b1);
  const mp::float128 asB1(R::sign(t1)*sqrt(t1*t1/(nu+t1*t1)));
  const mp::float128 a2(R::sign(t2)*sqrt(t2*t2/nu));
  const mp::float128 b2(nu/(nu+t2*t2));
  const mp::float128 sB2(sqrt(b2));
  const mp::float128 ab2(a2*b2);
  const mp::float128 asB2(R::sign(t2)*sqrt(t2*t2/(nu+t2*t2)));
  const mp::float128 R(sqrt(nu)*(delta1 - delta2)/(t1-t2));
  const int n = nu-1;
  std::vector<mp::float128> H(n);
  std::vector<mp::float128> M1(n);
  std::vector<mp::float128> M2(n);
  H[0] = -dnorm_mp(R) * (pnorm_mp(a2*R-delta2) - pnorm_mp(a1*R-delta1));
  M1[0] = asB1*dnorm_mp(delta1*sB1)*(pnorm_mp(delta1*asB1)-pnorm_mp((delta1*ab1-R)/sB1));
  M2[0] = asB2*dnorm_mp(delta2*sB2)*(pnorm_mp(delta2*asB2)-pnorm_mp((delta2*ab2-R)/sB2));
  if(nu >= 3){
    H[1] = R * H[0];
    M1[1] = delta1*ab1*M1[0] +
      ab1*dnorm_mp(delta1*sB1)*(dnorm_mp(delta1*asB1)-dnorm_mp((delta1*ab1-R)/sB1));
    M2[1] = delta2*ab2*M2[0] +
      ab2*dnorm_mp(delta2*sB2)*(dnorm_mp(delta2*asB2)-dnorm_mp((delta2*ab2-R)/sB2));
    if(nu >= 4){
      std::vector<mp::float128> A(n);
      std::vector<mp::float128> L1(n-2);
      std::vector<mp::float128> L2(n-2);
      A[0] = 1;
      A[1] = 1;
      L1[0] = ab1 * R * dnorm_mp(R) * dnorm_mp(a1*R-delta1) / 2;
      L2[0] = ab2 * R * dnorm_mp(R) * dnorm_mp(a2*R-delta2) / 2;
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          L1[k] = A[k+2] * R * L1[k-1];
          L2[k] = A[k+2] * R * L2[k-1];
        }
      }
      for(k=2; k<n; k++){
        H[k] = A[k] * R * H[k-1];
        M1[k] = (k-1.0)/k *
          (A[k-2] * delta1 * ab1 * M1[k-1] + b1*M1[k-2]) - L1[k-2];
        M2[k] = (k-1.0)/k *
          (A[k-2] * delta2 * ab2 * M2[k-1] + b2*M2[k-2]) - L2[k-2];
      }
    }
  }
  if(nu % 2 == 0){
    const mp::float128 sqrt2pi = 2.5066282746310005024157652848110452530069867406099Q;
    mp::float128 sum;
    int i;
    for(i=0; i<nu-1; i+=2){
      sum += M2[i]-M1[i]+H[i];
    }
    mp::float128 res = pnorm_mp(-delta2) - pnorm_mp(-delta1) + sqrt2pi * sum;
    return res.convert_to<double>();
  }else{
    mp::float128 sum;
    int i;
    for(i=1; i<nu-1; i+=2){
      sum += M2[i]-M1[i]+H[i];
    }
    return C+2*sum.convert_to<double>();
  }
}
