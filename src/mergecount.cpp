// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
using namespace Rcpp;

// large floats
#include <boost/multiprecision/float128.hpp>
// for large integers
#include <boost/multiprecision/cpp_int.hpp>
using namespace boost::multiprecision;

// This is the merging function for the fast win ratio computation algorithm.

// [[Rcpp::export]]
List mergeCountWins(const IntegerMatrix& group1,
                   const IntegerMatrix& group2) {
  int n1 = group1.nrow();
  int n2 = group2.nrow();
  float128 win = 0;
  float128 loss = 0;
  int i = 0;
  int j = 0;

  int grp1_tie_ct = 1;
  int grp2_tie_ct = 1;

  int iflag = 0;
  int jflag = 0;

  // version with intra and intergroup ties
  while (i < n1 && j < n2) {
    while (group1(i,1) == group1(i+1,1) ) {
      grp1_tie_ct += 1;
      i++;}
    while (group2(j,1) == group2(j+1,1) ){
      grp2_tie_ct += 1;
      j++;}
    if (group1(i,1) >= group2(j,1)) {
      win += (float128) grp1_tie_ct*group2(j,0);
      iflag = 1;
      }
    if (group1(i,1) <= group2(j,1)) {
      loss += (float128) grp2_tie_ct*group1(i,0);
      jflag = 1;
    }
    if (iflag == 1){
      i++;
      iflag = 0;
      grp1_tie_ct = 1;}
    if (jflag == 1){
      j++;
      jflag = 0;
      grp2_tie_ct = 1;}
    }

  float128 winprop_mp = (float128) win / (float128) (n1 * n2);
  float128 lossprop_mp = (float128) loss / (float128) (n1 * n2);

  double winprop = winprop_mp.convert_to<double>();
  double lossprop = lossprop_mp.convert_to<double>();

  return(List::create(Named("win", winprop),
                      Named("loss", lossprop))
           );

}

/*** R
lambda_D=0.2  # base rate for death
lambda_H=1    # base rate for nonfatal event
tau=4         # total length of followup
lambda_L=0.01
kappa=1

library(gumbel)
N=1000
outcome=rgumbel(N,alpha=kappa,dim=2)
D=-log(outcome[,1])/lambda_D
H=-log(outcome[,2])/lambda_H
Ca=tau
C = Ca
death=ifelse(D<C,1,0)
deathtime=D
deathtime[!death] = 0
nonfatalevent=ifelse(H<C,1,0)
nonfataleventtime=H
nonfataleventtime[!nonfatalevent] = 0
dat1=data.frame(death, deathtime, nonfatalevent, nonfataleventtime)
outcome=rgumbel(N,alpha=kappa,dim=2)
D=-log(outcome[,1])/(lambda_D)
H=-log(outcome[,2])/(lambda_H)
Ca=tau
C = Ca
death=ifelse(D<C,1,0)
deathtime=D
deathtime[!death] = 0
nonfatalevent=ifelse(H<C,1,0)
nonfataleventtime=H
nonfataleventtime[!nonfatalevent] = 0
dat2=data.frame(death, deathtime, nonfatalevent, nonfataleventtime)
rank1 = data.table::frankv(dat1, order=c(-1L, 1L, -1L, 1L), ties.method="first")
rank2 = data.table::frankv(dat2, order=c(-1L, 1L, -1L, 1L), ties.method="first")

newdf1 = cbind(dat1, group=1, groupwins=rank1)
newdf2 = cbind(dat2, group=2, groupwins=rank2)
totaldf = rbind(newdf1, newdf2)
data.table::setorderv(totaldf, order = c(1L, -1L, 1L, -1L, 1L, -1L))
totaldf$key = dim(totaldf)[1]:1
datatopass1 = totaldf[totaldf$group == 1, c("groupwins", "key")]
datatopass2 = totaldf[totaldf$group == 2, c("groupwins", "key")]

test1 = as.matrix(datatopass1)
test2 = as.matrix(datatopass2)
mergeCountWins(test1, test2)
*/
