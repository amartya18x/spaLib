#include <stdlib.h>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;
class GradeS
{
 public:
  long long int m, n, s;
  GradeS(long long int m, long long int n, long long int s,double gamma); //initialize
  //Encoder matrix is mxn and the vector is s-sparse
  void initValues(char *fX, char *fy, char *ftheta);//Initial Values
  void sparse_multiply(double *A, double*b, double *c);//Multiply matrix with vector
  sp_mat threshold(mat X);//returns a thresholded vector
  mat update(mat x);
  void train();//start training
  mat returnX(){return mat(this->X);}
  mat returny(){return mat(this->y);}
  double return_error();//return the error
  double return_error(mat x);//return the error
  double getBesterr(){return this->best;}
  int stopCriterion();
 private:
  double prodThresh;
  double decreThresh;
  double best;
  double lastErr;
  double gamma;
  mat Lambda;
 mat phiLambda;
  mat phi;
  mat phiT;
  sp_mat  X;
  mat y;
  sp_mat orig;
};
