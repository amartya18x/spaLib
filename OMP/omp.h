#include <stdlib.h>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;
class OMP
{
 public:
  long long int m, n, s;
  OMP(long long int m, long long int n, long long int s); //initialize
  //Encoder matrix is mxn and the vector is s-sparse
  void dictFind(); //Find next atom
  void initValues(char *fX, char *fy, char *ftheta);//Initial Values
  int findBestProj(mat r_curr);//Find next atom
  void augment(int lambda);//Augemnt the atom matrix
  mat solveLeastSquare(mat A, mat b);//Solve the leastSquareproblem
  void sparse_multiply(double *A, double*b, double *c);//Multiply matrix with vector
  void train();//start training
  mat returnX(){return mat(this->X);}
   mat returny(){return mat(this->y);}
 double return_error();//return the error
  double return_error(mat x);//return the error
 private:
  double lastErr;
  mat Lambda;
 mat phiLambda;
  mat phi;
  mat phiT;
  sp_mat  X;
  mat y;
  sp_mat orig;
  int atomCount;
};
