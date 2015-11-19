#include <stdlib.h>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;
class OMPR
{
 public:
  long long int m, n, s,l;
  double decreThresh;
  double learning_rate;
  OMPR(long long int m, long long int n, long long int s,long long int l, double learning_rate); //initialize
  //Encoder matrix is mxn and the vector is s-sparse
  void initValues(char *fX, char *fy, char *ftheta);//Initial Values
  mat updateStep();//Gradient Update Step
  vector<int> findTopL(mat r_curr);//Find next atom
  void augment(vector <int> lambda);//Augemnt the atom matrix
  mat solveLeastSquare(mat A, mat b);//Solve the leastSquareproblem
  void sparse_multiply(double *A, double*b, double *c);//Multiply matrix with vector
  void threshold(mat temp);//threshold back to s-sparsity
  void train();//start training
  mat returnX(){return mat(this->X);}
  mat returny(){return mat(this->y);}
  double return_error();//return the error
  double return_error(mat x);//return the error
  int stopCriterion();
  double getBesterr(){return this->best;}
 private:
  vector<int> permLambda ;
  double best;
  mat Lambda; //List of non-zero atoms
  mat phiLambda; //Matrix with only relevant atoms as column_index
  mat phi; //Compressor matrix
  mat phiT; //Transpose of phi
  mat  X; //Stores the current X
  mat y; //Transmitted signal
  sp_mat orig; //Original Vector
  int atomCount;
  double lastErr;
};
