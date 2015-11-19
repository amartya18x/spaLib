#include <stdlib.h>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;
class COSAMP
{
 public:
  long long int m, n, s,l, p;
  double learning_rate,decreThresh,prodThresh;
  COSAMP(long long int m, long long int n, long long int s,long long int l, long long int prune_param); //initialize
  //Encoder matrix is mxn and the vector is s-sparse
  void initValues(char *fx, char *fy, char *ftheta);//Initial Values
  mat calculateProxy(mat residue);//Gradient Update Step
  vector<int> findTopL(mat r_curr);//Find next atom
  void augment(vector <int> lambda);//Augemnt the atom matrix
  mat solveLeastSquare(mat A, mat b);//Solve the leastSquareproblem
  void sparse_multiply(double *A, double*b, double *c);//Multiply matrix with vector
  void prune(mat temp);//threshold back to s-sparsity
  void prune(mat temp, int p);//threshold back to s-sparsity
  void train();//start training
  mat returny(){return mat(this->y);}
  mat returnX(){return mat(this->X);}
  double return_error();//return the error
  double return_error(mat x);//return the error
  void printAtoms();
  double getBesterr(){return this->lastErr;}
  int stopCriterion();
  mat  X; //Stores the current X
  sp_mat orig; //Original Vector
 private:
  double lastErr;
  vector<int> permLambda ;
  mat Lambda; //List of non-zero atoms
  mat phiLambda; //Matrix with only relevant atoms as column_index
  mat phi; //Compressor matrix
  mat phiT; //Transpose of phi
  mat y; //Transmitted signal
  int atomCount;
  double best;
};
