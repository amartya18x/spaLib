#include "omp.h"
#include <string>
#include <dirent.h>
#include "vector"
#include <iostream>
#include <random>
#include <assert.h>
#include <algorithm>
#include <cblas.h>
#include <chrono>
#include <stdlib.h>
#include<fstream>
#include <iostream>
#include <armadillo>
#include <string.h>
#include <stdlib.h>
#include<fstream>
#include <time.h>
#include <sys/time.h>
using namespace std;
using namespace arma;


typedef unsigned long long timestamp_t;

    static timestamp_t
    get_timestamp ()
    {
      struct timeval now;
      gettimeofday (&now, NULL);
      return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
    }

OMP::OMP(long long int m, long long int n, long long int s)
{
  this->m = m;
  this->n = n;
  this->s = s;
  //  printf("m :[%d] n:[%d] s:[%d]\n",m,n,s);

  
}
void OMP::initValues(char* fX, char*fy, char *ftheta)
{
  mat denseX;
  this->phi.load(fX);
  denseX.load(ftheta);
  this->y.load(fy);
  this->orig =sp_mat(denseX);
  this->s = this->orig.n_nonzero;
  this->X = sp_mat(size(denseX));
  this->phiT = this->phi.t();
  //this->orig.print();
  this->phiLambda = mat(size(phi));
  this->Lambda = sp_mat(size(denseX));
  this->atomCount = 0;
  //printf("Number of rows in X: [%d]\n",this->X.n_rows);
  //printf("Number of cols in X: [%d]\n",this->X.n_cols);
  //printf("Number of rows in y: [%d]\n",this->y.n_rows);
  //printf("Number of cols in y: [%d]\n",this->y.n_cols);
  //printf("Number of rows in phi: [%d]\n",this->phi.n_rows);
  //printf("Number of cols in phi: [%d]\n",this->phi.n_cols);
}
double OMP::return_error()
{
  mat temp_X = this->orig - mat(this->X);
  double error = 0;
  for(int i=0;i<this->orig.n_rows;i++)
    error += temp_X(i)*temp_X(i);
  return error;
}

double OMP::return_error(mat X)
{
  mat temp_X = this->orig - X;
  double error = 0;
  for(int i=0;i<this->orig.n_rows;i++)
    error += temp_X(i)*temp_X(i);
  return error;
}
int OMP::findBestProj(mat r_curr)
{
  vec max_find = abs(phiT*r_curr);
  //max_find.print();
  uword max_ind;
  max_find.max(max_ind);
  return int(max_ind);
}
void OMP::augment(int lambda)
{
  this->Lambda[lambda] = 1;
  int i =0;
  for(i=0; i<this->phi.n_rows; i++)
    {
      this->phiLambda(i,lambda) = this->phi(i,lambda);
    }
}


mat OMP::solveLeastSquare(mat A, mat b)
{
  mat x;
  x = pinv(A)*b;
  return x;
  }
void OMP::train()
{
  mat r = mat(this->y.col(1));
  int lambda = -1;
  mat x_i = mat(size(this->X));
	//printf("%d and %d\n",this->s, atomCount);
    while(atomCount< this->s)
    {
    lambda = this->findBestProj(r);
    //printf("%d is selected\n",lambda);
    this->augment(lambda);
    this->atomCount++;
    x_i = this->solveLeastSquare(this->phiLambda,(mat)this->y.col(1));
    //printf("Error returned is %lf\n",this->return_error(x_i));
    this->X = sp_mat(x_i);
    r = mat(this->y.col(1)) - phiLambda*mat(X);
    }
    return ;
}
int main(int argc, char* argv[])
{
  OMP *solver = new OMP(500,1554,25);
  solver->initValues(argv[1],argv[2],argv[3]);
  timestamp_t t_0, t_1;
  t_0 = get_timestamp();
  solver->train();
  t_1 = get_timestamp();
  double time_taken = (t_1 - t_0)/1000000.0; // in seconds
  //printf("Time taken is %lf sec\n",time_taken);
  solver->returnX().save("result",raw_ascii);
  //printf("Error returned is %lf\n",solver->return_error());
  printf("OMP:%s:%d:%d:%lf:%lf\n",argv[1],solver->returnX().n_rows,solver->returny().n_rows,time_taken,solver->return_error());
  return 0;
}

