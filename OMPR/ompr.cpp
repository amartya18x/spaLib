#include "ompr.h"
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
FILE * logfile;
using namespace std;
using namespace arma;
int patience = 5;
double preLast=-1;

typedef unsigned long long timestamp_t;

    static timestamp_t
    get_timestamp ()
    {
      struct timeval now;
      gettimeofday (&now, NULL);
      return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
    }


OMPR::OMPR(long long int m, long long int n, long long int s, long long int l,double learning_rate)
{
  this->m = m;
  this->n = n;
  this->s = s;
  this->l = l;
  this->learning_rate = learning_rate;
  //printf("m :[%d] n:[%d] s:[%d] l:[%d]  learning_rate:[%lf] \n",m,n,s,l,learning_rate);
}
void OMPR::initValues(char* fX, char *fy, char *ftheta)
{
  this->lastErr = 9999999.99;
  mat denseX;
  this->phi.load(fX);
  denseX.load(ftheta);
  this->y.load(fy);
  this->orig =sp_mat(denseX);
  this->s = this->orig.n_nonzero;
  this->l = this->s/10;
  this->X = sp_mat(size(denseX));
  this->phiT = this->phi.t();
  //this->orig.print();
  this->phiLambda = mat(size(phi));
  this->Lambda = sp_mat(size(denseX));
  this->atomCount = 0;
  this->permLambda.resize(this->phi.n_cols);
  for(int i=0; i<this->phi.n_cols;i++)
  {
    this->permLambda[i] = 0;
  }
  for(int i=0;i<this->s;i++)
    {
      //this->X(i) = 1;
      this->permLambda[i] = 1;
      for(int j=0;j<this->phi.n_rows;j++)
        {
          this->phiLambda(j,i)= this->phi(j,i);
        }
    }
  // printf("Number of rows in X: [%d]\n",this->X.n_rows);
  // printf("Number of cols in X: [%d]\n",this->X.n_cols);
  // printf("Number of rows in y: [%d]\n",this->y.n_rows);
  // printf("Number of cols in y: [%d]\n",this->y.n_cols);
  // printf("Number of rows in phi: [%d]\n",this->phi.n_rows);
  //printf("Number of cols in phi: [%d]\n",this->phi.n_cols);
}
int OMPR::stopCriterion()
{
  double currErr = return_error();
  if(currErr>2*this->lastErr || currErr < 0.1 ) //If twice the laswt error or reached convergence exit
    return 1;
 if(patience ==0 || preLast == currErr)
                return 1;
        if(abs(currErr - this->lastErr) <0.5)
                patience --;

  if(abs(currErr-this->lastErr)<this->decreThresh || currErr > 1.2*this->lastErr) // If decrease is very less or error has started increasing change threshold
    {
      if(this->l >1 &&  this->learning_rate > 0.05)
        {
	        this->l = this->l - 1;
          this->learning_rate = 0.5 * this->learning_rate;
          this->best = currErr;
        }
      else
        {
          return 1;
        }
    }
  else
    {
      preLast = this->lastErr;
      this->lastErr = return_error();
      this->best = this->lastErr;
      return 0;
    }
}
double OMPR::return_error()
{
  mat temp_X = this->orig - mat(this->X);
  double error = 0;
  for(int i=0;i<this->orig.n_rows;i++)
    error += temp_X(i)*temp_X(i);
  return error;
}
mat OMPR::updateStep()
{
  mat temp =  this->learning_rate*this->phiT*((mat)(this->y.col(1)) - this->phi*(this->X));    
  mat z = this->X + temp;
  return z;
}
double OMPR::return_error(mat X)
{
  mat temp_X = this->orig - X;
  double error = 0;
  for(int i=0;i<this->orig.n_rows;i++)
    error += temp_X(i)*temp_X(i);
  return error;
}
vector <int> OMPR::findTopL(mat r_curr)
{
  uvec max_find = sort_index(abs(r_curr),"descend");
  vector <int> indices ;
  int count = 0;
  for(int i=0; i<this->phi.n_cols && count <this->l; i++)
    {
      if(this->permLambda[(int)max_find(i)]!=1)
        {
          count ++;
          indices.push_back((int)max_find(i));
        }
      }
  return indices;
}
void OMPR::augment(vector <int> lambda)
{
  for(int j=0; j<this->l;j++)
    {
      this->permLambda[lambda[j]] = 1;
      for(int i=0; i<this->phi.n_rows; i++)
        {
          this->phiLambda(i,lambda[j]) = this->phi(i,lambda[j]);
        }
    }
}


mat OMPR::solveLeastSquare(mat A, mat b)
{
  mat x;
  x = pinv(A)*b;
  return x;
  }
void OMPR::threshold(mat temp)
{
  uvec max_find = sort_index(abs(temp),"descend");
  for(int j=this->s;j<this->s+this->l;j++)
    {
      temp((int)max_find(j))=0;
      this->permLambda[(int)max_find(j)] = 0;
      for(int i=0; i<this->phi.n_rows; i++)
        {
          this->phiLambda(i,(int)max_find(j)) = 0;
        }
    }
}
void OMPR::train()
{
  mat r = mat(this->y.col(1));
  mat x_i = mat(size(this->X));
while(!(this->stopCriterion()))
    {
      vector <int> lambda;
      mat r = this->updateStep();
      lambda = this->findTopL(r);
      this->augment(lambda);
      this->atomCount++;
      this->threshold(r);
      x_i = this->solveLeastSquare(this->phiLambda,(mat)this->y.col(1));
      this->X=sp_mat(x_i);
 //printf("Error:[%lf]\n",return_error()); 
    }
    return ;
}
int main(int argc, char* argv[])
{
  char foldername[100]="/home/amartya/data/samp/";
  //printf("%s\n",name.c_str());
  OMPR *solver = new OMPR(100,100,25,5,0.2);
  solver->initValues(argv[1],argv[2],argv[3]);
  timestamp_t t0, t1;
  t0 = get_timestamp();
  solver->train();
  t1 = get_timestamp();
  double time_taken = (t1 - t0)/1000000.0; // in seconds
  //printf("Time taken is %lf sec\n",time_taken);
  solver->returnX().save("result",raw_ascii);
  //printf("Error returned is %lf\n",solver->getBesterr());
  printf("OMPR:%s:%d:%d:%lf:%lf\n",argv[1],solver->returnX().n_rows,solver->returny().n_rows,time_taken,solver->getBesterr());
  return 0;
}

