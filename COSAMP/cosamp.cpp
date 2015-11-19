#include "cosamp.h"
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
using namespace std;
using namespace arma;

COSAMP::COSAMP(long long int m, long long int n, long long int s, long long int l,long long int prune_param)
{
  this->m = m;
  this->n = n;
  this->s = s;
  this->l = l;
  this->p = prune_param;
  //  printf("m :[%d] n:[%d] s:[%d] l:[%d]  learning_rate:[%lf] \n",m,n,s,l,learning_rate);
}
void COSAMP::initValues(char* fY, char *fX, char *ftheta)
{
  this->decreThresh = 0.5;
  this->prodThresh = 1.2;
  this->lastErr = 99999.999;
  mat denseX;
  this->phi.load(fX);
  denseX.load(ftheta);
  this->y.load(fY);
  this->orig =sp_mat(denseX);
  this->X = sp_mat(size(denseX));
  this->s = this->orig.n_nonzero;
  this->l = 2*this->s;
  this->p = this->l;
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
      this->X(i) = 1;
      this->permLambda[i] = 1;
      for(int j=0;j<this->phi.n_rows;j++)
        {
          this->phiLambda(j,i)= this->phi(j,i);
        }
    }
  //printf("Sparsity is %d\n", this->s);
  //printf("Number of rows in X: [%d]\n",this->X.n_rows);
  //printf("Number of cols in X: [%d]\n",this->X.n_cols);
  //printf("Number of rows in y: [%d]\n",this->y.n_rows);
  //printf("Number of cols in y: [%d]\n",this->y.n_cols);
  //printf("Number of rows in phi: [%d]\n",this->phi.n_rows);
  //printf("Number of cols in phi: [%d]\n",this->phi.n_cols);
}

double COSAMP::return_error()
{
  mat temp_X = this->orig - mat(this->X);
  double error = 0;
  for(int i=0;i<this->orig.n_rows;i++)
    error += temp_X(i)*temp_X(i);
  return error;
}
mat COSAMP::calculateProxy(mat residue)
{
  mat z =  this->phiT*residue;
  return z;
}
double COSAMP::return_error(mat X)
{
  mat temp_X = this->orig - X;
  double error = 0;
  for(int i=0;i<this->orig.n_rows;i++)
    error += temp_X(i)*temp_X(i);
  return error;
}
int COSAMP::stopCriterion()
{
  double currErr = return_error();
  if(currErr>1.3*this->lastErr || currErr < 0.1 ) //If twice the laswt error or reached convergence exit
    return 1;
  if(abs(currErr-this->lastErr)<this->decreThresh || currErr > prodThresh*this->lastErr) // If decrease is very less or error has started increasing change threshold
    {
      if(this->l >1 &&  this->learning_rate > 0.05)
        {
          
          this->best = currErr;
          this->l = this->l - 1;
          this->p = this->p - 1;
        }
      else
        {
          return 1;
        }
    }
  else
    {
      this->lastErr = return_error();
      this->best = this->lastErr;
      return 0;
    }
}
/*vector <int> COSAMP::findTopL(mat r_curr)
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
}*/
vector <int> COSAMP::findTopL(mat r_curr)
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
void COSAMP::augment(vector <int> lambda)
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


mat COSAMP::solveLeastSquare(mat A, mat b)
{
  mat x;
  x = pinv(A)*b;
  return x;
  }
void COSAMP::printAtoms()
{
     for(int i=0;i<this->n;i++)
      {
        if(this->permLambda[i]==1)
          {
            printf("(%d):[%lf]\n",i,this->X(i));
          }
      }
      printf("################################\n");
 }
void COSAMP::prune(mat temp)
{
  uvec max_find = sort_index(abs(temp),"descend");
  for(int j=this->s;j<this->s+this->p;j++)
    {
      this->X((int)max_find(j))=0;
      this->permLambda[(int)max_find(j)] = 0;
      for(int i=0; i<this->phi.n_rows; i++)
        {
          this->phiLambda(i,(int)max_find(j)) = 0;
        }
    }
}

void COSAMP::prune(mat temp, int pr)
{
  uvec max_find = sort_index(abs(temp),"descend");
  for(int j=this->s;j<this->s+this->l;j++)
    {
      this->X((int)max_find(j))=0;
      this->permLambda[(int)max_find(j)] = 0;
      for(int i=0; i<this->phi.n_rows; i++)
        {
          this->phiLambda(i,(int)max_find(j)) = 0;
        }
    }
}

void COSAMP::train()
{
  mat residue = mat(this->y.col(1));
  mat x_i = mat(size(this->X));
  double lastErr = 1000000000000;
  vector <int> lambda;
  while(!(stopCriterion()))
    {
      residue = this->calculateProxy(residue);
      lambda = this->findTopL(residue);
      //printf("Number of non-zero elements in X are %d\n", sp_mat(this->X).n_nonzero);
      this->augment(lambda);
      this->atomCount++;
      //printf("Number of non-zero elements in X are %d\n", sp_mat(this->X).n_nonzero);
      x_i = this->solveLeastSquare(this->phiLambda,(mat)this->y.col(1));
      sp_mat Xt=sp_mat(size(x_i));
      for(int i=0;i< permLambda.size();i++)
        {
          if(permLambda[i] == 1) Xt(i) = x_i(i);
        }
      this->X = Xt;
      //sp_mat(this->X).print();
      //printf("Number of non-zero elements in X are %d\n", sp_mat(this->X).n_nonzero);
      this->prune(x_i);
      //this->printAtoms();
      //printf("Number of non-zero elements in X are %d\n", sp_mat(this->X).n_nonzero);
      residue = this->y.col(1) - this->phi*residue;
      //printf("Error returned is %lf\n",this->return_error(this->X));
    }
    for(int i=0;i<this->n;i++)
      {
        if(this->permLambda[i]!=1)
          {
            this->X(i) = 0;
          }
      }
    //sp_mat(this->X).print();
    return ;
}
int main(int argc, char* argv[])
{
  string name = argv[1];
  char foldername[100]="/home/amartya/data/samp/";
  //printf("%s\n",name.c_str());
  COSAMP *solver = new COSAMP(1554,500,25,25,25);
  solver->initValues(argv[2], argv[1], argv[3]);
  clock_t t;
  t = clock();
  solver->train();
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
  //printf("Time taken is %lf sec\n",time_taken);
  solver->returnX().save("result",raw_ascii);
  //printf("Error is %lf \n", solver->getBesterr());
  int corrRow = 0;
  for(int i=0; i< solver->X.n_rows; i++)
    if(abs(solver->X(i)) > 0.1 && abs(solver->orig(i)) > 0.1)
      corrRow ++;
  printf("COSAMP:%s:%d:%d:%lf:%lf:%lf\n",argv[1],solver->returnX().n_rows,solver->returny().n_rows,time_taken,solver->getBesterr(),corrRow*1.0/solver->s);
return 0;
}

