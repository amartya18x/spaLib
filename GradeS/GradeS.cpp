#include "GradeS.h"
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
int patience = 10;

typedef unsigned long long timestamp_t;

    static timestamp_t
    get_timestamp ()
    {
      struct timeval now;
      gettimeofday (&now, NULL);
      return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
    }


GradeS::GradeS(long long int m, long long int n, long long int s,double gamma)
{
  this->m = m;
  this->n = n;
  this->s = s;
  this->gamma = gamma;
  //printf("m :[%d] n:[%d] s:[%d] gamma:[%lf]\n",m,n,s,gamma); 
}
void GradeS::initValues(char* filenameX, char*filenamey, char*filenametheta)
{
  
  this->decreThresh = 0.5;
  this->prodThresh = 1.05;
  this->lastErr = 9999999.99;
  mat denseX;
  this->phi.load(filenameX);
  denseX.load(filenametheta);
  this->y.load(filenamey);
  this->orig =sp_mat(denseX);
  
  this->s = this->orig.n_nonzero;
  this->X = sp_mat(size(denseX));
  this->phiT = this->phi.t();
  //this->orig.print();
  this->phiLambda = mat(size(phi));
  this->Lambda = sp_mat(size(denseX));
  //printf("Number of rows in X: [%d]\n",this->X.n_rows);
  //printf("Number of cols in X: [%d]\n",this->X.n_cols);
  //printf("Number of rows in y: [%d]\n",this->y.n_rows);
  //printf("Number of cols in y: [%d]\n",this->y.n_cols);
  //printf("Number of rows in phi: [%d]\n",this->phi.n_rows);
  //printf("Number of cols in phi: [%d]\n",this->phi.n_cols);
}
int GradeS::stopCriterion()
{
  double currErr = return_error();
  if(currErr>2*this->lastErr || currErr < 0.2 || (currErr < 1 && abs(lastErr - currErr)< 0.1) ) //If twice the laswt error or reached convergence exit
    return 1;
else
{
if(abs(currErr - lastErr)<0.001)
{
	if(patience == 0)
		return 1;
	patience --;
	this->gamma = this->gamma * 2;
}
      this->lastErr = return_error();
      this->best = this->lastErr;
      return 0;
}
}
double GradeS::return_error()
{
  mat temp_X = this->orig - mat(this->X);
  double error = 0;
  for(int i=0;i<this->orig.n_rows;i++)
    error += temp_X(i)*temp_X(i);
  return error;
}

double GradeS::return_error(mat X)
{
  mat temp_X = this->orig - X;
  double error = 0;
  for(int i=0;i<this->orig.n_rows;i++)
    error += temp_X(i)*temp_X(i);
  return error;
}
mat GradeS::update(mat x)
{
  mat temp =  this->phiT*((mat)(this->y.col(2)) - this->phi*(this->X))/this->gamma;    
  mat z = this->X + temp;
  return z;
}
sp_mat GradeS::threshold(mat temp)
{
  uvec max_find = sort_index(abs(temp),"descend");
  vector <int> indices ;
  int count = 0;
  sp_mat z = sp_mat(size(this->X));
  for(int i = 0;i<this->s;i++)
    {
      z((int)max_find(i)) = temp((int)max_find(i));
    }
  return z;
 
}
void GradeS::train()
{
  mat r = mat(this->y.col(1));
  int lambda = -1;
  mat x_i = mat(size(this->X));
  while(!(this->stopCriterion()))
    {
      mat temp = update((mat)this->X);
      this->X = threshold(temp);
  //    printf("Error returned is %lf\n",this->return_error());
    }
    return ;
}
int main(int argc, char* argv[])
{
  //char *filenamey = "/home/amartya/data/Samples_p500_s25_e0.1_o10_Ceps0.1/y";
  //char *filenametheta = "/home/amartya/data/Samples_p500_s25_e0.1_o10_Ceps0.1/X";
  //char *filenameX = "/home/amartya/data/Samples_p500_s25_e0.1_o10_Ceps0.1/theta";
  char *filenameX = argv[1];
  char *filenamey = argv[2];
  char *filenametheta = argv[3];
  string name = argv[1];
  char foldername[100]="/home/amartya/data/samp/";
  //printf("%s\n",name.c_str());
  GradeS *solver = new GradeS(100,100,25,2);
  solver->initValues(filenameX,filenamey,filenametheta);
  timestamp_t t0, t1;
  t0 = get_timestamp();
  solver->train();
  t1 = get_timestamp();
  double time_taken = (t1 - t0)/1000000.0; // in seconds
  //  printf("Time taken is %lf sec\n",time_taken);
  solver->returnX().save("result",raw_ascii);
  // mat temp_pr;
  // temp_pr.load("result");
  // sp_mat print_x = sp_mat(temp_pr);
  //print_x.print();
  //printf("Error returned is %lf\n",solver->getBesterr());
printf("GRADES:%s:%d:%d:%lf:%lf\n",argv[1],solver->returnX().n_rows,solver->returny().n_rows,time_taken,solver->getBesterr());
 return 0;
}

