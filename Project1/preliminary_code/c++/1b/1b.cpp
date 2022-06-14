/*
Same 1b code as written in python, just translated to C++: Brute force non-interacting bosons VMC.
Feel free to improve and clean this code.

Comment: The variance rounds to zero, but it should probably be lower with number of MC samples M=100000,
perhaps change the variable type of the energy to something with more decimals or increase M.
*/


#include <iostream>
#include <cmath> //For pow() and exp()
#include <cstdlib> //Includes rand() function
#include <fstream> //Write to file
#include <stdlib.h> //Type commands to terminal

using namespace std;


//Create an NxM  matrix
double **mat(int N,int M){

  double ** A = new double * [N]; //Array with N elements
  for (int i=0; i<N; i++){
    A[i] = new double [M];

    }
    return A;
  }

//Create wave function
//Input: R; positions, N; number of particles, alpha;
double waveFunction(double** R,int N,double alpha){

  double r2 = 0;
  for (int i=0;i<N;i++){
    r2 += pow(R[i][0],2.0)+pow(R[i][1],2.0)+pow(R[i][2],2.0);
  }

  return exp(-alpha*r2);
}


//Local energy
double localEnergy(double** R,int N,int dim,double alpha,double omega){

  double r2 = 0;
  for (int i=0;i<N;i++){
    r2 += pow(R[i][0],2.0)+pow(R[i][1],2.0)+pow(R[i][2],2.0);
  }

  return dim * N * alpha - 2*pow(alpha,2.0) * r2 + 0.5 * pow(omega,2.0) * r2;
}


//Returns Nx3 matrix with the inital positionss of the particles
double** initialize(int N,int dim, double stepSize){

  double** mat(int N,int M) ; //First declare a function that creates a (N,dim) matrix with zeros
  double** R = mat(N,3); //Declare R to be this matrix
  for (int i=0;i<N;i++){
    for (int j=0;j<dim;j++){
      double randNum = rand()/double(RAND_MAX); //random number in element (0,1], same seed everytime
      R[i][j] = stepSize*(randNum-0.5); //Initial position
    }
  }

  return R;
}



//Monte Carlo sampling function
double* monteCarlo(int M,int N,int dim,double stepSize,double alpha,double omega){

  //Initialize functions needed for the Monte Carlo sampling
  double* energySamplesArray = new double [2]; //[0]: energySample, [1] energy2Sample

  energySamplesArray[0] = energySamplesArray[1] = 0 ;

  energySamplesArray[0] =  energySamplesArray[1] = 0; //Set the elements to zero
  
  double** initialize(int N,int dim, double stepSize);
  double** R = initialize(N,dim,stepSize); //matrix
  double localEnergy(double** R,int N,int dim,double alpha,double omega);
  double E = localEnergy(R,N,dim,alpha,omega);
  double waveFunction(double** R,int N, double alpha);
  double WF = waveFunction(R,N,alpha);

  //Start Monte Carlo sampling
  for (int m=0;m<M;m++){
    for (int i=0;i<N;i++){ //Iterate thorugh the particles, one particle at a time

	double* savePrevPos = new double [3]; //Array for saving previous position

	for(int u = 0;u<dim;u++){
		savePrevPos[u] = R[i][u]; //Saving previous position
		}

	//cout << savePrevPos[0] << endl;

        for (int j=0;j<dim;j++){ //Move one particle at a time
          double randNum = rand()/double(RAND_MAX);
           R[i][j] = R[i][j] + stepSize*(randNum-0.5); //Displacement
        }
	//cout << savePrevPos[0] << endl;


      double WF_Trial = waveFunction(R,N,alpha); //Trial wavefunction

      //Metropolis test
      if ((rand()/double(RAND_MAX)) < pow(WF_Trial/WF,2.0)){
        E = localEnergy(R,N,dim,alpha,omega);
        WF = WF_Trial;
        }
      else{
	for (int k = 0;k<dim;k++){
		R[i][k] = savePrevPos[k]; //Rejecting movement, return to previous position
	}

    }
    delete savePrevPos;
}
    //Sample energy here
    energySamplesArray[0] += E;
    energySamplesArray[1] += E*E;

}
  delete R;

  //Get average
  energySamplesArray[0] /= M;
  energySamplesArray[1] /= M;
  //double e = energySamplesArray[0];
  //double e2 = energySamplesArray[1];
  //delete energySamplesArray;

  //return e,e2;
  return energySamplesArray;

}

//Generates an array of alpha values
//a: start value
//n: number of values
//incr: how much the value increments each time
double* alphaValues(double a,int n, double incr){

  double* listAlpha = new double [n];
  listAlpha[0] = a;
  double j = a;
  for (int i=1;i<n;i++){
    j += incr;
    listAlpha[i] = j;
  }
  return listAlpha;
}


//First argument: Number of Monte Carlo cycles, M
//Second argument: Number of particles, N
//Third argument: Number of dimentions, dim
int main(int argc, char** argv){

  double* monteCarlo(int M,int N,int dim,double stepSize,double alpha,double omega);

  //DECLARE VARIABLES HERE
  int M = atoi(argv[1]); //1000; //Number of Monte Carlo samples
  int N = atoi(argv[2]);//5; //Number of particles
  int dim = atoi(argv[3]);//2; //Number of dimentions
  double stepSize = 1.; //Step size
  double omega =1.0; //Omega (frequency of oscillator)

  int n = 10; //Number of alpha values
  double a = 0.1;
  double incr = 0.1;
  double* alphaValues(double a,int n, double incr);
  double* alpha = alphaValues(a,n,incr);

  //Calculate exact energy
  double exactEnergy = dim*N*(omega/2);


  //Calculate and plot results
  cout << "************* START ***************" << endl;

  //Inputs
  cout << "Number of Monte Carlo samples: " << M << endl;
  cout << "Number of particles N: " << N << endl;
  cout << "Dimention(s): " << dim << endl;
  cout << "Omega: " << omega << endl;
  cout << "**************************" << endl;
  cout << "**************************" << endl;

  //Write to file data.dat
  const char* outputFileName = "data.dat"; //data.dat has to be in the same directory!
  ofstream my_out(outputFileName);
  for (int i=0;i<n;i++){
    double* MC = monteCarlo(M,N,dim,stepSize,alpha[i],omega);
    //averageEnergy[i] = MC[0];
    my_out << alpha[i] << " " << MC[0] << " " << exactEnergy << endl ;

    //Stats.
    cout << "Alpha: " << alpha[i] << endl;
    cout << "Average energy: " << MC[0] << endl;
    cout << "Variance: " << MC[1]-MC[0]*MC[0] << endl;
    cout << "Exact energy: " << exactEnergy << endl;
    cout << "*******************************" << endl;

    delete MC;
  }
  delete alpha;
  system("python plot.py data.dat"); //Plot with plot.py, python code

  return 0;
}
