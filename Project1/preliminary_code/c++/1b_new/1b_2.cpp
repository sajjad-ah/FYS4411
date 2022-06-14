
#include <iostream>
#include <cmath> //exp
#include <cstdlib>
#include <fstream>
#include <random>
#include <ctime>

using namespace std;




//create matrix
double **mat(int N, int M){

  double **A = new double * [N];
  for (int i=0; i<N; i++){
    A[i] = new double [M];
    }

  return A;
  }


//Wavefunction for non-interacting bosons
double waveFunction (double** R, int N, int dim, double alpha){

  double sum_squaredPos = 0.0;
  for (int i=0;i<N;i++){
    for (int j=0;j<dim;j++){
      sum_squaredPos += R[i][j]*R[i][j];
      }
    }
  return exp(-alpha*sum_squaredPos);
  }


double localEnergy(double** R,int N,int dim, double alpha){

  double sum_squaredPos = 0.0;
  for (int i=0; i<N;i++){
    for (int j=0;j<dim;j++){
      sum_squaredPos += R[i][j]*R[i][j];
      }
    }

  return dim*N*alpha+(0.5-2*alpha*alpha)*sum_squaredPos;
  }

double numericalLocalEnergy(double psi,double** R,int N,int dim, double alpha){

  double    h = 1.E-4;
  double    KE_sum , sum_squaredPos;
  double**  hPlus; double** hMinus;
  double    psiPlus, psiMinus;

  hPlus  = mat(N,dim);
  hMinus = mat(N,dim);

      //Copy position to hPlus and hMinus
      for (int n=0;n<N;n++){
	for (int m=0;m<dim;m++){
	  hPlus[n][m] = hMinus[n][m] = R[n][m];
	  }
	}

  //Kinetic energy part
  KE_sum = 0.;
  for (int i=0;i<N;i++){
    for (int j=0;j<dim;j++){



      hPlus[i][j]  += h;
      hMinus[i][j] -= h;


      psiPlus  = waveFunction(hPlus,N,dim,alpha);
      psiMinus = waveFunction(hMinus,N,dim,alpha);

      KE_sum += (psiPlus+psiMinus-2.0*psi)/(h*h);

      hPlus[i][j]  = R[i][j];
      hMinus[i][j] = R[i][j];


      }
    }


  delete hPlus; delete hMinus;

  //Harmonic oscillator part
  sum_squaredPos = 0.0;
  for (int i=0; i<N;i++){
    for (int j=0;j<dim;j++){
      sum_squaredPos += R[i][j]*R[i][j];
      }
    }


  return -0.5*KE_sum/(psi)+0.5*sum_squaredPos;
  }



void MonteCarloSampling(int M, int N, int dim, double stepSize, double alpha, mt19937_64 &randNum, double &energySample,double &energySquaredSample){

  //Declare variables
  double WF,testWF,relWF;
  double P;
  double E;

  double timeInitial;
  double timeFinal;

  //Mersenne twister
  double randomNumber;
  uniform_real_distribution<double> RNG(0,1);

  //Declare matrices for positions
  double** R      = mat(N,dim);
  double** testR  = mat(N,dim);

  energySample = 0.0; energySquaredSample = 0.0;

  //Initial position
  for (int i = 0; i<N;i++)
  {
     for (int j=0; j<dim;j++)
     {
        randomNumber = RNG(randNum);
        testR[i][j] = R[i][j] = stepSize*(randomNumber-0.5);
     }
  }

   WF = waveFunction(R,N,dim,alpha); //Initial wavefunction

   //Number of equilibration steps
   int equilibration  = 100;

   //Equlibration
   for (int eq=0;eq<equilibration;eq++){
     for (int i=0;i<N;i++){
       //Displace one particle
       for (int j=0;j<dim;j++){
	  randomNumber = RNG(randNum);
	  testR[i][j] = R[i][j] + stepSize*(randomNumber-0.5);
	 }

	//Evaluate wavefunction at test position
        testWF = waveFunction(testR,N,dim,alpha);

        relWF = testWF/WF;
	P     = relWF*relWF;
	//Metropolis test
	if (RNG(randNum) < P ){
	  for (int j=0;j<dim;j++){
	    R[i][j] = testR[i][j];
	    }
	  WF = testWF;
	  }else{
	    for (int j=0;j<dim;j++){
	      testR[i][j] = R[i][j];
	      }
	    }
        }
    }



   //Initial energy after equilibration
  //  E = localEnergy(R,N,dim,alpha); //Analytic local energy
   E = numericalLocalEnergy(WF,R,N,dim,alpha); //Numeric local energy

   int counter_print = 0 ;

   timeInitial = clock();
   //Start Monte Carlo Sampling
   for (int m=0;m<M;m++){
     for (int i=0;i<N;i++){
       //Displace one particle
       for (int j=0;j<dim;j++){
	  randomNumber = RNG(randNum);
	  testR[i][j] = R[i][j] + stepSize*(randomNumber-0.5);
	 }

	//Evaluate wavefunction at test position
        testWF = waveFunction(testR,N,dim,alpha);

        relWF = testWF/WF;
	P     = relWF*relWF;
	//Metropolis test
	if (RNG(randNum) < P ){
	  for (int j=0;j<dim;j++){
	    R[i][j] = testR[i][j];
	    }
	  //E = localEnergy(R,N,dim,alpha);
	  WF = testWF;
          E = numericalLocalEnergy(WF,R,N,dim,alpha);
	  }else{
	    for (int j=0;j<dim;j++){
	      testR[i][j] = R[i][j];
	      }
	    }

     }

     energySample        += E;
     energySquaredSample += E*E;

     //if (m % 50 == 0){ cout << m << endl;}
     cout << m << endl;
  }

  energySample        /= M;
  energySquaredSample /= M;

  timeFinal = clock();

  cout << "CPU time: " << timeFinal-timeInitial << endl;


  delete R; delete testR;

}

int main(){

  //Variables for random number
  double setSeed = clock();
  mt19937_64 randNum;
  randNum.seed(setSeed);

  //Variables for MC
  int M;
  int N;
  int dim;
  int a = 10;
  //double alpha;
  double alpha [a];
  //double alphaValues [8];
  double stepSize = 1.0;
  double energySample, energySquaredSample;
  double exactEnergy;

  cout << "Enter number of Monte Carlo Samples: ";
  cin  >> M;
  cout << "Enter number of particles: ";
  cin  >> N;
  cout << "Enter number of dimension(s): ";
  cin  >> dim;
  //cout << "Enter value of alpha: ";
  //cin  >> alpha;

  exactEnergy = 0.5*N*dim;

  alpha[0] = 0.1; alpha[1] = 0.2; alpha[2]  = 0.3; alpha[3] = 0.4;
  alpha[4] = 0.5; alpha[5] = 0.6; alpha[6]  = 0.7; alpha[7] = 0.8;
  alpha[8] = 0.9; alpha[9] = 1.0; alpha[10] = 1.1;


   const char* outputFileName = "data.dat";
   ofstream my_out(outputFileName);


for (int i=4; i<5;i++){

   //timeInitial = clock();
   MonteCarloSampling(M,N,dim,stepSize,alpha[i],randNum,energySample,energySquaredSample);
   //timeFinal = clock();


   cout << "alpha: " << alpha[i] << endl;
   cout << "Average energy: " << energySample << endl;
   cout << "exact energy: "   << exactEnergy << endl;
   cout << "Variance: " ;
   cout << energySquaredSample - energySample*energySample << endl;
  // cout << timeFinal - timeInitial << endl;
   cout << "************" << endl;

   my_out << alpha[i] << " " << energySample << " " << exactEnergy << endl;
}

  system("python plot.py data.dat");

  return 0;
  }
