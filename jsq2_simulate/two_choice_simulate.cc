#include<iostream>
#include<random>

std::mt19937 mt_rand(0);
std::uniform_real_distribution<double> d(0.0,1.0); 

double rand_u(){
  return d(mt_rand);
}

   
class two_choice {
  double *X; // X[i] = number of queues with i or more jobs 
  int N; // total number of queues
  double rho; // system's load
  double oneOverN;
  int K; // queues capacities
  double q; // fraction of jobs that are served by a centralized server
  int withoutReplacement;  // to pick jobs without replacement 
  
 public:
  two_choice(int N, double rho, int K, double q, bool withoutReplacement)
    : N(N), rho(rho), K(K), q(q), withoutReplacement(withoutReplacement) {
    oneOverN = 1./(double)N;
    X = new double[K+1];
    for(int i=0;i<K+1;i++) X[i] = 0;
    X[0] = 1;
  }
  void print() {
    for(int i=0;i<K;i++) std::cout << X[i] << " ";
    std::cout << "\n";
  }
  void simulate(int T){
    if (q==0) simulate_twoChoice(T);
    else if (q==-1) simulate_oneChoice(T);
    else      simulate_powerCentralized(T);
  }
  void simulate_twoChoice(int T){
    if (withoutReplacement)
      simulate_twoChoice_withoutReplacement(T);
    else
      simulate_twoChoice_withReplacement(T);
  }
  void simulate_twoChoice_withReplacement(int T) {
    for(int t=0;t<T;t++){
      double u = rand_u();
      //std::cerr << u << "\n";
      if (u < 1/(1+rho) ) { // departure
	//std::cout << "departure\n";
	double u = rand_u();
	int i=0;
	//while( i<K && X[i]-X[i+1] < u) {u += X[i+1]-X[i]; i++;}
	while( i<K && u < X[i+1] ) {i++; }
	if (i>0) X[i]-=oneOverN;
      }
      else{                 // arrival
	//std::cout << "arrival @ ";
	double u = rand_u();
	
	int i=0;
	while( i<K && u < X[i+1]*X[i+1]) { i++;}
	if (i<K) X[i+1]+=oneOverN;
      }
    }
  }
   void simulate_oneChoice(int T) {
    for(int t=0;t<T;t++){
      double u = rand_u();
      if (u < 1/(1+rho) ) { // departure
	double u = rand_u();
	int i=0;
	while( i<K && u < X[i+1] ) {i++; }
	if (i>0) X[i]-=oneOverN;
      }
      else{                 // arrival
	//std::cout << "arrival @ ";
	double u = rand_u();
	double v = rand_u(); (void) v;
	int i=0;
	//while( i<K && X[i]-X[i+1] < u) {u += X[i+1]-X[i]; i++;}
	while( i<K && u < X[i+1] ) { i++; }
	if (i<K) X[i+1]+=oneOverN;
	//std::cout << i<< " (" << u<<")\n";
      }
    }
  }
  void simulate_twoChoice_withoutReplacement(int T) {
    for(int t=0;t<T;t++){
      double u = rand_u();
      if (u < 1/(1+rho) ) { // departure
	//std::cout << "departure\n";
	double u = rand_u();
	int i=0;
	while( i<K && X[i]-X[i+1] < u) {u += X[i+1]-X[i]; i++;}
	if (i>0) X[i]-=oneOverN;
      }
      else{                 // arrival
	//std::cout << "arrival @ ";
	double u = rand_u();
	int i=0;
	while( i<K && (X[i]-X[i+1])*(X[i+1]+X[i]+oneOverN)/(1-oneOverN) < u) {
	  u += -(X[i]-X[i+1])*(X[i+1]+X[i]+oneOverN)/(1-oneOverN); i++;}
	if (i<K) X[i+1]+=oneOverN;
	//std::cout << i<< " (" << u<<")\n";
      }
    }
  }
  void simulate_powerCentralized(int T) {
    for(int t=0;t<T;t++){
      double u = rand_u();
      if (u < 1/(1+rho) ) { // departure
	double u = rand_u();
	if (u<1-q) {
	  double u = rand_u();
	  int i=0;
	  while( i<K && X[i]-X[i+1] < u) {u += X[i+1]-X[i]; i++;}
	  if (i>0) X[i]-=oneOverN;
	} else {
	  int i=0;
	  while( i<K && X[i]>0){i++;}
	  if (i-1>0) X[i-1]-=oneOverN; 
	}
      }
      else{                 // arrival
	double u = rand_u();
	int i=0;
	while( i<K && X[i]-X[i+1] < u) {u += X[i+1]-X[i]; i++;}
	if (i<K) X[i+1]+=oneOverN;
	//std::cout << i<< " (" << u<<")\n";
      }
    }
  }
  void steady_state(){
    int startup_time = 10000*N; // sufficient for rho = .95 at least 
    simulate(startup_time);
    double *average = new double[20];
    for(int i=0;i<20;i++) average[i] = 0;
    int nb_samples = 20000;
    for(int t=0;t<nb_samples;t++){
      simulate(100);
      for(int i=0;i<30;i++) average[i] += X[i]/nb_samples;
    }
    for(int i=0;i<20;i++) std::cout << average[i] << " ";
    std::cout << std::endl;
  }
  void test_convergence(){
    simulate(1000*N); // seems to suffices for rho = 0.99
    for(int t=0;t<10000;t++){
      simulate(1000);
      print();
    }
  }
  void print_one_trajectory(){
    X[1]=1; X[2] = 1;
    for(int t=0;t<10*N;t++){
      simulate(1);
      print();
    }
  }
};


int main(int argc, char ** argv) {
  int N=1000;
  double rho=0.8;
  double q=0.0;
  int nb_experiments = -1;
  bool only_print_one_trajectory = false;
  bool withoutReplacement = false;
  for(int i=1;i<argc;i++){
    switch(* (argv[i]) ){
    case 'N': N=atoi(argv[i]+1); break;
    case 'r': rho=atof(argv[i]+1); break;
    case 'e': nb_experiments=atoi(argv[i]+1); break;
    case 'c': {two_choice simu(N,rho,40,q,withoutReplacement); simu.test_convergence(); exit(1);} break;
    case 't': only_print_one_trajectory=true; break;
    case 'q': q=atof(argv[i]+1); break;
    case 'W': withoutReplacement = true; break;
    }
  }
  time_t t;
  time(&t);
  mt_rand.seed(t+N);
  //mt_rand.seed(2);

  //std::cerr << q <<"\n";
  if (only_print_one_trajectory){
    two_choice simu(N,rho,30,q, withoutReplacement);
    simu.print_one_trajectory();
  }
  else
    {
      if (nb_experiments<0) nb_experiments = 5+(50 + N*N)/100;
      //std::cerr << N << " " << nb_experiments << "\n";
      for(int i=0;i< nb_experiments;i++) {
	two_choice simu(N,rho,40,q, withoutReplacement);
	simu.steady_state();
      }
    }
}
