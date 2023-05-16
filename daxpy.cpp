#include <iostream>
#include <vector>
#include <omp.h>
#include <cstdint>

using namespace std;

bool checkEqualityDaxpy(const vector<double>& v1, const vector<double>& v2) {
  if (v1.size() != v2.size()) {
    cerr << "Vectors are of unequal size." << endl;
    return false;
  }

  for (int i = 0; i < v1.size(); i++) {
    if ((abs(v1[i] - v2[i]) / v1[i]) > 0.005) {
      cerr << "Error at " << i << ":" << v1[i] << " vs. " << v2[i] << endl;
      return false;
    }
  }

  cout << "Correctness test passed." << endl;
  return true;
}

void printVector(const vector<double>& v) {
  if (v.size() > 25) {
    cout << "Vectors too long to print." << endl;
    return;
  }

  for (double value : v) {
    cout << value << " ";
  }
  cout << endl;
}

int fillVectorWithNValues(vector<double>& v, const int64_t n) {
	v.clear();
	v.resize(n);
	for (int i = 0; i < n; i++) {
		v[i] = i;
	}
	return 0;
}

/*
* y = y + a*x
*/
int daxpy_ser(vector<double>& y, const vector<double>& x, double a) {

  // iterate over all elements in y
  for (int i = 0; i < y.size(); i++){

    // update each value in the vector y
    y[i] = y[i] + x[i]*a;
  }
	return 0;
}

/*
* y = y + a*x
*/
int daxpy_par(vector<double>& y, const vector<double>& x, double a) {

  // d)
  // int i; // define index variable before the for loop to declare it as private in OpenMP construct
  // #pragma omp parallel for default(none) private(i) shared(y, x, a)
  // for (i = 0; i < y.size(); i++){
  //   // update each value in the vector y
  //   y[i] = y[i] + x[i]*a;
  // }

  // f)
  int i; // define index variable before the for loop to declare it as private in OpenMP construct
  #pragma omp parallel for default(none) private(i) shared(y, x, a) schedule(static,3)
  for (i = 0; i < y.size(); i++){
    // update each value in the vector y
    y[i] = y[i] + x[i]*a;
  }
  return 0;
}

int main(int argc, char** argv) {
  const int nt = argc > 1 ? atoi(argv[1]) : 1;

  // set the number of threads to the value of input parameter nt
  omp_set_num_threads(nt);

  // print out the maximum number of threads that OpenMP can use
  cout << "Number of used threads: " << omp_get_max_threads() << endl;

  // Initialization
  vector<double> x, y, y2;
  const double a  = 0.5;
  const int64_t n = 500000000;
  fillVectorWithNValues(x, n);
  fillVectorWithNValues(y, n);
  fillVectorWithNValues(y2, n);

  // Serial daxpy
  const auto t1 = omp_get_wtime();
  daxpy_ser(y, x, a);
  const auto t2 = omp_get_wtime();
  cout << "Serial daxpy: " << t2 - t1 << "s." << endl;

  // Parallel daxpy
  const auto t3 = omp_get_wtime();
  daxpy_par(y2, x, a);
  const auto t4 = omp_get_wtime();
  cout << "Parallel daxpy: " << t4 - t3 << "s." << endl;

  // Print results
  printVector(y);
  printVector(y2);

  if (checkEqualityDaxpy(y, y2)) {
    return 0;
  }
  return 1;
}
