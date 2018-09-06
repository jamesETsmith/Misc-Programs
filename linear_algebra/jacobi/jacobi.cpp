#include <math.h> /* sqrt */
#include <mpi.h>
#include <array>
#include <fstream>
#include <iostream>
#include <string>

const size_t N = 1000;

typedef std::array<std::array<double, N>, N> Mat2d;

// Helper function do not modify
void readMatrix(Mat2d& A) {
  std::string line;
  std::ifstream f("jacobi_matrix.txt");

  if (f.is_open()) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        f >> A[i][j];
      }
    }
    f.close();
  }
}

// Helper function do not modify
void readVector(std::array<double, N>& b, std::string name) {
  std::string line;
  std::ifstream f(name);

  if (f.is_open()) {
    for (int i = 0; i < N; i++) {
      f >> b[i];
    }
    f.close();
  }
}

// Parallelize this
std::array<double, N> jacobi(const size_t maxiter) {
  // Set up
  Mat2d A;
  std::array<double, N> b;
  std::array<double, N> ans;
  std::array<double, N> x_old = {0};
  readMatrix(A);
  readVector(b, "jacobi_vector.txt");
  readVector(ans, "jacobi_ans.txt");

  // Actual Jacobi algorithm
  for (int it = 0; it < maxiter; it++) {
    std::array<double, N> x_new = {0.0};  // Initialize with all 0.0
    double abs_error = 0;
    double abs_diff = 0;

    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        if (i != j) {
          x_new[i] -= A[i][j] * x_old[j];
        }
      }

      x_new[i] += b[i];
      x_new[i] /= A[i][i];

      // Copying and error calc.
      abs_diff += (x_new[i] - x_old[i]) * (x_new[i] - x_old[i]);
      abs_error += (x_new[i] - ans[i]) * (x_new[i] - ans[i]);
      x_old[i] = x_new[i];
      // std::cout << x_old[i] << " ";
    }
    // Print out error and difference from previous iterarion
    // std::cout << std::endl;

    abs_error = sqrt(abs_error);
    abs_diff = sqrt(abs_diff);
    if (it % 50 == 0) {
      std::cout << "Iteration " << it << " " << abs_error << " " << abs_diff
                << std::endl;
    }

    if (abs_diff < 1e-8) {
      std::cout << "CONVERGED" << std::endl;
      return x_old;
    }
  }
  std::cout << "NOT CONVERGED" << std::endl;
  return x_old;
}

int main() {
  jacobi(1000);

  return 0;
}
