#include <boost/mpi.hpp>
#include <iostream>
#include <string>
#include <boost/serialization/string.hpp>
#include <Eigen/Dense>
#include <ctime>

using namespace Eigen;
using namespace std;
namespace mpi = boost::mpi;
//mpicxx test2.cpp -I/home/james/Documents/Apps/boost_1_65_1/stage/lib/ -I/home/james/Documents/Apps/eigen -L/home/james/Documents/Apps/boost_1_65_1/stage/lib/ -lboost_mpi -lboost_serialization -o test2


int main()
{

  std::clock_t start;
  double duration;

  start = std::clock();
  mpi::environment env;
  mpi::communicator world;

  int rank = world.rank();
  int commsize = world.size();
  int matrix_size = 9;

  MatrixXd A (matrix_size, matrix_size);
  MatrixXd B (matrix_size, matrix_size);
  MatrixXd C (matrix_size, matrix_size);

  for (int r=0; r<A.rows(); r++)
    for (int c=0; c<A.cols(); c++) {
      A(r,c) = 1;
      B(r,c) = 2;
    }



  cout << rank << endl;
  world.barrier();

  // if (rank == 1) {
  for (int i=0; i<C.rows(); i++) {
      if (i%commsize != rank) { continue; }
      for (int j=0; j<C.cols(); j++) {
        C(i,j) = 0;
        for (int k=0; k<A.rows(); k++) {
          C(i,j) += A(i,k)*B(k,j);
        }
      }
    }
  // }
  MPI_Allreduce(MPI_IN_PLACE, &C(0,0), C.rows()*C.cols(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (rank == 0) {
    if (C.rows() < 15 ) { cout << C << endl; }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    cout<<"printf: "<< duration <<'\n';
  }
  return 0;
}
