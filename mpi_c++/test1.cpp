#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <iostream>
namespace mpi = boost::mpi;

//mpicxx test1.cpp -I/home/james/Documents/Apps/boost_1_65_1/stage/lib/ -L/home/james/Documents/Apps/boost_1_65_1/stage/lib/ -lboost_mpi -lboost_serialization -o test1

int main()
{
  mpi::environment env;
  mpi::communicator world;
  std::cout << "I am process " << world.rank() << " of " << world.size()
            << "." << std::endl;
  return 0;
}
