/*
Test file for plain cholesky decomposition of matrix.
*/
// #include <string.h>
// #include <stdlib.h>
// #include <complex.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>


using namespace Eigen;

// Helper Functions //


// Primary Functions //
void ICCholesky ( MatrixXd& A, MatrixXd& L ) {
  int n = A.rows();
  assert( n == A.cols() );

  for (int i=0; i < n; i++) {
    for (int j=0; j < i + 1; j++) {

      if ( i == j ) {
        L(i,i) = A(i,i);
        for (int k = 0; k < j; k++) {
          L(i,i) -= L(j,k) * L(j,k);
        }
        L(i,i) = sqrt( L(i,i) );
      }

      else {
        L(i,j) = A(i,j);
        for (int k=0; k < j; k++) {
          L(i,j) -= L(i,k)*L(j,k);
        }
        L(i,j) = L(i,j)/L(j,j);
      }
    }
  }
  std::cout << "Cholesky Decomposition of A" << std::endl << L << std::endl;
  std::cout << std::endl;
  return;
}

void OCCholesky ( MatrixXd M, double d, double s, double tau, MatrixXd L ) {
  /*
  Using the algorithm in Aquilante et al. 2011. Linear Scaling Techniques in
  Computational Chemistry and Physics. Chapter 13: Cholesky Decomposition
  Techniques in Electronic Structure Theory.
  NOTE: There are some notation changes, but I tried to keep to the notation
  used in the reference above as much as possible.
  */
  int n = M.rows();
  assert( n == M.cols() );

  // Step 1
  MatrixXd D (n,0);
  for ( int p=0; p < n; p++)
    D(p,0) = M(p,p);
  double D_max = D.maxCoeff();

  // Step 2
  std::vector<int> Ll (0);
  for (int p=0; p<n; p++) {
    if ( d * sqrt(D_max*D(p,0)) > tau ) {
      Ll.push_back( p );
    }
  }

  // Step 3
  int nv = 0;

  // Step 4
  int i = 0;

  // Step 5
  while ( D_max > tau ) {
    //// a
    i++;
    //// b
    double D_min = std::max(s*D_max,tau);
    //// c
    std::vector<int> Q (0);
    for (int q=0; q<Ll.size(); q++) {
      if ( D(Ll[q], 0) > D_min ) {
        Q.push_back( q );
      }
    }
    //// d
    MatrixXd M_pq (Ll.size(),Q.size());
    for (int p=0; p<Ll.size(); p++)
      for (int q=0; Q.size(); q++) {
        M_pq(p,q) = M( Ll[p], Q[q] );
      }

    //// e
    MatrixXd delta_pq (Ll.size(),Q.size());

    //// f
    double Q_max = 0;
    int q_j = 0;
    for (int q=0; q<Q.size(); q++) {
      if ( D(Q[q], 0) > Q_max ) {
        q_j = q;
        Q_max = D(Q[q], 0);
      }
    }

    //// g
    int j = 0;

    //// h
    while ( j < Q.size() && Q_max > D_min ) {
      ////// i
      j++;
      int J = nv + j;

      ////// ii
      // ?
      ////// iii
      for (int p=0; p<Ll.size(); p++) {
        L(p,J) = delta_pq(p,q_j)/sqrt(Q_max); // TODO Not sure if q_j is correct here
      }

      ////// iv
      for (int p=0; p<Ll.size(); p++) {
        D(p,0) -= (L(p,J) * L(p,J));
        for (int q=0; q<Q.size(); q++) {
          delta_pq(p,q) -= L(p,J)*L(q,J);
        }
      }

      Q_max = 0;
      for (int q=0; q<Q.size(); q++) {
        if ( D(Q[q],0) > Q_max ){
          q_j = q;
          Q_max = D(Q[q],0);
        }
      }

      //// i
      nv += j;
      D_max = 0;
      for (int p=0; p<Ll.size(); p++) {
        if ( D(Ll[p],0) > D_max ) {
          D_max = D(Ll[p],0);
        }
      }

      //// k
      Ll.clear();
      for (int p=0; p<n; p++) {
        if ( d*sqrt(D_max*D(p,0)) > tau ) {
          Ll.push_back(p);
        }
      }
    }
  }


  return;
}

void TestFromEigen ( MatrixXd A ) {
  std::cout << "The matrix A is" << std::endl << A << std::endl;
  LLT<MatrixXd> lltOfA(A); // compute the Cholesky decomposition of A
  MatrixXd L = lltOfA.matrixL(); // retrieve factor L  in the decomposition
  // The previous two lines can also be written as "L = A.llt().matrixL()"
  std::cout << "The Cholesky factor L is" << std::endl << L << std::endl;
  std::cout << "To check this, let us compute L * L.transpose()" << std::endl;
  std::cout << L * L.transpose() << std::endl;
  std::cout << "This should equal the matrix A" << std::endl;
}

/*******************************************************************************
******************************** Main ******************************************
*******************************************************************************/
int main() {
  // Read in 2RDM
  char line[255];

  char const *fIn = "o2_2rdm.txt";
  size_t chkorbs = 8;
  // printf("%zi\n",chkorbs);

  FILE *fp = fopen( fIn, "r" );
  fgets(line, sizeof(line), fp);
  int norb = atoi( strtok(line, " ,\t\n") );
  assert (norb==chkorbs);
  // printf("%i\n", norb);

  MatrixXd m2(norb*norb,norb*norb);

  while ( fgets(line, sizeof(line), fp) != NULL ) {
    int i = atoi( strtok(line, " ,\t\n") );
    int j = atoi( strtok(NULL, " ,\t\n") );
    int k = atoi( strtok(NULL, " ,\t\n") );
    int l = atoi( strtok(NULL, " ,\t\n") );
    double val = atof( strtok(NULL, " ,\t\n") );
    m2(i*norb+j,k*norb+l) += val;
  }

  // Create 1RDM
  MatrixXd m1(norb,norb);
  m1 = MatrixXd::Zero(norb,norb);

  for (int i=0; i < norb; i++ ) {
    for (int j=0; j < norb; j++ ) {
      // m1(i,j) = 0;
      // std::cout << i << " " << j << std::endl;
      for (int k=0; k < norb; k++ ) {
        m1(i,j) += m2(i*norb+k,j*norb+k)/11;
      }
    }
  }

  double chkelec = 0;
  for (int i=0; i < norb; i++ ) {
    chkelec += m1(i,i);
  }

  // std::cout << chkelec << std::endl; //TODO

  // Cholesky Decomposition
  MatrixXd l(norb,norb);
  l = MatrixXd::Zero(norb,norb);

  // ICCholesky( m1, l );

  // std::ofstream oput ("CD_o2_1RDM.txt");
  // if ( oput.is_open() ) {
  //   for (int i=0; i<norb; i++) {
  //     for (int j=0; j<norb; j++) {
  //       oput << l(i,j) << "         ";
  //     }
  //     oput << std::endl;
  //   }
  // }

  // Testing
  MatrixXd A(3,3);
  A << 4,-1,2, -1,6,0, 2,0,5;
  std::cout << "Matrix to decompose: \n" << A << std::endl << std::endl;

  MatrixXd l2(3,3);
  ICCholesky(A,l2);
  std::cout << "Check that recombining returns original \n";
  std::cout << l2 * l2.transpose() << std::endl << std::endl;

  MatrixXd l3(3,3);
  OCCholesky(A,1,0,0,l3);

  std::cout << l3 << std::endl;

  // TestFromEigen(A);

  return 0;
}
