/*
   Test file for plain cholesky decomposition of matrix.

   Compile with:
   g++ -std=c++11 -g cholesky.cpp -I/usr/local/include/eigen3/ -o cholesky
   or
   g++ -std=c++11 -g cholesky.cpp -I/home/james/Documents/Apps/eigen/ -o cholesky
   or
   g++ -std=c++11 cholesky.cpp -I/projects/jasm3285/eigen/ -o cholesky
 */
// #include <string.h>
#include <stdlib.h>
// #include <complex.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>

// using namespace std;
using namespace Eigen;

// Helper Functions
bool descend (int i,int j) {
	return (i>j);
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) { //TODO Correct spelling in name

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
	  [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

	return idx;
}

void reorderBasis (MatrixXd& L, std::vector<size_t>& idx) {
	MatrixXd lcopy = L.replicate(L.rows(),L.cols());

	for (int i=0; i<L.rows(); i++)
		for (int j=0; j<L.cols(); j++) {
			L(i,j) = lcopy(idx[i],idx[j]);
		}
	return;
}

void TestFromEigen ( MatrixXd& A ) {
	std::cout << "The matrix A is\n" << std::endl << A << std::endl;
	LLT<MatrixXd> lltOfA(A); // compute the Cholesky decomposition of A
	MatrixXd L = lltOfA.matrixL(); // retrieve factor L  in the decomposition
	// The previous two lines can also be written as "L = A.llt().matrixL()"
	std::cout << "\nThe Cholesky factor L is" << std::endl << L << std::endl;
	std::cout <<"\nTo check this, let us compute L * L.transpose() - A\n";
	std::cout << "\n"<<L * L.transpose() - A << std::endl;
	std::cout << "\nThis should equal the 0 matrix" << std::endl;
}


void checkDecomposition(MatrixXd& M, std::vector<std::vector<double> >& Lv,
  int nv) {

	MatrixXd L = MatrixXd::Zero(M.rows(),M.cols());
	for (int c=0; c<M.rows(); c++) {
		for (int r=c; r<M.rows(); r++) {
			L(r,c) = Lv[c][r-c];
		}
	}

	// std::cout.precision(2);

	// Eigenvalues to test that it's positive semidefinite
	EigenSolver<MatrixXd> es;
	es.compute(M, /* computeEigenvectors = */ false);
	bool semiDef = true;
	for (int i=0; i<M.cols(); i++) {
		if ( es.eigenvalues()[i].real() < 0 ) {semiDef = false; }
		if ( es.eigenvalues()[i].imag() < 0 ) {semiDef = false; }

	}

	if (!semiDef) std::cout << "Matrix isn't semidefinite!!"<<"\n\n";

	if ( M.rows()<21 ) {
		std::cout << "The eigenvalues of M are:\n";
		std::cout << es.eigenvalues().transpose() <<"\n\n";

		std::cout << "M\n\n" << M << "\n\n\n";
		std::cout << "L\n\n" << L << "\n\n\n";
	}

	//Comparison to original
	MatrixXd comp = L * L.transpose() - M;
	if (comp.rows() < 21) {
		std::cout << "Comparison of L*L^dag - M\n\n";
		std::cout << comp << "\n\n";
	}

	// Calculate Norm of Comparison
	double norm = 0;
	for (int r=0; r<L.rows(); r++)
		for (int c=0; c<L.cols(); c++)  {
			norm += comp(r,c)*comp(r,c);
		}
	norm = sqrt(norm);
	std::cout << "Norm of comparison: " <<norm << "\n";
	std::cout << "Number of rows truncated: " << M.rows() - nv << "\n\n";
}

void checkDecomposition(MatrixXd& M, MatrixXd& L) {
	// std::cout.precision(2);

	// Eigenvalues to test that it's positive semidefinite
	EigenSolver<MatrixXd> es;
	es.compute(M, /* computeEigenvectors = */ false);
	bool semiDef = true;
	for (int i=0; i<M.cols(); i++) {
		if ( es.eigenvalues()[i].real() < 0 ) {semiDef = false; }
	}

	if (!semiDef) std::cout << "Matrix isn't semidefinite!!"<<"\n\n";

	if ( M.rows()<21 ) {
		std::cout << "The eigenvalues of M are:\n";
		std::cout << es.eigenvalues().transpose() <<"\n\n";

		std::cout << "M\n\n" << M << "\n\n\n";
		std::cout << "L\n\n" << L << "\n\n\n";
	}

	//Comparison to original
	MatrixXd comp = L * L.transpose() - M;
	if (comp.rows() < 21) {
		std::cout << "Comparison of L*L^dag - M\n\n";
		std::cout << comp << "\n\n";
	}

	// Calculate Norm of Comparison
	double norm = 0;
	for (int r=0; r<L.rows(); r++)
		for (int c=0; c<L.cols(); c++)  {
			norm += comp(r,c)*comp(r,c);
		}
	norm = sqrt(norm);
	std::cout << "Norm of comparison: " <<norm << "\n";
	std::cout << "Number of rows truncated: " << L.cols()-M.cols() << "\n\n";

}


////////////////////////////////////////////////////////////////////////////////
// Reading RDM Functions ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int gen3Idx(int i, int j, int k, int norb ) {
  return i*norb*norb + j*norb + k;
}

int gen4Idx(int i, int j, int k, int l, int norb ) {
  return i*norb*norb*norb+j*norb*norb+k*norb+l;
}

void r2RDM (const char* fIn, MatrixXd& m2) {
  char line[255];

  FILE *fp = fopen( fIn, "r" );
  fgets(line, sizeof(line), fp);
  int norb = atoi( strtok(line, " ,\t\n") );

  while ( fgets(line, sizeof(line), fp) != NULL ) {
    int i = atoi( strtok(line, " ,\t\n") );
    int j = atoi( strtok(NULL, " ,\t\n") );
    int k = atoi( strtok(NULL, " ,\t\n") );
    int l = atoi( strtok(NULL, " ,\t\n") );
    double val = atof( strtok(NULL, " ,\t\n") );
    m2(i*norb+j,k*norb+l) += val;
  }

  return;
}



void r3RDM (const char* fIn, MatrixXd& m3) {
  char line[255];
  FILE *fp = fopen(fIn, "r");
  fgets(line,sizeof(line),fp);
  int norb = atoi( strtok(line," ,\t\n") );

  while (fgets(line, sizeof(line), fp) != NULL) {
    int i = atoi( strtok(line, " ,\t\n") );
    int j = atoi( strtok(NULL, " ,\t\n") );
    int k = atoi( strtok(NULL, " ,\t\n") );
    int l = atoi( strtok(NULL, " ,\t\n") );
    int m = atoi( strtok(NULL, " ,\t\n") );
    int n = atoi( strtok(NULL, " ,\t\n") );
    double val = atof( strtok(NULL, " ,\t\n") );
    m3(gen3Idx(i,j,k,norb),gen3Idx(n,m,l,norb)) += val;
  }
  return;
}

void r4RDM (const char* fIn, MatrixXd& m4) {
  char line[255];
  FILE *fp = fopen(fIn, "r");
  fgets(line,sizeof(line),fp);
  int norb = atoi( strtok(line," ,\t\n") );

  while (fgets(line, sizeof(line), fp) != NULL) {
    int i = atoi( strtok(line, " ,\t\n") );
    int j = atoi( strtok(NULL, " ,\t\n") );
    int k = atoi( strtok(NULL, " ,\t\n") );
    int l = atoi( strtok(NULL, " ,\t\n") );
    int m = atoi( strtok(NULL, " ,\t\n") );
    int n = atoi( strtok(NULL, " ,\t\n") );
    int o = atoi( strtok(NULL, " ,\t\n") );
    int p = atoi( strtok(NULL, " ,\t\n") );
    double val = atof( strtok(NULL, " ,\t\n") );
    m4(gen4Idx(i,j,k,l,norb),gen4Idx(p,o,m,n,norb)) += val;
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Primary Functions ///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void ICCholesky ( MatrixXd& A, MatrixXd& L ) {
	int n = A.rows();
	assert( n == A.cols() );

	for (int i=0; i < n; i++) {
		for (int j=0; j < i + 1; j++) {
			// logg(L); logg("\n");

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
	checkDecomposition(A,L);
	return;
}

void ICCScreened ( MatrixXd& M, std::vector<std::vector<double> >& L,
  double tau ) {
	int n = M.rows();
	int nv = n; //# of Cholesky Vectors
	// Keep in mind that the first index is the column in L!

	// Cholesky-Crout (calculating one column at a time)
	std::vector<double> Lc;
	for (int c=0; c<n; c++) { //Loops over indices of M
		Lc.resize(n-c); std::fill(Lc.begin(), Lc.end(), 0);
		for (int r=c; r<n; r++) {
			// Diagonals
			if (c==r) {
				Lc[r-c] = M(c,c);
				for (int i=0; i<c; i++)
					Lc[r-c] -= L[i][r-i] * L[i][r-i];

				Lc[r-c] = sqrt(Lc[r-c]);
				if ( Lc[r-c] < tau ) { nv-=1; /*std::cout<<"c "<<c<<"\n\n";*/ break; }
			}

			// Off-diagonals
			else {
				Lc[r-c] = M(r,c);
				for (int i=0; i<c; i++)
					Lc[r-c] -= L[i][r-i]*L[i][c-i];

				Lc[r-c] = Lc[r-c]/Lc[c-c];
			}
		}
		L.push_back(Lc);
	}

	checkDecomposition(M,L,nv);
	return;
}


void PCD ( MatrixXd& A, MatrixXd& L, std::vector<size_t>& P, double tau ) {
  /*!
    Adapted Algorithm 3.1 from:	LAPACK-Style Codes for Level 2 and 3 Pivoted
    Cholesky Factorizations by C. Lucas. 2004.

    Here we copy the original matrix the original matrix into the new matrix L
    that will contain our output Cholesky matrix.

    .. note:: Ordering of RDM matrix must be c0 c1 c2 d2 d1 d0 where c0 is paired with d2, c1 is paired with d1, and c2 is paired with d0.

    :Inputs:

      MatrixXd& A:
        Matrix to decompose.
      MatrixXd& L:
         Matrix that stores the decomposition.
      vector<size_t>& P:
        Vector to store the permutation indices for reordering. The permutation
        matrix is P_mat( P[i], i ). TODO merge this with piv.
      double tau:
        Tolerance for truncation of the decomposition.
  */

	size_t n = A.rows();
	std::vector<size_t> piv (n); // Pivot vector
	std::vector<double> dots(n); // Store accumulated dot products
	std::iota(piv.begin(), piv.end(), 0); // Fill pivot vector
	std::fill(dots.begin(), dots.end(), 0); // Fill dots vector W/ ZEROS

	for (int c=0; c<n; c++) {
		for (int r=0; r<n; r++) { if (true) { L(r,c) = A(r,c); } } // TODO
	}


	size_t q; // Used to store the index of the maximum diagonal
	double qmax;
	size_t lrank = n;

	// Loop over cholesky vectors
	for (int j=0; j<n; j++) {
		// Update dots
		if (j>0) {
			for (int i=j; i<n; i++) {
				dots[i] += L(i,j-1)*L(i,j-1);
			}
		}

		// Find qmax and q
		qmax = 0;
		for (int p=j; p<n; p++) {
			if ( L(p,p) - dots[p] > qmax ) { q = p; qmax = L(p,p)-dots[p]; }
		}

		// Stopping criterion
		if (qmax < tau) {
			lrank = j-1;
			//std::cout<<"Exiting decomposition early...   ";
			//std::cout<<"Qmax = " <<qmax << "   ";
			//std::cout<< "Rank of L: " << lrank << "   ";
			break;
		}

		// Swapping rows and columns
		if ( q != j ){
			L.row(j).swap(L.row(q));
			L.col(j).swap(L.col(q));

			// Swap dots and pivot
			std::swap(dots[j],dots[q]);
			std::swap(piv[j],piv[q]);
		}

		// Diagonals
		L(j,j) -= dots[j];
		L(j,j) = sqrt(L(j,j));

		// Update the jth column
		for (int k=j+1; k<n; k++) {
			if (j>0 && j<n) {
				for (int l=0; l<j; l++) {
					L(k,j) -= L(k,l)*L(j,l);
				}
			}

			if (j<n) {
				L(k,j) = L(k,j)/L(j,j);
			}
		}
	} // end main loop over Cholesky vectors

	// Build Permutation matrix
	MatrixXd Pm = MatrixXd::Zero(n,n);
	for (int i=0; i<piv.size(); i++) { Pm(piv[i],i) = 1;  }

	// Setting values of L above diagonal to zero (TODO)
	for (int c=0; c<n; c++) {
		for (int r=0; r<n; r++) {
			if (c>r) { L(r,c) = 0; }
			if (c>lrank) { L(r,c) = 0; }
		}
	}

	// Printing Results
	MatrixXd comp = L*L.transpose()-Pm.transpose()*A*Pm;
	std::cout << tau << "\t";
	std::cout <<  comp.squaredNorm() << "\t";
	std::cout << A.rows() << "\t";
	std::cout << lrank << "\n";
	return;
}


void Test2RDM ( MatrixXd& m2, int& nelec, int& norb ) {
  MatrixXd m1 = MatrixXd::Zero(norb,norb);

  for (int i=0; i < norb; i++ ) {
    for (int j=0; j < norb; j++ ) {
      for (int k=0; k < norb; k++ ) {
	m1(i,j) += m2(i*norb+k,j*norb+k)/(nelec-1);
      }
    }
  }

  double chkelec = 0;
  for (int i=0; i < norb; i++ ) {
    chkelec += m1(i,i);
  }

  std::cout << "Check Electrons " << chkelec << std::endl; //TODO
}

void Test3RDM ( MatrixXd& m3, int& nelec, int& norb ) {
  MatrixXd m2 = MatrixXd::Zero(norb*norb,norb*norb);

  for (int i=0; i<norb; i++) {
    for (int j=0; j<norb; j++) {
      for (int k=0; k<norb; k++) {
	for (int l=0; l<norb; l++) {
	  for (int m=0; m<norb; m++) {
	    m2(i*norb+j,l*norb+k)+= m3(gen3Idx(i,j,m,norb),gen3Idx(m,k,l,norb));
	    m2(i*norb+j,l*norb+k) /= (nelec-2);
	  }
	}
      }
    }
  }

  MatrixXd m1;
  Test2RDM(m2,nelec,norb);

}


/*******************************************************************************
******************************** Main ******************************************
*******************************************************************************/
int main(int argc, char** argv) {


	if (true) {
	  if (argc == 1) { return -1; }

	  int n = atoi(argv[1]);
	  int nelec = atoi(argv[2]);
	  const char *fIn = "./spatialRDM.0.0.txt";
	  //	  int n = 22;
	  //	  const char *fIn = "../testing/cholesky/5cene_2rdm.txt";

	  //int n=29;
	  //const char *fIn = "fep_2rdm.txt";

	  // 2RDM
	  std::cout <<"2RDM Test\n===========================\n\n";
	  std::cout << "Full number of rows = " << n*n << "\n\n";
	  std::cout << "Error Tol.\tNorm\tA Rank\tL Rank\n";
	  MatrixXd m = MatrixXd::Zero(n*n,n*n);
	  r2RDM(fIn,m);

	  MatrixXd L = MatrixXd::Zero(n*n,n*n);
	  std::vector<size_t> P (0);
	  for (int taue=1; taue<11; taue++) {
	    //std::cout<< "Tau: 1e-"<<taue<<"\n";
	    PCD(m, L, P, pow(10,-taue) );
	  }

	  // 3RDM
	  std::cout <<"3RDM Test\n===========================\n\n";
	  std::cout << "Error Tol.\tNorm\tA Rank\tL Rank\n";
	  MatrixXd m3 = MatrixXd::Zero(n*n*n,n*n*n);
	  const char *fIn3 = "./spatial3RDM.0.0.txt";
	  r3RDM(fIn3,m3);
	  Test3RDM(m3,nelec,n);

	  MatrixXd L3 = MatrixXd::Zero(n*n*n,n*n*n);
	  std::cout << L3.rows() << "\n\n";
	  //std::cout << m3(7*64+7*8+6,7*64+7*8+6) << "\n\n";
	  //std::cout <<

	  std::vector<size_t> P3;
	  for (int taue=1; taue<11; taue++) {
	    //std::cout<< "Tau: 1e-"<<taue<<"\n";
	    PCD(m3, L3, P3, pow(10,-taue) );
	  }


	  // 4RDM
	  /*std::cout <<"4RDM Test\n===========================\n\n";
	  MatrixXd m4 = MatrixXd::Zero(n*n*n*n,n*n*n*n);
	  r4RDM("../testing/cholesky/spatial4RDM.0.0.txt",m4);
	  MatrixXd L4 = MatrixXd::Zero(n*n*n*n,n*n*n*n);
	  std::cout << L4.rows() << "\n\n";
	  std::cout << m4(1,7) << "\n\n";
	  //std::cout <<

	  std::vector<size_t> P4;
	  for (int taue=1; taue<10; taue++) {
	    std::cout<< "Tau: 1e-"<<taue<<"\n";
	    PCD(m4, L4, P4, pow(10,-taue) );
	  }*/


	}



	if (false) {
		int n = 400;
		// Initialize Matrix
		MatrixXd B = MatrixXd::Random(n,n); B = B*B; B = B * B.transpose();

		//Calculate Eigenvectors/values
		EigenSolver<MatrixXd> es(B);

		//Change Eigenvalues to turn on degeneracy
		MatrixXd lam = MatrixXd::Zero(n,n);
		for (int r=0; r<n; r++)
			for (int c=0; c<n; c++) {
				if (r==c) {lam(c,c)=es.eigenvalues()[c].real(); }
			}

		// Add degeneracy
		for (int d=0; d<n/2; d++) {
		  lam(2*d+1,2*d+1) = 0;//2 * es.eigenvalues()[2*d].real();
			//es.eigenvectors().col(2*d+1) = 2*es.eigenvectors().col(2*d);
		}

		// Recombine
		B = es.eigenvectors().real()*lam*es.eigenvectors().transpose().real();

		// Run Tests
		std::vector<size_t> piv(n);
		MatrixXd L = MatrixXd::Zero(n,n);

		for (int taue=2; taue<16; taue++) {
			std::cout << "Tau = " << pow(10,-taue) << "    ";
			if (n<21) { std::cout << "\n";}
			PCD(B,L,piv,pow(10,-taue));
		}
	}

	if (false) {
		for (int n=100; n<150; n++) {
			MatrixXd B = MatrixXd::Random(n,n); B = B*B; B = B * B.transpose();
			// std::cout << "B\n" << B << "\n\n";
			std::cout << "n = " << n << "   ";//<< "\n";

			std::vector<size_t> piv(n);
			MatrixXd L = MatrixXd::Zero(n,n);
			PCD(B,L,piv,1e-16);
		}
	}

	if (false) {

		int n = 6;
		MatrixXd B = MatrixXd::Random(n,n); B = B*B; B = B * B.transpose();
		// std::cout << "B\n" << B << "\n\n";
		std::cout << "n = " << n << "\n\n";//<< "\n";

		std::vector<size_t> piv(n);
		MatrixXd L = MatrixXd::Zero(n,n);
		PCD(B,L,piv,1e-16);

		// Eigen Tests
		std::cout << "==================================================\n\n";
		LDLT<MatrixXd> ldltOfB(B); // compute the Cholesky decomposition of A
		MatrixXd res (n,n); res.setIdentity();
		MatrixXd Leigen = ldltOfB.matrixL();
		MatrixXd diag = ldltOfB.vectorD().real().asDiagonal();
		// std::cout << diag << "\n\n";
		for (int d=0; d<n; d++) {diag(d,d) = sqrt(diag(d,d)); }
		// std::cout << diag << "\n\n";
		MatrixXd test = Leigen*diag;

		res.setIdentity();
		res = ldltOfB.transpositionsP() * res;
		// L^* P
		res = ldltOfB.matrixU() * res;
		// D(L^*P)
		res = ldltOfB.vectorD().real().asDiagonal() * res;
		// L(DL^*P)
		res = ldltOfB.matrixL() * res;
		// P^T (LDL^*P)
		res = ldltOfB.transpositionsP().transpose() * res;


		std::cout << "The Cholesky factor L is" << std::endl << Leigen << "\n\n";
		std::cout << "To check this, let us compute L * L.transpose()" << "\n\n";
		std::cout << res - B << "\n\n";
		std::cout << "TEST\n";
		std::cout << test << "\n\n";
	}


	return 0;
}
