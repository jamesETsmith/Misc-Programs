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
// #include <stdlib.h>
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

// Logging Functions //
void logg ( char const * msg ) {
	std::cout << msg << std::endl;
	return;
}

void logg ( int msg ) {
	std::cout << msg << std::endl;
	return;
}

void logg ( double d ) {
	std::cout << d << std::endl;
}

void logg ( MatrixXd& m ) {
	std::cout << m << std::endl;
}

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

void makeLDSMatrix (MatrixXd& M) {
	size_t n = M.rows();
	std::cout << M << "\n\n";

	for (int i=0; i<n; i++) {
		M(i,n-2) = 1;
		M(i,n-1) = 1;
		M(n-2,i) = 1;
		M(n-1,i) = 1;
	}
	std::cout<< M << "\n\n";
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

void checkDecomposition(MatrixXd& M, MatrixXd& L, std::vector<double>& D,
  std::vector<int>& Ll, std::vector<size_t>& idx, int n, int nv) {

	if(n<21) {
		EigenSolver<MatrixXd> es;
		es.compute(M, /* computeEigenvectors = */ false);
		std::cout << "The eigenvalues of M are:\n";
		std::cout << es.eigenvalues().transpose() <<"\n\n";

		// Print out ordering
		std::cout << "Index Dictionary\n";
		for (int i=0; i<n; i++) std::cout<<idx[i]<<" ";
		std::cout << "\n\n";

		// Print D
		std::cout << "Diagonals\n";
		for (int i=0; i<n; i++) std::cout<<D[idx[Ll[i]]]<<" ";
		std::cout << "\n\n";
	}

	MatrixXd Mp (n,n); Mp.setZero(n,n);
	for (int r=0; r<L.rows(); r++)
		for (int c=0; c<L.cols(); c++) Mp(idx[r],idx[c]) = L(r,c);    //Mp(r,c)=L(r,c); //

	// reorderBasis(Mp,idx);

	// Debugging //TODO
	// Standard Output
	if (Mp.rows() < 21) {
		std::cout << "M\n\n" << M << "\n\n\n";
		std::cout << "L\n\n" << L << "\n\n\n";
	}
	std::cout.precision(2);

	//Comparison to original
	MatrixXd comp = Mp * Mp.transpose() - M;
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
	std::cout << "Number of rows truncated: " << n-nv << "\n\n";
}


////////////////////////////////////////////////////////////////////////////////
// Reading RDM Functions ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int gen3Idx(int i, int j, int k, int norb ) {
  return i*norb*norb+j*norb+k;
}

int gen4Idx(int i, int j, int k, int l, int norb ) {
  return i*norb*norb*norb+j*norb*norb+k*norb+l;
}

void r2RDM (char* fIn, MatrixXd& m2) {
	char line[255];

	//char const *fIn = "o2_2rdm.txt";
	//size_t chkorbs = 8;
	// printf("%zi\n",chkorbs);

	FILE *fp = fopen( fIn, "r" );
	fgets(line, sizeof(line), fp);
	int norb = atoi( strtok(line, " ,\t\n") );
	//assert (norb==chkorbs);
	// printf("%i\n", norb);

	//MatrixXd m2(norb*norb,norb*norb);

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



void r3RDM (char* fIn, MatrixXd& m3) {
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
    m3(gen3Idx(i,j,k,norb),gen3Idx(l,m,n,norb)) += val;
  }
  return;
}

void r4RDM (char* fIn, MatrixXd& m4) {
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
    m4(gen4Idx(i,j,k,l,norb),gen4Idx(m,n,o,p,norb)) += val;
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


void OCC ( MatrixXd& M, double d, double s, double tau,
  std::vector< std::vector<double> >& L, std::vector<size_t>& idx  ) {
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
	std::vector<double> D(n);
	for (int p=0; p < n; p++) {
		D[p] = M(p,p);
	}
	double D_max = *std::max_element(D.begin(),D.end()); // TODO
	idx = sort_indexes(D);
	// std::iota(idx.begin(), idx.end(), 0);
	// for (int i=0; i < idx.size(); i++)
	// std::cout << D[idx[i]] << std::endl;

	// Step 2
	std::vector<int> Ll (0);
	for (int p=0; p<n; p++) {
		if ( d * sqrt(D_max*D[p]) > tau ) {
			Ll.push_back( p );
		}
	}

	// Step 3
	int nv = 0;

	// Step 4
	int i = 0;

	// Step 5
	while ( D_max > tau && nv < n ) {
		//// a
		i++;
		//// b
		double D_min = std::max(s*D_max,tau);
		//// c
		std::vector<int> Q (0);
		for (int q=0; q<Ll.size(); q++) {
			if ( D[Ll[q]] > D_min ) {
				Q.push_back( q );
			}
		}

		//// d
		MatrixXd M_pq (Ll.size(),Q.size());
		for (int p=0; p<Ll.size(); p++)
			for (int q=0; q<Q.size(); q++) {
				M_pq(p,q) = M( Ll[p], Q[q] );
			}

		//// e
		MatrixXd delta_pq (Ll.size(),Q.size());
		for (int p=0; p<Ll.size(); p++)
			for (int q=0; q<Q.size(); q++) {
				delta_pq(p,q) = M_pq(p,q);
				for (int J=0; J<nv; J++) {
					delta_pq(p,q) -= L[J][ idx[Ll[p]] - J ]*L[J][ idx[Q[q]] - J ];
				}
			}

		//// f
		double Q_max = 0;
		int q_j = 0;
		for (int q=0; q<Q.size(); q++) {
			if ( D[Q[q]] > Q_max ) {
				q_j = q;
				Q_max = D[Q[q]];
			}
		}
		//// g
		int j = 0;
		int J;
		std::vector<double> Lj;
		int LjSize;
		//// h
		while ( j < Q.size() && Q_max > D_min) {
			////// i
			J = nv + j;
			LjSize = n-J;
			// std::cout<< "\n" << "LjSize " <<LjSize << "\n";//TODO
			Lj.resize(LjSize);
			L.push_back(Lj);

			////// ii
			// logg(L); logg("\n");

			////// iii
			for (int p=0; p<Ll.size(); p++) {
				if (idx[Ll[p]] >= J ) {
					// std::cout << "Delta(p,q) " << delta_pq(p,q_j)/sqrt(Q_max) << "\n"; //TODO
					// std::cout << "Idx " <<idx[Ll[p]]<<"\n";
					L[J][ idx[Ll[p]] - J ] = delta_pq(p,q_j)/sqrt(Q_max);
				}
			}

			////// iv
			for (int p=0; p<Ll.size(); p++) {
				D[p] -= (L[J][idx[Ll[p]]-J] * L[J][idx[Ll[p]]-J]);
				for (int q=0; q<Q.size(); q++) {
					if ( idx[Ll[p]] >= J && idx[Q[q]] >= J) {
						delta_pq(p,q) -= L[J][ idx[Ll[p]] - J ]*L[J][ idx[Q[q]] - J ];
					}
				}
			}
			Q_max = 0;
			for (int q=0; q<Q.size(); q++) {
				if ( D[Q[q]] > Q_max ) {
					q_j = q;
					Q_max = D[Q[q]];
				}
			}

			// std::cout << "Lj\n\n";
			// for (int i=0; i<L[J].size(); i++)
			//  std::cout<<L[J][i]<<std::endl;
			// std::cout << "\n";
			// L.push_back(Lj);
			//TODO Remove
			// std::vector<double> tst = {0,1,1};
			// L.push_back(tst);
			// Lj.erase(Lj.begin(),Lj.end());
			// L[J] = Lj;
			// std::cout << "j " << j << "\n";
			// for (int a=0; a<n-j; a++) {
			//  std::cout<<L[J][a]<<" ";
			//  std::cout << "\n";
			// for (int b=0; b<n; b++) {
			//  if ( a < L[b].size() && b<= J ) std::cout << L[b][a] << " ";
			//  else std::cout << 0 << " ";
			// }
			// std::cout << "\n";
			// }

			j++;

		}

		//// i
		nv += j;
		D_max = 0;
		for (int p=0; p<Ll.size(); p++) {
			if ( D[Ll[p]] > D_max ) {
				D_max = D[Ll[p]];
			}
		}

		//// k
		Ll.clear();
		for (int p=0; p<n; p++) {
			if ( d*sqrt(D_max*D[Ll[p]]) > tau ) {
				Ll.push_back(p);
			}
		}
	}
	// std::cout.precision(6);
	// std::cout << std::fixed;
	// std::cout << "Cholesky Decomposition of M" << std::endl << L << std::endl;
	// std::cout << std::endl;
	return;
}


void OCCholesky ( MatrixXd& M, double d, double s, double tau, MatrixXd& L,
  std::vector<size_t>& idx  ) {
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
	std::vector<double> D(n);
	for (int p=0; p < n; p++) {
		D[p] = M(p,p);
	}
	double D_max = *std::max_element(D.begin(),D.end()); // TODO
	idx = sort_indexes(D);
	// std::iota(idx.begin(), idx.end(), 0);
	// for (int i=0; i < idx.size(); i++)
	// std::cout << D[idx[i]] << std::endl;

	// Step 2
	std::vector<int> Ll (0);
	for (int p=0; p<n; p++) { //iterate over sorted indices
		if ( d * sqrt(D_max*D[idx[p]]) > tau ) {
			Ll.push_back( p ); //store sorted indices
		}
	}

	// Step 3
	int nv = 0;

	// Step 4
	int i = 0;

	// Step 5
	while ( D_max > tau && nv < n ) {
		//// a
		i++;
		//// b
		double D_min = std::max(s*D_max,tau);
		//// c
		std::vector<int> Q (0);
		for (int q=0; q<Ll.size(); q++) { // iterate over sorted indices
			if ( D[idx[Ll[q]]] > D_min ) { //convert back to unsorted to access
				Q.push_back( q ); // store sorted indices
			}
		}

		//// d
		MatrixXd M_pq (Ll.size(),Q.size());
		for (int p=0; p<Ll.size(); p++)
			for (int q=0; q<Q.size(); q++) { // iterate over sorted
				M_pq(Ll[p],Q[q]) = M( idx[Ll[p]], idx[Q[q]] ); // conver back to unsorted
			}

		//// e
		MatrixXd delta_pq (Ll.size(),Q.size());
		delta_pq.setZero(Ll.size(),Q.size());
		for (int p=0; p<Ll.size(); p++)
			for (int q=0; q<Q.size(); q++) { // sorted indices
				// if (Ll[p] >= Q[q] ) {
				delta_pq(Ll[p],Q[q]) = M_pq(Ll[p],Q[q]);
				for (int J=0; J<nv; J++) {
					delta_pq(Ll[p],Q[q]) -= L(Ll[p],J)*L(Q[q],J);   //sorted indices
				}
				// }
			}

		//// f
		double Q_max = 0;
		int q_j = 0;
		for (int q=0; q<Q.size(); q++) { //sorted indices
			if ( D[idx[Q[q]]] > Q_max ) {
				q_j = q;
				Q_max = D[idx[Q[q]]];
			}
		}
		//// g
		int j = 0;
		//// h
		while ( j < Q.size() && Q_max > D_min) {
			////// i
			int J = nv + j;

			////// iii
			std::cout << delta_pq << "\n\n"; //TODO
			for (int p=0; p<Ll.size(); p++) { //sorted
				L(Ll[p],J) = delta_pq(Ll[p],q_j)/sqrt(Q_max);
			}

			////// iv
			for (int p=0; p<Ll.size(); p++) { //sorted
				D[idx[p]] -= (L(Ll[p],J) * L(Ll[p],J)); //convert back to unsorted
				for (int q=0; q<Q.size(); q++) { //sorted
					// if (Ll[p] >= Q[q] ) delta_pq(Ll[p],Q[q]) -= L(Ll[p],J)*L(Q[q],J);//TODO
					delta_pq(Ll[p],Q[q]) -= L(Ll[p],J)*L(Q[q],J);
				}
			}
			Q_max = 0;
			for (int q=0; q<Q.size(); q++) {
				if ( D[idx[Q[q]]] > Q_max ) {
					q_j = q;
					Q_max = D[idx[Q[q]]];
				}
			}
			j++;
		}

		//// i
		nv += j;
		D_max = 0;
		for (int p=0; p<Ll.size(); p++) { // sorted
			if ( D[idx[Ll[p]]] > D_max ) {
				D_max = D[idx[Ll[p]]];
			}
		}

		//// k
		Ll.clear();
		for (int p=0; p<n; p++) { //sorted
			if ( d*sqrt(D_max*D[idx[Ll[p]]]) > tau ) {
				Ll.push_back(p);
			}
		}
	}
	checkDecomposition(M,L,D,Ll,idx,n,nv);
	return;
}

void PCD ( MatrixXd& A, MatrixXd& L, std::vector<size_t>& P, double tau ) {
	/*!
	Adapted Algorithm 3.1 from:	LAPACK-Style Codes for Level 2 and 3 Pivoted
	Cholesky Factorizations by C. Lucas. 2004.

	Here we copy the original matrix the original matrix into the new matrix L
	that will contain our output Cholesky matrix.

	Inputs:
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


	// MatrixXd Acopy = A; // Keep a copy for analysis
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
			std::cout<<"Exiting decomposition early...   ";
			std::cout<< "Rank of L: " << lrank << "   ";
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

	MatrixXd Pm = MatrixXd::Zero(n,n);
	for (int i=0; i<piv.size(); i++) { Pm(piv[i],i) = 1;  }
	for (int c=0; c<n; c++) {
		for (int r=0; r<n; r++) {
			if (c>r) { L(r,c) = 0; }
			if (c>lrank) { L(r,c) = 0; }
		}
	}

	// Printing Results
	// std::cout << "==================================================\n\n";
	MatrixXd comp = L*L.transpose()-Pm.transpose()*A*Pm;

	if (n<21) {
		std::cout << L << "\n\n";
		std::cout << comp << "\n\n";
	}

	std::cout << "Norm of comparison: " << comp.squaredNorm() << "\n\n";
	return;
}


void Calculate1RDM ( MatrixXd& m2, MatrixXd& m1, int& nelec ) {
	int norb = m1.rows();

	for (int i=0; i < norb; i++ ) {
		for (int j=0; j < norb; j++ ) {
			// m1(i,j) = 0;
			// std::cout << i << " " << j << std::endl;
			for (int k=0; k < norb; k++ ) {
				m1(i,j) += m2(i*norb+k,j*norb+k)/(nelec-1);
			}
		}
	}

	double chkelec = 0;
	for (int i=0; i < norb; i++ ) {
		chkelec += m1(i,i);
	}

	std::cout << chkelec << std::endl; //TODO
}

void testForRandomMatrix ( int rn ) {
	MatrixXd r = MatrixXd::Random(rn,rn);
	r = r * r.transpose();

	std::vector< std::vector<double> > rdv( 0, std::vector<double> (0) );
	std::vector<size_t> ridx(rn);
	MatrixXd rd (rn,rn); rd.setZero(rn,rn);
	OCC(r,1,0,1e-16,rdv,ridx);

	std::cout << "\n\n";

	for (int col=0; col < rdv.size(); col++) {
		std::cout << "Size of Cholesky vector " << rdv[col].size() << std::endl;
		for (int row=col; row < rdv[0].size(); row++) {
			rd(row,col) = rdv[col][row-col];
		}
	}

	return;
}

/*******************************************************************************
******************************** Main ******************************************
*******************************************************************************/
int main() {
	// Read in 2RDM
	/*
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
	*/

  size_t norb = 8;
  MatrixXd m2 = MatrixXd::Zero(norb*norb,norb*norb);

	// Create 1RDM
	MatrixXd m1(norb,norb);
	m1 = MatrixXd::Zero(norb,norb);

	int nelec = 12;
	//Calculate1RDM(m2,m1,nelec);

	// Test several random matrices
	// for (int i=20; i < 21; i++) {
	//  std::cout << "Size of random matrix " << i << "\n\n";
	//  testForRandomMatrix(i);
	// }

	if (true) {
	  int n = 8;
	  char const *fIn = "../testing/cholesky/spatialRDM.0.0.txt";
	  MatrixXd m = MatrixXd::Zero(n*n,n*n);
	  r2RDM(fIn,m);

	  std::cout << m(0,0);
	  /*
	  MatrixXd L = MatrixXd::Zero(n*n,n*n);
	  std::vector<size_t> P (0);
	  for (int taue=1; taue<11; taue++) {
	    std::cout<< "Tau: 1e-"<<taue<<"\n";
	    PCD(m2, t, tidx, pow(10,-taue) );
	  }
	  */
	}

	if (false) {
		std::cout <<"2RDM Test\n============\n\n";
		MatrixXd t (norb*norb,norb*norb);
		t.setZero(norb*norb,norb*norb);

		std::vector<size_t> tidx (0);
		for (int taue=1; taue<11; taue++) {
			std::cout<< "Tau: 1e-"<<taue<<"\n";
			PCD(m2, t, tidx, pow(10,-taue) );
		}
	}


	if (false) {
		int n = 100;
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
			lam(2*d+1,2*d+1) = 2 * es.eigenvalues()[2*d].real();
			es.eigenvectors().col(2*d+1) = 2*es.eigenvectors().col(2*d);
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
		for (int n=80; n<100; n++) {
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

	// 4X4
	if (false) {
		int n = 6;
		MatrixXd B = MatrixXd::Random(n,n); B = B*B; B = B * B.transpose();

		// MatrixXd B(4,4);
		// B << 10, 4, 4, -4, 4, 16, 4, 2, 4, 4, 6, -2, -4, 2, -2, 4;
		// std::cout << B << "\n\n";

		std::vector<size_t> bidx (0);
		MatrixXd L (n,n); L.setZero(n,n);
		OCCholesky(B,1,0,1e-10,L,bidx);
		// OCCScreened(B,L,bidx,1e-10);
	}

	if (false) {
		int n = 100;
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
		// std::cout << lam << "\n\n"; // Check that LD was introduces correctly
		lam(n-2,n-2)=es.eigenvalues()[n-3].real();
		lam(n-1,n-1)=es.eigenvalues()[n-4].real();
		// lam(n-5,n-5)=es.eigenvalues()[n-4].real();
		// std::cout << lam << "\n\n"; // Check that LD was introduces correctly
		B = es.eigenvectors().real()*lam*es.eigenvectors().transpose().real();

		//Decompose
		// MatrixXd L (n,n); L.setZero(n,n);
		std::vector< std::vector<double> > L;
		for (int taue=15; taue<16; taue++) {
			std::cout<<"================================"<<"\n";
			std::cout<< "Tau: 1e-"<<taue<<"\n";
			ICCScreened(B,L,pow(10,-taue));
		}
		MatrixXd LL(n,n); LL.setZero(n,n); //TODO
		ICCholesky(B,LL);
		// std::cout << LL-L << "\n\n";

	}


	// 2RDM
	if (false) {
		std::cout <<"2RDM Test\n============\n\n";
		MatrixXd t (norb*norb,norb*norb);
		t.setZero(norb*norb,norb*norb);

		std::vector<size_t> tidx (0);
		for (int taue=1; taue<11; taue++) {
			std::cout<< "Tau: 1e-"<<taue<<"\n";
			OCCholesky(m2, 1,0, pow(10,-taue), t, tidx);
		}
	}

	return 0;
}
