/*
   Test file for plain cholesky decomposition of matrix.

   Compile with:
   g++ -std=c++11 -g cholesky.cpp -I/usr/local/include/eigen3/ -o cholesky
   or
   g++ -std=c++11 -g cholesky.cpp -I/home/james/Documents/Apps/eigen/ -o cholesky

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

// Primary Functions //
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
	// std::cout << "Cholesky Decomposition of A" << std::endl << L << std::endl;
	// std::cout << std::endl;
	return;
}

void OCCScreened (MatrixXd& M, MatrixXd& L, std::vector<size_t>& bidx,
	double tau ) {

		// Calculate and sort diagonals
		int n = M.rows();
		std::vector<double> D(n);
		for (int p=0; p < n; p++) {
			D[p] = M(p,p);
		}
		bidx = sort_indexes(D); // Get sorted ordering
		std::sort(D.begin(),D.end(),descend); // Sort vector // TODO combine these

		// Screen the number of columns we are going to look at
		int maxJ = 0;
		for (int c=0; c<D.size(); c++) {
			std::cout << sqrt(D[0]*D[c])<< " tau: " <<tau<<"\n\n";
			if (sqrt(D[0]*D[c]) > tau) maxJ ++;
		}

		// Generate screened matrix that we'll decompose with standard cholesky.
		// TODO this should get deleted in the true out of core version
		MatrixXd Mpq(maxJ,maxJ);
		for (int r=0; r<maxJ; r++)
			for (int c=0; c<maxJ; c++)
				Mpq(r,c) = M(bidx[r],bidx[c]);

		// Initialize matrix to store decomposition
		// MatrixXd L(maxJ,maxJ); L.setZero(maxJ,maxJ);

		// Call Standard cholesky
		ICCholesky(Mpq,L);

		/////////////////////////
		// Check decomposition //
		/////////////////////////
		MatrixXd Mp (n,n); Mp.setZero(n,n);
		for (int r=0; r<maxJ; r++)
			for (int c=0; c<maxJ; c++) Mp(bidx[r],bidx[c]) = L(r,c);

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
		for (int r=0; r<maxJ; r++)
			for (int c=0; c<maxJ; c++)  {
				 norm += comp(r,c)*comp(r,c);
			 }
		norm = sqrt(norm);
		std::cout << "Norm of comparison: " <<norm << "\n\n";
		std::cout << "Number of rows truncated: " << n-maxJ << "\n\n";

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
				if (idx[Ll[p]] >= J )	{
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
			// 	std::cout<<L[J][i]<<std::endl;
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
					delta_pq(p,q) -= L(idx[Ll[p]],J)*L(idx[Q[q]],J);
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
		//// h
		while ( j < Q.size() && Q_max > D_min) {
			////// i
			int J = nv + j;

			////// ii
			// logg(L); logg("\n");

			////// iii
			for (int p=0; p<Ll.size(); p++) {
				L(idx[Ll[p]],J) = delta_pq(p,q_j)/sqrt(Q_max); // TODO Not sure if q_j is correct here
			}

			////// iv
			for (int p=0; p<Ll.size(); p++) {
				D[p] -= (L(idx[Ll[p]],J) * L(idx[Ll[p]],J));
				for (int q=0; q<Q.size(); q++) {
					delta_pq(p,q) -= L(idx[Ll[p]],J)*L(idx[Q[q]],J);
				}
			}
			Q_max = 0;
			for (int q=0; q<Q.size(); q++) {
				if ( D[Q[q]] > Q_max ) {
					q_j = q;
					Q_max = D[Q[q]];
				}
			}
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
	// std::cout << "R\n" << r << std::endl;

	// Eigenvalues
	// EigenSolver<MatrixXd> es;
	// es.compute(r, /* computeEigenvectors = */ false);
	// std::cout << "The eigenvalues of A are:\n";
	// std::cout << es.eigenvalues().transpose() <<"\n\n";

	// Decomposition
	// rd = MatrixXd::Zero(rn,rn);
	std::vector< std::vector<double> > rdv( 0, std::vector<double> (0) );
	std::vector<size_t> ridx(rn);
	MatrixXd rd (rn,rn); rd.setZero(rn,rn);
	OCC(r,1,0,1e-16,rdv,ridx);

	// std::cout << "Indices to Reorder\n";
	// for (int i=0; i<rn; i++)
	//  std::cout<< ridx[i] << " , " << std::endl;
	std::cout << "\n\n";

	// Convert vector<vector> to MatrixXd

	for (int col=0; col < rdv.size(); col++) {
		std::cout << "Size of Cholesky vector " << rdv[col].size() << std::endl;
		for (int row=col; row < rdv[0].size(); row++) {
			rd(row,col) = rdv[col][row-col];
		}
	}
	std::cout << "Decomposed Matrix\n\n" << rd <<"\n";
	reorderBasis(rd,ridx);


	MatrixXd rnew = rd*rd.transpose() - r;
	std::cout << "R_{new}\n" << rnew << std::endl;

	double rnorm =0;
	for (int i=0; i<rn; i++)
		for (int j=0; j<rn; j++) {
			rnorm += rnew(i,j) * rnew(i,j);
		}

	std::cout << "Norm of random matrix " << sqrt(rnorm) << std::endl;

	MatrixXd rdnew (rn,rn); rdnew.setZero(rn,rn);
	OCCholesky(r,1,0,1e-16,rdnew,ridx);
	std::cout << rdnew << std::endl;
	reorderBasis(rdnew,ridx);
	std::cout<< "OCC Comparison\n\n" << rdnew << std::endl;
	return;
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

	int nelec = 12;
	Calculate1RDM(m2,m1,nelec);

	// Test several random matrices
	// for (int i=20; i < 21; i++) {
	//  std::cout << "Size of random matrix " << i << "\n\n";
	//  testForRandomMatrix(i);
	// }

	// 4X4
	// MatrixXd B(4,4);
	// B << 10, 4, 4, -4, 4, 16, 4, 2, 4, 4, 6, -2, -4, 2, -2, 4;
	if (false) {
		int n = 20;
		MatrixXd B = MatrixXd::Random(n,n); B = B*B; B = B * B.transpose();
		// std::vector<size_t> Bidx = {3,2,0,1};
		std::vector<size_t> bidx (0);
		MatrixXd Bcopy = B.replicate(B.rows(),B.cols());
		MatrixXd L (n,n); L.setZero(n,n);
		OCCScreened(B,L,bidx,1e-10);
	}


	// 2RDM
	MatrixXd t (norb*norb,norb*norb);
	t.setZero(norb*norb,norb*norb);

	std::vector<size_t> tidx (0);
	OCCScreened(m2, t, tidx, 0.5);
	// reorderBasis(t,tidx);


	return 0;
}
