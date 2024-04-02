

#ifndef _FTS_MOLEC_PARTICLE
#define _FTS_MOLEC_PARTICLE
#include "fts_molecule.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/complex.h>

class FTS_Box;

class ParticleMolec : public FTS_Molec {
	protected:

	public:
		int particleNum;
		float phiNP;
		float Vnptot;
		std::string particleSpecies;
		thrust::device_vector<double> d_Rp; // radius of each particle
		thrust::host_vector<double> Rp; // radius of each particle
		
		thrust::device_vector<double> d_xi; //xi of each particle
		thrust::host_vector<double> xi; //xi of each particle 

		thrust::device_vector<double> d_center; //center x,y,z position of each particle
		thrust::host_vector<double> center; //center xyz position of each particle
		
		thrust::device_vector<thrust::complex<double>> d_NPdensity;
		
		thrust::host_vector<int> intSpecies;
		thrust::device_vector<int> d_intSpecies;
		~ParticleMolec();
		ParticleMolec(std::istringstream& iss, FTS_Box*);

		void calcDensity() override;
		void computeLinearTerms() override;
		void calcHamiltonian();
};

#endif
