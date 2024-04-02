#ifndef _FTS_POTEN_PARTICLE
#define _FTS_POTEN_PARTICLE 


#include "fts_potential.h"


class FTS_Box;

class PotentialParticle : public FTS_Potential {

	private:

	public:
		PotentialParticle();
		PotentialParticle(std::istringstream &iss, FTS_Box*);
		~PotentialParticle();
		
		void updateFields() override;
		std::complex<double> calcHamiltonian() override;
		void writeFields(int) override;
		void initLinearCoeffs() override;

		double chiN; //Strength of the potential
		std::string typeI, typeJ; //Types involved in the potential
		int intTypeI, intTypeJ; //species integer for types I, J 
		thrust::host_vector<int> intSpecies;
};

#endif
