#include "include_libs.h"
#include "fts_potential_particle.h"
#include "FTS_Box.h"
#include "fts_species.h"




PotentialParticle::PotentialParticle(std::istringstream& iss, FTS_Box* p_box) : FTS_Potential(iss, p_box) {
	iss.seekg(0);
	std::string s1;
	iss >> s1;
	iss >> s1;

	potentialStyle = "Particle";
	
	iss >> typeI;
	iss >> typeJ;
	std::cout << "Particle Flory-like potential acting on species " << typeI << ", " << typeJ << std::endl;

	actsOn.push_back(typeI);
	actsOn.push_back(typeJ);

	iss >> chiN;
	
	//Assign species I and J
	intSpecies.resize(2); 
	for (int i=0 ; i<mybox->Species.size() ; i++ ) {
		if ( actsOn[0] == mybox->Species[i].fts_species ) {
			intSpecies[0] = i;
		}
		else if (actsOn[1] == mybox->Species[i].fts_species ) {
			intSpecies[1] = i;
		}
	}
	wpl.resize(mybox->M);
	d_wpl.resize(mybox->M);
	//update w[X] to incorporate NP interaction term
	// wX += exp_nr_chiXPN / N * rho_exp_nr * C * N
	//     = exp_nr_chiXPN * rho_exp_nr * C (N's cancel out)
	int is1 = intSpecies[0] ;
	int is2 = intSpecies[1] ;
	

	// fill vector with C*chiN

	thrust::device_vector<thrust::complex<double>> chiC(mybox->M);
	thrust::fill(chiC.begin(), chiC.end(), -1.0 * chiN * mybox->C);


	thrust::device_vector<thrust::complex<double>> dtmp(mybox->M);

	// multiply density by chiN*C, store in temp variable
	thrust::transform(chiC.begin(), chiC.end(),  mybox->Species[is2].d_density.begin(), dtmp.begin(), thrust::multiplies<thrust::complex<double>>());
	
	
	//assign wpl
	thrust::transform(d_wpl.begin(), d_wpl.end(), dtmp.begin(), d_wpl.begin(), thrust::plus<thrust::complex<double>>());

}

void PotentialParticle::initLinearCoeffs() {
}

void PotentialParticle::writeFields(int potInd) {
	char nm[30];
	sprintf(nm, "wpl_Particle%d.dat", potInd);
	
	wpl = d_wpl;
	mybox->writeTComplexGridData(nm, wpl);

}

void PotentialParticle::updateFields() {
}

std::complex<double> PotentialParticle::calcHamiltonian() {
	Hterm = 0 ;
	return Hterm ;
}

PotentialParticle::~PotentialParticle() {}
