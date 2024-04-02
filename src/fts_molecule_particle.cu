#include "globals.h"
#include "fts_molecule_particle.h"
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/complex.h>
#include <thrust/copy.h>
#include "FTS_Box.h"
#include "fts_species.h"
#include "fts_potential.h"
#include <sstream>
#include <string>
#include <iostream>
void die(const char*);

ParticleMolec::~ParticleMolec(){}

ParticleMolec::ParticleMolec(std::istringstream& iss, FTS_Box* p_box) : FTS_Molec(iss, p_box) {
 // iss comes into this routine have already passed "molecule" and "particle"
 // structure of iss should be: 
 // particleNum, int, total number of particles 
 // particleSpecies, char, should be a species defined earlier in input file
 // Rp, float, for particle 1
 // xi, float, for particle 1
 // center x, y, z positions, Dim floats, for particle 1. Should be a number between 0 and 1 which will be scaled to the box length. 
 // Repeat for number of particles

 //we are having the user optionally put in a file which has 
 // column 1 = particle number
 // column 2 = NP radius 
 // column 3 = NP xi 
 //column 4 = x center position 
 // column 5 = y center position 
 // column 6 = z center position 
  
 
 //Later we can implement a particleType for field-based or nanorods
 // Right now I'm just goin to do explicit NPs

	iss >> particleNum;
	std::cout << "particleNum=" << particleNum <<std::endl;
	iss >> particleSpecies;
	std::cout << "particleSpecies=" << particleSpecies <<std::endl;  
	// determine integer species
	intSpecies.resize(1);
	d_intSpecies.resize(1);
	for (int i = 0; i<mybox->Species.size(); i++ ) {
		if ( particleSpecies == mybox->Species[i].fts_species ) {
			intSpecies[0] = i;
		}
	}
	std::cout << "integer species= " << intSpecies[0] <<std::endl; 
	d_intSpecies = intSpecies;

	//resize density arrays
	density.resize(mybox->M);
	std::cout << "density array resized" << std::endl;
	d_density.resize(mybox->M);
        d_NPdensity.resize(mybox->M); 
	 
	//need to resize R, xi, center arrays based on particleNum
	Rp.resize(particleNum);
	std::cout << "rp array resized" << std::endl;
	d_Rp.resize(particleNum);

	xi.resize(particleNum);
	std::cout << "xi array resized" << std::endl;   
	d_xi.resize(particleNum);

        center.resize(particleNum*mybox->returnDimension());
	std::cout << "center array resized" << std::endl;   
	d_center.resize(particleNum*mybox->returnDimension());
	
	Vnptot = 0; //zero total NP volume 
	//read input file 
	std::string s1;
	iss >> s1;
	std::cout << "looking for file" <<std::endl;
	if (s1 == "file") {
		iss >> s1;
		std::cout << "Filename: " << s1 <<std::endl;
		std::ifstream in2(s1);
		if (not in2.is_open()){
			std::cout << "File" << s1 << " does not exist."<<std::endl;
			die("");
		}
		std::cout << "opening file: " << s1 <<std::endl;  
		// Store the contents into a vector of strings
		int npCounter = 0;
			
		std::string line;
		while (std::getline(in2, line)) {
			std::istringstream iss(line);
			std::string word;
    			std::vector<std::string> outputs;
			while (iss >> word) {
				outputs.push_back(word);
			}
			Rp[npCounter] = stof(outputs[1]);
			xi[npCounter] = stof(outputs[2]);
			center[3*npCounter] = stof(outputs[3]);
			center[(3*npCounter)+1] = stof(outputs[4]);
			center[(3*npCounter)+2] = stof(outputs[5]);
			npCounter += 1;
			}
		}

	//loop to compute volume 
	for (int j=0; j<particleNum; j++ ) {
		std::cout << "Particle Number " << j << std::endl;
		std::cout << "R = " << Rp[j] << std::endl;
		std::cout << "xi = " << xi[j] << std::endl;
		std::cout << "x center = " << center[(3*j)] << std::endl;
		std::cout << "y center = " << center[(3*j)+1] << std::endl;
		std::cout << "z center = " << center[(3*j)+2] << std::endl;
		float Vnp = ( 4 / 3) * PI * Rp[j] * Rp[j] * Rp[j]; //volume of 1 particle
	
		Vnptot += Vnp; // sum volume of all particles		
	}
	// copy center positions, radii, xi to device (necessary?)
	d_center = center;
	d_Rp = Rp;
	d_xi = xi;
	
	//calculate particle volume fraction
	phiNP = Vnptot / mybox->V;
	
	std::cout << "phiNP  = " << phiNP << std::endl;
	//update free volume available in box
	mybox->Vfree -= Vnptot;			

	// here we want to use our erfc function to calculate the density

	thrust::fill(d_NPdensity.begin(), d_NPdensity.end(), 0.0);
	// Loop over particles
	for (int j = 0; j < particleNum; j++ ) {

		// Loop over grid points
		for (int i=0; i < mybox->M; i++) {
			double r[mybox->returnDimension()];
			mybox->get_r(i, r); // gives position in each direction based on grid point
			double mdr2 = 0;
			double dr_abs;
			double dr[mybox->returnDimension()];
			// Loop over dimensions
			for (int k = 0; k < mybox->returnDimension(); k++) {
				//calculate distance from NP center
				dr[k] = center[(3*j)+k] - r[k];
				//take into account periodic boundaries
				if (dr[k] >= 0.5 * mybox-> L[k]) dr[k] -= mybox->L[k];
				else if (dr[k] < -0.5 * mybox -> L[k]) dr[k] += mybox->L[k];
				mdr2 += dr[k] * dr[k]; 
			}
			dr_abs = sqrt(mdr2);
			density[i] += mybox->Nr * mybox->rho0 * 0.5 * erfc( ( dr_abs-Rp[j] ) / xi[j] );	
}	//Zero NP density field
}
	//transfer density to device
	d_NPdensity = density;
	int is = intSpecies[0];
	// Also need to accumulate density onto the relevant species fields
	thrust::transform(d_NPdensity.begin(), d_NPdensity.end(),  mybox->Species[is].d_density.begin(), mybox->Species[is].d_density.begin(), thrust::plus<thrust::complex<double>>());
}
void ParticleMolec::calcDensity() {
	d_density=d_NPdensity;			
 	int is = intSpecies[0]; 
	thrust::transform(d_density.begin(), d_density.end(),  mybox->Species[is].d_density.begin(), mybox->Species[is].d_density.begin(), thrust::plus<thrust::complex<double>>());
	calcHamiltonian();
} 


// here we will calculate the Hamiltonian term which incorporates the NP density
// - I * C * int (wpl(r) * ( - rhoNP(r)))

void ParticleMolec::calcHamiltonian() {
	thrust::device_vector<thrust::complex<double>> dtmp(mybox->M);
	thrust::complex<double> I(0.0, 1.0);
	thrust::device_vector<thrust::complex<double>> dtmp2(mybox->M);

	//dtmp(r) = wpl(r) * (- rhoNP(r))

	// first, create  - rhoNP(r)....

	//filling vector with 1
	thrust::device_vector<float> V1(mybox->M);
	thrust::fill(V1.begin(), V1.end(), -1.0);

	//multiplying -1*rhoNP, storing in dtmp2
	thrust::transform(V1.begin(), V1.end(), d_NPdensity.begin(), dtmp2.begin(), thrust::multiplies<thrust::complex<double>>()); 
	//then multiply wpl(r) * (- rhoNP(r)) = d_wpl * dtmp2, storing in dtmp
	int ip;
	// find Helfand potential to get wpl
	for (int i = 0; i < mybox->Potentials.size(); i++ ) {
		if ( mybox->Potentials[i]->printStyle() == "Helfand" ) {
			ip = i;
		}
	}		
	thrust::transform(mybox ->Potentials[ip]->d_wpl.begin(),mybox -> Potentials[ip]->d_wpl.end(), dtmp2.begin(), dtmp.begin(), thrust::multiplies<thrust::complex<double>>());
	
	// integrate int (wpl(r) * (rhoNP(r)))
	thrust::complex<double> integral = thrust::reduce(dtmp.begin(), dtmp.end()) * mybox->gvol;

	// -i*C*int
	Hterm = -I * mybox->C * integral;


} 



void ParticleMolec::computeLinearTerms() {

}
