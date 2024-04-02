// Copyright (c) 2023 University of Pennsylvania
// Part of MATILDA.FT, released under the GNU Public License version 2 (GPLv2).


#include "include_libs.h"
#include "fts_potential_helfand.h"
#include "FTS_Box.h"
#include "fts_species.h"
//#include "globals.h"

void die(const char*);
double ran2();


PotentialHelfand::PotentialHelfand(std::istringstream& iss, FTS_Box* p_box) : FTS_Potential(iss, p_box) {


    // Set stringstream to be ready to read kappa
    iss.seekg(0);
    std::string s1;
    iss >> s1;
    iss >> s1;

    potentialStyle = "Helfand";

    iss >> kappaN;
    iss >> delt;

    double ivalue = 0.0;
    wpl.resize(mybox->M,ivalue);
    d_wpl.resize(mybox->M, ivalue);
    d_Akpl.resize(mybox->M, ivalue);

    // Set default update scheme
    updateScheme = "EM";

    while (iss.tellg() != -1 ) {
        iss >> s1;
        if ( s1 == "initialize" ) {
            std::cout << "caught initialize!" << std::endl;
            iss >> s1;
            if ( s1 == "value" ) {
                double rVal, iVal;
                iss >> rVal;
                iss >> iVal;
                thrust::fill(wpl.begin(), wpl.end(), std::complex<double>(rVal, iVal));
                d_wpl = wpl;
            }
	    // User expected to input name of file containing Helfand wpl field to resume from
	    else if (s1 == "resume") { 
		iss >> s1;
		std::ifstream in2(s1);
		if (not in2.is_open()){
			std::cout << "File " << s1 << " does not exist." << std::endl;
			die("");
	    	}
	    // reads input file to set wpl
	    // expects x, y, z, rVal, iVal in file columns
	    // currently assumes each row is a new grid point going from 0->M
	    // does not actually use the xyz values
	    	double x;
	    	double y;
	    	double z;
	    	double rVal;
	    	double iVal;
	    	for (int i = 0; i<mybox->M; i++) {
			in2 >> x >> y >> z >> rVal >> iVal;
			wpl[i] = std::complex<double>(rVal, iVal);
		}
		d_wpl = wpl;
	    } 
            // Two floats expected: amplitude of noise on real part and imag part
            else if ( s1 == "random" ) {
                double rAmp, iAmp;
                iss >> rAmp;
                iss >> iAmp;
                // Fill host field with random noise
                for ( int i=0 ; i<mybox->M ; i++ ) {
                    wpl[i] = std::complex<double>(rAmp * ran2(), iAmp * ran2() );
                }
                
                // transfer to device
                d_wpl = wpl;
            }
            
            // Expects an int and two doubles [int dir] [double amplitude] [double period]
            else if ( s1 == "sin" || s1 == "sine" ) {
                double amp, period;
                int dir;
                iss >> dir;
                iss >> amp;
                iss >> period; 

                std::complex<double> I(0.0,1.0);
                for ( int i=0 ; i<mybox->M ; i++ ) {
                    double r[3];
                    mybox->get_r(i, r);
                    wpl[i] = I * amp * sin(2.0 * PI * r[dir] * period / mybox->L[dir]);
                }

                // transer to device
                d_wpl = wpl;
            }

            else {
                die("Invalid initialize option on potential helfand");
            }
        }

        else if ( s1 == "updateScheme" ) {
            iss >> updateScheme;
        }

        else if ( s1 == "modify" ) {
            iss >> s1;
            if ( s1 != "zeromean" && s1 != "zeroMean" ) { die("Invalid modify option on Helfand potential"); }
            zeroMean = true;
        }
    }// optional arguments



}

void PotentialHelfand::updateFields() {

    // Construct total density field
    // This field should contain the *smeared* density fields
    thrust::device_vector<thrust::complex<double>> d_rho_total(mybox->M);

    // Initialize to zero
    thrust::fill(d_rho_total.begin(), d_rho_total.end(), 0.0);

    
    // Loop over species, adding them to the field
    for ( int i=0 ; i<mybox->Species.size() ; i++ ) {
        thrust::transform(mybox->Species[i].d_density.begin(), mybox->Species[i].d_density.begin()+mybox->M,
            d_rho_total.begin(), d_rho_total.begin(), thrust::plus<thrust::complex<double>>());
    }

    // Vector to store the force term
    thrust::device_vector<thrust::complex<double>> d_dHdw(mybox->M);

    // cast thrust vectors to cuDoubleComplex for use in kernel
    cuDoubleComplex* _d_dHdw = (cuDoubleComplex*)thrust::raw_pointer_cast(d_dHdw.data());
    cuDoubleComplex* _d_wpl = (cuDoubleComplex*)thrust::raw_pointer_cast(d_wpl.data());
    cuDoubleComplex* _d_rho_total = (cuDoubleComplex*)thrust::raw_pointer_cast(d_rho_total.data());

    // Make the force in real space
    d_makeHelfandForce<<<mybox->M_Grid, mybox->M_Block>>>(_d_dHdw, _d_wpl, _d_rho_total, mybox->C,
        kappaN, mybox->Nr, mybox->M);


    // Update the fields
    if ( updateScheme == "EM" ) {
        d_fts_updateEM<<<mybox->M_Grid, mybox->M_Block>>>(_d_wpl, _d_dHdw, delt, mybox->M);
    }


    else if ( updateScheme == "1S" ) {
        // Put the force and potential into k-space
        mybox->cufftWrapperDouble(d_dHdw, d_dHdw, 1);
        mybox->cufftWrapperDouble(d_wpl, d_wpl, 1);

        // Pointer to linear coefficient
        cuDoubleComplex* _d_Ak = (cuDoubleComplex*)thrust::raw_pointer_cast(d_Akpl.data());

        // Call updater
        d_fts_update1S<<<mybox->M_Grid, mybox->M_Block>>>(_d_wpl, _d_dHdw, _d_Ak, delt, mybox->M);

        // Bring potential back to r-space
        mybox->cufftWrapperDouble(d_wpl, d_wpl, -1);
    }

    // Check for modifiers
    if ( zeroMean == true ) {
        thrust::complex<double> mean = thrust::reduce(d_wpl.begin(), d_wpl.end()) / double(mybox->M);
        
        // dtmp = mean
        thrust::device_vector<thrust::complex<double>> dtmp(mybox->M, mean);

        // wpl(r) = wpl(r) - mean
        thrust::transform(d_wpl.begin(), d_wpl.end(), dtmp.begin(), d_wpl.begin(), 
            thrust::minus<thrust::complex<double>>());
    }

}


// This routine is currently written to deal with dHdw in real space
// the I*rho0 term should change if it changed to k-space updating
__global__ void d_makeHelfandForce(
    cuDoubleComplex* dHdw,              // Field holding dHdw
    const cuDoubleComplex* w,           // current d_wpl
    const cuDoubleComplex* rho_total,   // current total density
    const double C,                     // Polymer concentration, based on Nr
    const double kN,                    // kappa * N
    const double Nr,                    // Reference chain length
    const int M                         // number of grid points
    ) {

    const int ind = blockIdx.x * blockDim.x + threadIdx.x;
    if (ind >= M)
        return;

    // dHdw = C * w / kN - I * C + I * rho_total / Nr ;
    dHdw[ind].x = C * w[ind].x / kN - rho_total[ind].y / Nr;
    dHdw[ind].y = C * w[ind].y / kN + rho_total[ind].x / Nr - C;
}


void PotentialHelfand::writeFields(int potInd ) { 
    char nm[30];
    sprintf(nm, "wpl_Helfand%d.dat", potInd);

    // Transfer field to the host;
    wpl = d_wpl;
    mybox->writeTComplexGridData(nm, wpl);
}


// Computes this potential's contribution to the effective Hamiltonian
std::complex<double> PotentialHelfand::calcHamiltonian() {
    thrust::device_vector<thrust::complex<double>> dtmp(mybox->M);
    thrust::complex<double> I(0.0,1.0);

    // dtmp(r) = wpl(r)^2
    thrust::transform(d_wpl.begin(), d_wpl.end(), d_wpl.begin(), dtmp.begin(), 
        thrust::multiplies<thrust::complex<double>>());

    thrust::complex<double> integral = thrust::reduce(dtmp.begin(), dtmp.end()) * mybox->gvol;

    Hterm = integral * mybox->C / 2.0 / kappaN;

    // -i C * int(wpl)
    integral = thrust::reduce(d_wpl.begin(), d_wpl.end()) * mybox->gvol;

    Hterm += -I * mybox->C * integral;

    //std::cout << Hterm << std::endl;

    return Hterm;
    
}

void PotentialHelfand::initLinearCoeffs() {
    // Akpl = C / kappaN
    thrust::fill(d_Akpl.begin(), d_Akpl.end(), mybox->C/kappaN);
}

PotentialHelfand::~PotentialHelfand() {}
