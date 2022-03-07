#include <iostream>
#include <math.h>
#include <gwat/fisher.h>
#include <gwat/IMRPhenomP.h>
#include <gwat/util.h>
#include <gwat/waveform_util.h>
#include <gwat/io_util.h>
#include <gwat/ortho_basis.h>
#include <gwat/detector_util.h>
#include <gsl/gsl_rng.h>

int fisher_batch(int argc, char *argv[]);
void ERROR_MESSAGE();

int main(int argc, char *argv[]){
	// if(argc != 2){
	// 	ERROR_MESSAGE();
	// 	return -1;
	// }
	int option = atoi(argv[1]);	
	if(option == 0){
		return fisher_batch(argc, argv);
	}
}

int fisher_batch(int argc, char *argv[])
{	
	int n_inj = atoi(argv[2]);
	std::string path_injection = "input/injection.csv";
	std::string path_fisher = "output/fisher.csv";
	std::string path_snr = "output/snr.csv";
	std::cout << "Num of Injections: " << n_inj << std::endl;

	/*--------------------*/
	int n_det = 3;
	std::string detector_moniker[n_det] = {"H","L","V"};
	std::string detectors[n_det] = {"Hanford","Livingston","Virgo"};
	std::string psd_names[n_det] = {"AdLIGOAPlus_smoothed","AdLIGOAPlus_smoothed","AdVIRGOPlus2_opt"};
	std::string wf_model_name = "IMRPhenomD_NRT";

	double FLOW = 5;
	double FHIGH = 2048;

	gen_params gp;
	gp.tidal_love = false;
	gp.shift_time = false;
	gp.shift_phase = true;
	gp.f_ref = FLOW;
	std::cout << "Waveform: " << wf_model_name << std::endl;
	std::cout << "Ref Freq: " << gp.f_ref << std::endl;


	/*--------------------*/
	int n_len = 5000;
	double frequency[n_len];
	double weights[n_len];
	gauleg(log10(FLOW), log10(FHIGH), frequency, weights, n_len);
	for (int i = 0; i < n_len; i++){
		frequency[i] = pow(10., frequency[i]);
	}

	double psds[n_det][n_len];
	for (int i = 0; i < n_det; i++) {
		populate_noise(frequency, psd_names[i], psds[i], n_len, 48);
		for (int k = 0; k < n_len; k++) {
			psds[i][k] *= psds[i][k];
		}
	}

	/*--------------------*/
	// Injection array: ra, cos(dec), psi, cos(iota), phic, tc, ln(DL), ln(Mc), eta, chi1z, chi2z, lam1, lam2
	// This is the same order as in the output Fisher mat
	// The output Fisher mat is flattened
	int n_dim = 13;
	double **input_inj = new double* [n_inj];
	for (int i = 0; i < n_inj; i++) {
		input_inj[i] = new double [n_dim];
	}
	read_file(path_injection, input_inj, n_inj, n_dim);

	double **output_fisher = new double* [n_inj];
	for (int i = 0; i < n_inj; i++) {
		output_fisher[i] = new double [n_dim*n_dim];
	}
	double *output_snr = new double [n_inj];

	/*--------------------*/
	for (int i = 0; i < n_inj; i++) {

		std::cout << "Processing " << to_string(i) << " / " << to_string(n_inj);

		gp.RA = input_inj[i][0];
		gp.DEC = asin(input_inj[i][1]);
		gp.psi = input_inj[i][2];
		gp.incl_angle = acos(input_inj[i][3]);
		gp.phiRef = input_inj[i][4];
		gp.tc = 0;
		gp.gmst = gps_to_GMST_radian(input_inj[i][5]);
		gp.Luminosity_Distance = exp(input_inj[i][6]);
		double Mc = exp(input_inj[i][7]);
		double eta = input_inj[i][8];
		gp.mass1 = calculate_mass1(Mc, eta);
		gp.mass2 = calculate_mass2(Mc, eta);
		gp.spin1[2] = input_inj[i][9];
		gp.spin2[2] = input_inj[i][10];
		gp.tidal1 = input_inj[i][11];
		gp.tidal2 = input_inj[i][12];
		gp.tidal_s = .5 * (gp.tidal1 + gp.tidal2);

		output_snr[i] = 0;
		for (int j = 0; j < n_det; j++) {
			double times;
			std::complex<double> *response = new std::complex<double>[n_len];
			fourier_detector_response(
				frequency, n_len, response, 
				detectors[j], wf_model_name, &gp, &times);
			double snr = calculate_snr_internal(
				psds[j], response, frequency, n_len, 
				"GAUSSLEG", weights, true);
			delete [] response;
			output_snr[i] += snr * snr;
		}
		output_snr[i] = sqrt(output_snr[i]);

		std::cout << ", SNR = " << to_string(output_snr[i]) << " ...";

		for(int k = 0; k < n_dim*n_dim; k++){
			output_fisher[i][k] = 0;
		}
		for (int j = 0; j < n_det; j++) {

			double **fisher = new double* [n_dim];
			for (int k = 0 ; k < n_dim; k++) {
				fisher[k] = new double [n_dim];		
			}

			fisher_autodiff(
				frequency, n_len, 
				wf_model_name, 
				detectors[j], detectors[0], 
				fisher, n_dim, &gp, 
				"GAUSSLEG", weights, true, psds[j], 
				(int *) NULL, (int *)NULL);

			double fac1;
			double fac2;
			for (int k = 0; k < n_dim; k++){
				fac1 = 1;
				if (k == 1 || k == 3) {
					fac1 /= sqrt(1 - input_inj[i][k]*input_inj[i][k]);
					if (k == 3) fac1 = -fac1;
				}
				for(int l = 0; l < n_dim; l++){
					fac2 = 1;
					if (l == 1 || l == 3) {
						fac2 /= sqrt(1 - input_inj[i][l]*input_inj[i][l]);
						if (l == 3) fac2 = -fac2;
					}
					output_fisher[i][k*n_dim+l] += fisher[k][l] * fac1 * fac2;
				}
			}

			for(int k = 0 ; k < n_dim; k++){
				delete [] fisher[k];
			}	
			delete [] fisher;
		}

		std::cout << " Fisher done!" << std::endl;
	}

	/*--------------------*/
	std::cout << "Writing Output...";

	write_file(path_snr, output_snr, n_inj);
	write_file(path_fisher, output_fisher, n_inj, n_dim*n_dim);

	std::cout << "done!" << std::endl;

	for (int i = 0; i < n_inj; i++) {
		delete [] input_inj[i];
		delete [] output_fisher[i];
	}
	delete [] input_inj;
	delete [] output_fisher;
	delete [] output_snr;
	
	return 0;

}

void ERROR_MESSAGE(){
	std::cout<<"ERROR -- 1 option:"<<std::endl;
	std::cout<<"0 -- Fully random"<<std::endl;
	std::cout<<"1 -- targeted "<<std::endl;
	std::cout<<"2 -- DL/iota investigation "<<std::endl;
	std::cout<<"3 -- Single injection "<<std::endl;
	return;
}
