//Cahn-Hilliard Phase-Field Model to investigate decomposition of two-dimensional domain
//Grid spacing, h = 1 in both directions. Periodic boundary conditions
//Square domain of 128 x 128 size on Cartesian grid

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <ctime>
#include <omp.h>

int main(int argc, char* argv[])
{
	srand(time(NULL));
	//Defining the maximum domain size in X and Y direction
	int imax = 128, jmax = 128;

	//Defining the Grid Spacing
	int h = 1;

	//Defining the fixed parameters
	int a = 10, b = 1;

	//Defining time step
	double dt = 2 * pow(10, -3);

	//The simulation is to be executed to tmax = 500
	int tmax = 500;
	double t = tmax / dt;
	int k;

	std::vector<double> m(t, 0.0);

	for (int i = 1; i < t; i++)
	{
		m[i] = m[i - 1] + dt;
	}

	std::ofstream write_m("m.dat");
	write_m.setf(std::ios::scientific);
	write_m.precision(16);
	assert(write_m.is_open());
	for (int i = 0; i < t; i++)
	{
		write_m << m[i] << "\n";
	}

	double p[] = { 0.49, 0.51 };

	//Defining the domain phi
	std::vector<std::vector<double> > phi(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double> > dphi_x(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double> > dphi_y(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double> > d2phi_x(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double> > d2phi_y(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double> > gdash(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double> > mu(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double> > d2mu_x(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double> > d2mu_y(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double> > phi_new(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double> > div_phi(imax, std::vector<double>(jmax, 0.0));
	std::vector<double> rc(t, 0.0);
	std::vector<double> IE(t, 0.0);
	std::vector<double> ME(t, 0.0);
	std::vector<double> TE(t, 0.0);
	std::vector<double> int_area(t, 0.0);
	std::vector<double> int_length(t, 0.0);

	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
			phi[i][j] = (p[1] - p[0])*((double)rand()/RAND_MAX) + p[0];
		}
	}

#pragma omp parallel for private(k) shared(t) schedule(dynamic, 10000)
	//Running the simulation and marching in time
	for (k = 0; k < 250000; k++)
	{
		//Finding out the Laplacian of phi
		for (int j = 0; j < jmax; j++)
		{
			for (int i = 1; i < imax - 1; i++)
			{
				d2phi_x[i][j] = (phi[i + 1][j] - 2 * phi[i][j] + phi[i - 1][j]) / (h*h);
			}
			d2phi_x[0][j] = (phi[1][j] - 2 * phi[0][j] + phi[imax-1][j]) / (h*h);
			d2phi_x[imax - 1][j] = (phi[0][j] - 2 * phi[imax - 1][j] + phi[imax - 2][j]) / (h*h);
		}


		for (int i = 0; i < imax; i++)
		{
			for (int j = 1; j < jmax - 1; j++)
			{
				d2phi_y[i][j] = (phi[i][j + 1] - 2 * phi[i][j] + phi[i][j - 1]) / (h*h);
			}
			d2phi_y[i][0] = (phi[i][1] - 2 * phi[i][0] + phi[i][jmax - 1]) / (h*h);
			d2phi_y[i][jmax - 1] = (phi[i][0] - 2 * phi[i][jmax - 1] + phi[i][jmax - 2]) / (h*h);
		}

		//Calculating out g'(phi)
		for (int i = 0; i < imax; i++)
		{
			for (int j = 0; j < jmax; j++)
			{
				gdash[i][j] = 32 * b*(2 * pow(phi[i][j], 3) - 3 * pow(phi[i][j], 2) + phi[i][j]);
			}
		}

		//Calculating the values of mu
		for (int i = 0; i < imax; i++)
		{
			for (int j = 0; j < jmax; j++)
			{
				mu[i][j] = gdash[i][j] - a*(d2phi_x[i][j] + d2phi_y[i][j]);
			}
		}

		//Calculating the Laplacian of mu
		for (int j = 0; j < jmax; j++)
		{
			for (int i = 1; i < imax - 1; i++)
			{
				d2mu_x[i][j] = (mu[i + 1][j] - 2 * mu[i][j] + mu[i - 1][j]) / (h*h);
			}
			d2mu_x[0][j] = (mu[1][j] - 2 * mu[0][j] + mu[imax - 1][j]) / (h*h);
			d2mu_x[imax - 1][j] = (mu[0][j] - 2 * mu[imax - 1][j] + mu[imax - 2][j]) / (h*h);
		}
		for (int i = 0; i < imax; i++)
		{
			for (int j = 1; j < jmax - 1; j++)
			{
				d2mu_y[i][j] = (mu[i][j + 1] - 2 * mu[i][j] + mu[i][j - 1]) / (h*h);
			}
			d2mu_y[i][0] = (mu[i][1] - 2 * mu[i][0] + mu[i][jmax - 1]) / (h*h);
			d2mu_y[i][jmax - 1] = (mu[i][0] - 2 * mu[i][jmax - 1] + mu[i][jmax - 2]) / (h*h);
		}

		//Calculating the new values of phi
		for (int i = 0; i < imax; i++)
		{
			for (int j = 0; j < jmax; j++)
			{
				phi_new[i][j] = phi[i][j] + dt*(d2mu_x[i][j] + d2mu_y[i][j]);
			}
		}

		//Calculating the Interfacial Energy, Mixing Energy and Total Energy
		//Finding out the divergence of phi
		for (int j = 0; j < jmax; j++)
		{
			for (int i = 1; i < imax - 1; i++)
			{
				dphi_x[i][j] = (phi[i + 1][j] - phi[i - 1][j]) / (2 * h);
			}
			dphi_x[0][j] = (phi[1][j] - phi[imax - 1][j]) / (2 * h);
			dphi_x[imax - 1][j] = (phi[0][j] - phi[imax - 2][j]) / (2 * h);
		}
		for (int i = 0; i < imax; i++)
		{
			for (int j = 1; j < jmax - 1; j++)
			{
				dphi_y[i][j] = (phi[i][j + 1] - phi[i][j - 1]) / (2 * h);
			}
			dphi_y[i][0] = (phi[i][1] - phi[i][jmax - 1]) / (2 * h);
			dphi_y[i][jmax - 1] = (phi[i][0] - phi[i][jmax - 2]) / (2 * h);
		}

		for (int i = 0; i < imax; i++)
		{
			for (int j = 0; j < jmax; j++)
			{
				div_phi[i][j] = dphi_x[i][j] + dphi_y[i][j];
			}
		}

		for (int j = 0; j < jmax; j++)
		{
			for (int i = 0; i < imax; i++)
			{
				IE[k] = IE[k] + 0.5*a*pow(div_phi[i][j], 2);
				ME[k] = ME[k] + 16 * b*((pow((phi[i][j] - 1), 2))*pow(phi[i][j], 2));
				TE[k] = TE[k] + 0.5*a*pow(div_phi[i][j], 2) + 16 * b*((pow((phi[i][j] - 1), 2))*pow(phi[i][j], 2));
				int_area[k] = int_area[k] + phi[i][j];
				int_length[k] = int_length[k] + abs(sqrt(div_phi[i][j]));
			}
		}

		//Calculating the characteristic length 
		rc[k] = int_area[k] / int_length[k];

		//Updating the values of phi
		for (int j = 0; j < jmax; j++)
		{
			for (int i = 0; i < imax; i++)
			{
				phi[i][j] = phi_new[i][j];
			}
		}

	}

	//Writing the data for plots

	std::ofstream write_phi("phi.dat");
	write_phi.setf(std::ios::scientific);
	write_phi.precision(16);
	assert(write_phi.is_open());
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
			write_phi << phi[i][j] << " ";
		}
		write_phi << "\n";
	}
	write_phi.close();

	std::ofstream write_IE("IE.dat");
	write_IE.setf(std::ios::scientific);
	write_IE.precision(16);
	assert(write_IE.is_open());
	for (int i = 0; i < t; i++)
	{
		write_IE << IE[i] << " ";
	}
	write_IE.close();

	std::ofstream write_rc("rc.dat");
	write_rc.setf(std::ios::scientific);
	write_rc.precision(16);
	assert(write_rc.is_open());
	for (int i = 0; i < t; i++)
	{
		write_rc << rc[i] << " ";
	}
	write_rc.close();

	std::ofstream write_ME("ME.dat");
	write_ME.setf(std::ios::scientific);
	write_ME.precision(16);
	assert(write_ME.is_open());
	for (int i = 0; i < t; i++)
	{
		write_ME << ME[i] << " ";
	}
	write_ME.close();

	std::ofstream write_TE("TE.dat");
	write_TE.setf(std::ios::scientific);
	write_TE.precision(16);
	assert(write_TE.is_open());
	for (int i = 0; i < t; i++)
	{
		write_TE << TE[i] << " ";
	}
	write_TE.close();

	return 0;
}
