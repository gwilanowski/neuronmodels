#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

double z = 0.035; //injected current
double c[7]={0, 7.5e-3, 7.5e-3, 7.5e-3, 4.98e-4, 3.14e-5, 0}; //total capacitance
double gl[7]={0, 1.52e-3, 1.52e-3, 1.52e-3, 9.96e-5, 3.14e-5, 0}; //total leak conductance
double gam[7][7] = {
	{ 0,        0,        0,       0,        0,        0,       0},		
	{ 0,        0,  2.68e-3,       0,        0,        0,       0},		
	{ 0,        0,        0, 2.68e-3,        0,        0,       0},	//copuling conductances	
	{ 0,        0,        0,       0,  4.47e-3,        0,       0},		
	{ 0,        0,        0,       0,        0,  2.07e-3,       0},		
	{ 0,        0,        0,       0,        0,        0,       0},		
	{ 0,        0,        0,       0,        0,        0,       0}		
};
double omega[7]; //compartment-specific currents

ofstream out_file("traub.dat"); //file with results

typedef boost::array< double , 13 > state_type;

//large MN conductances
double gks_tot = 0.0075;
double gk_s_tot = 0.01;
double gna_s_tot = 0.05;
double gk_is_tot = 0.0031;
double gna_is_tot = 0.019;

double gks, gk_s, gna_s, gk_is, gna_is, m_s, h_s, n_s, q, m_is, h_is, n_is, a_m_s, a_h_s, a_n_s, a_m_is, a_h_is, a_n_is, b_m_s,
b_h_s, b_n_s, b_m_is, b_h_is, b_n_is, a_q, b_q;


/* // state vector
x[0]  
x[1] v1
x[2] v2
x[3] v3
x[4] v4
x[5] v5
x[6] m_s
x[7] h_s
x[8] n_s
x[9] q
x[10] m_is
x[11] h_is
x[12] n_is
*/

//Right rand side function of ode x'=f(x)
void traub( const state_type &x , state_type &dxdt , double t )
{
	gks = gks_tot*pow(x[9],2);
	gk_s = gk_s_tot*pow(x[8],4);
	gna_s = gna_s_tot*pow(x[6],3)*x[7];
	gk_is = gk_is_tot*pow(x[12],4);
	gna_is = gna_is_tot*pow(x[10],3)*x[11];

	omega[1] = 0;
	omega[2] = 0;
	omega[3] = 0;
	omega[4] = z -(gks+gk_s)*(x[4]+5)-gna_s*(x[4]-115);
    omega[5] = -gk_is*(x[5]+5)-gna_is*(x[5]-115);

	for (int ix = 1; ix <= 5; ++ix) 
			dxdt[ix]=(- gl[ix]*x[ix] + gam[ix-1][ix]*(x[ix-1]-x[ix]) + gam[ix][ix+1]*(x[ix+1]-x[ix]) + omega[ix])/c[ix];

	a_m_s = 0.4*(17.5-x[4])/(exp((17.5-x[4])/5)-1);
	b_m_s = 0.4*(x[4]-45)/(exp((x[4]-45)/5)-1);
	dxdt[6] = a_m_s*(1-x[6])-b_m_s*x[6];	

	a_h_s = 0.28*exp((25-x[4])/20);
	b_h_s = 4/(exp((40-x[4])/10)+1);
	dxdt[7] = a_h_s*(1-x[7])-b_h_s*x[7];	
	
	a_n_s = 0.02*(20-x[4])/(exp((20-x[4])/10)-1);
	b_n_s = 0.25*exp((10-x[4])/80);
	dxdt[8] = a_n_s*(1-x[8])-b_n_s*x[8];
	
	a_q = 3.5/(exp((55-x[4])/4)+1);
	b_q = 0.03;
	dxdt[9] = a_q*(1-x[9])-b_q*x[9];	
	
	a_m_is = 0.4*(10-x[5])/(exp((10-x[5])/5)-1);
	b_m_is = 0.4*(x[5]-35)/(exp((x[5]-35)/5)-1);
	dxdt[10] = a_m_is*(1-x[10])-b_m_is*x[10];	

	a_h_is = 0.28*exp((30-x[5])/20);
	b_h_is = 4/(exp((30-x[5])/10)+1);
	dxdt[11] = a_h_is*(1-x[11])-b_h_is*x[11];	
	
	a_n_is = 0.02*(10-x[5])/(exp((10-x[5])/10)-1);
	b_n_is = 0.25*exp(-x[5]/80);
	dxdt[12] = a_n_is*(1-x[12])-b_n_is*x[12];	

}
//write results to file
void write_traub( const state_type &x , const double t )
{
	static bool first = true;
	if (first)
	{
		out_file << "t" << '\t' << "v_s" << endl;
		first = false;
	}
	else
		out_file << t << '\t' << x[4] << endl;
}

int main(int argc, char **argv)
{
	//setting initial conditions
	state_type x = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	a_m_s = 0.4*(17.5-x[4])/(exp((17.5-x[4])/5)-1);
	b_m_s = 0.4*(x[4]-45)/(exp((x[4]-45)/5)-1);
	x[6] = a_m_s/(a_m_s+b_m_s);	

	a_h_s = 0.28*exp((25-x[4])/20);
	b_h_s = 4/(exp((40-x[4])/10)+1);
	x[7] = a_h_s/(a_h_s + b_h_s);	
	
	a_n_s = 0.02*(20-x[4])/(exp((20-x[4])/10)-1);
	b_n_s = 0.25*exp((10-x[4])/80);
	x[8] = a_n_s/(a_n_s + b_n_s);
	
	a_q = 3.5/(exp((55-x[4])/4)+1);
	b_q = 0.03;
	x[9] = a_q/(a_q + b_q);	
	
	a_m_is = 0.4*(10-x[5])/(exp((10-x[5])/5)-1);
	b_m_is = 0.4*(x[5]-35)/(exp((x[5]-35)/5)-1);
	x[10] = a_m_is/(a_m_is + b_m_is);

	a_h_is = 0.28*exp((30-x[5])/20);
	b_h_is = 4/(exp((30-x[5])/10)+1);
	x[11] = a_h_is/(a_h_is + b_h_is);	
	
	a_n_is = 0.02*(10-x[5])/(exp((10-x[5])/10)-1);
	b_n_is = 0.25*exp(-x[5]/80);
	x[12] = a_n_is/(a_n_is+b_n_is);
	
   	//integrate system form 0 ms to 50 ms with Runge-Kutta methods with time step 0.01 
	integrate( traub , x , 0.0 , 50.0 , 0.01 , write_traub );
}
