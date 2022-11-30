#ifndef GLOBAL_VARIABLE_CONSTANTS_HPP
#define GLOBAL_VARIABLE_CONSTANTS_HPP


#include<cmath>
#include<gsl/gsl_const_cgs.h>
#include<boost/math/interpolators/pchip.hpp>


using namespace std;
using boost::math::interpolators::pchip;
typedef std::vector< double > state_type;



//----------------------------------------------------------------------------------------------
//                                            Global Variables 
//----------------------------------------------------------------------------------------------


//=======parameterization method(method of joint low density table and high density)=======


class EoS;
int param_method = 0;
EoS *EOS;


//=======Control and storage related=======

//verbose level
bool verbose = 0;///< (command: -v).
bool vverbose = 0;///< (command: -vv).
bool vvverbose = 0;///< (command: -vvv).
//tov mass allowed
double minm_tov = 0.0;///<minimum limit of TOV mass
double maxm_tov = 0.0;///<maximum limit of TOV mass
bool check_causal = 1;
//efault border between high density and low density
double h_0 = 2.161113031691454245e-02;///< border of two eos approximation methods.
//important global volatile variables
double mrl_result[3], eos_props[6];


//=======internal structure related
bool cal_internal_structure=false;///< whether you want to find the internal structure of the neutron star given central density
vector< pchip<state_type> > structure_function_h_base;///< Store the internal structure as interpolated functions, reachable only if cal_internal_structure=true
double structure_mr[2];



//=======important variables for integration control=======

bool consid_const_inter_step = false;///< constant integration step
double sg_const_step = 0.;///< single step in constant integration step



//----------------------------------------------------------------------------------------------
//                                         Global Constants 
//----------------------------------------------------------------------------------------------


//=======unit transformation related constants=======

//meta unit
const double G = GSL_CONST_CGS_GRAVITATIONAL_CONSTANT;
const double C = GSL_CONST_CGS_SPEED_OF_LIGHT;
const double Hb = GSL_CONST_CGS_PLANCKS_CONSTANT_HBAR;
const double eV = GSL_CONST_CGS_ELECTRON_VOLT;
const double Ms = GSL_CONST_CGS_SOLAR_MASS;
//natural to cgs
const double MeV_to_ifm = eV/(Hb*C*1.e7);//original: 1e6*eV/(Hb*C*1e13)
const double MeV3_to_ifm3 = pow(MeV_to_ifm, 3);
const double MeVEifm3_to_ergEicm3 = 1.e45*eV;//original: 1e6*eV/(1e-13)^3
const double MeVEifm3_to_gEicm3 = 1.e45*eV/pow(C, 2);//original: MeVEifm3_to_ergEicm3/C^2
const double m_neutron = 939.56542052; // mass of neutron in MeV
//cgs to cactus
double gEicm3_to_dmls = pow(G, 3)*pow(Ms, 2)/pow(C, 6);//original: 1/Ms*(G*Ms/C^2)^3
double ergEicm3_to_dmls = pow(G, 3)*pow(Ms, 2)/pow(C, 8);//original: gEicm3_to_dmls/C^2
//natural to cactus
double MeV4_to_dmls = MeV3_to_ifm3*MeVEifm3_to_ergEicm3*ergEicm3_to_dmls;
double MeVEifm3_to_dmls = MeVEifm3_to_ergEicm3*ergEicm3_to_dmls;
//constants (transform factor with CGS unit)
double length_trans = G*Ms/pow(C, 2);
double rho_trans = pow(C, 6)/(pow(G, 3)*pow(Ms, 2));
double p_trans = rho_trans*pow(C, 2);
double r_trans = length_trans/1.e5;
double rho_sat = 2.7e+14/rho_trans;

//important constants for integration control
const double int_abs_err = 1.0e-18, int_rel_err = 1.0e-12, int_a_x = 1.0, int_a_dxdt = 1.0;/// integrate constants
const double sg_step = 1.e-10;///< single step.
const double h_surface = 1e-10;///< surface of the star.
const double h_ig=1e-7;///< where to start integrate TOV equation.
//border parameters between high density and low density
double e_0 = 9.075720079516068750e+13/rho_trans;
double p_0 = 2.979148040306152863e+32/p_trans;
double rho_0 = (e_0+p_0)/exp(h_0);


//=======struct/enum used frequently=======


enum EoS_type {
    TABULATED_ONLY = 1, ///< interpolation only.
    PIECEWISE_GAMMA = 2, ///< piecewise polytropic expansion, conflict with causal.
    PIECEWISE_PRESSURE = 3, ///< piecewise polytropic expansion with pressure at different densities. (see Jiang20)
    SPECTRAL_HBASE = 4, ///< h-based spectral expansion and interpolation.
    SPECTRAL_HBASE_CAUSAL = 5, ///< note that Gamma will have different meaning, see Lindblom_18.
    PWG_PHASE_TRANS = 6, ///< piecewise gamma eos+phase transition(https://arxiv.org/pdf/1811.10929.pdf, Model2).
    MITBAG_QUARK_STAR = 7, ///< quark star model(Alford_2005_ApJ_629_969), this model will not have atmosphere implemented.
    ADAPT_PIECEWISE = 8, ///< adaptive rho border, with parameters: delta_rho*n+speed_of_sound*n.
    CONS_CS = 9, ///< constant speed of sound.
    PIECE_SPEC_PHTR_CSS = 10, ///<piecewise+spectral+phase transition+constant speed of sound (see Tang20)
    PHYSICAL_SPEC_PHTR_CSS = 11, ///<physical+spectral+phase transition+constant speed of sound (see Tang21)
    PARAM_MU_CS = 12, ///< parameterize the speed of sound, see https://doi.org/10.1038/s41567-020-0914-9
    PARAM_MU_CS_PT = 13 ///< parameterize the speed of sound with phase transition, see ?
};


struct push_back_state_and_time {
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;
    push_back_state_and_time( std::vector< state_type > &states, std::vector< double > &times ): m_states( states ), m_times( times ) { }
    void operator()( const state_type &x, double t ){
        if (m_states.size()%20==0){
            if (vvverbose) cout<<"        t: "<<t<<", x: ";
            for(int i=0; i<x.size();i++){
                if (vvverbose) cout<<x[i]<<"    ";
                if (x[i]>1e6 or std::isnan(x[i]) or m_times.size()>1e5) {
                    if (verbose) cout<<"Invalid value x["<<i<<"]="<<x[i]<<" encountered at t="<<t<<", with x size: "<<x.size()<<", iteration times: "<<m_times.size()<<endl;
                    throw exception();
                }
            }
            if (vvverbose) cout<<", x size: "<<x.size()<<endl;
        }
        m_states.push_back( x ); m_times.push_back( t );
    }
};


#endif
