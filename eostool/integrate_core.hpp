#ifndef INTEGRATE_CORE_HPP
#define INTEGRATE_CORE_HPP


#include<boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include<boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#include<boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include<boost/numeric/odeint/integrate/integrate_const.hpp>
#include<iomanip>
#include"eos.hpp"


using namespace std;
using namespace boost::numeric::odeint;
typedef std::vector< double > state_type;


//----------------------------------------------------------------------------------------------
//basic integrate controller initiate
//----------------------------------------------------------------------------------------------


typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
controlled_stepper_type controlled_stepper(default_error_checker< double, range_algebra, default_operations >( int_abs_err, int_rel_err, int_a_x, int_a_dxdt ) );//controlled_stepper_type
state_type integrate_func(state_type x0, void func(const state_type &, state_type &, double), double start_t, double end_t, bool reverse = true);


//----------------------------------------------------------------------------------------------
// Lambda related functions
//----------------------------------------------------------------------------------------------


inline double cal_lambda(double C, double Y){return 16*pow(1-2*C, 2.0)*(2+2*C*(Y-1)-Y)/(15*
            (4*pow(C, 3.0)*(13-11*Y+C*(3*Y-2)+2*pow(C, 2.0)*(1+Y))+3*pow(1-2*C, 2.0)*(2-Y+2*C*(Y-1))*log(1-2*C)+2*C*(6-3*Y+3*C*(5*Y-8))));}
inline double lambda_tilde(double m1, double m2, double l1, double l2){return 16.0/13.0*(((12.0*m2+m1)/pow((m1+m2), 5))*pow(m1, 4)*l1+((12.0*m1+m2)/pow((m1+m2), 5))*pow(m2, 4)*l2);}


//----------------------------------------------------------------------------------------------
//ODE group to be solved, x = (p, e, m, r, y)
//----------------------------------------------------------------------------------------------


inline double dpdh(double p, double e){return e+p;}

inline double dedh(double p, double e, double Gamma){return pow((e+p), 2.0)/(p*Gamma);}

inline double dmdh(double p, double e, double m, double r){return -(4*M_PI*e*pow(r, 3.0)*(r-2*m))/(m+4*M_PI*pow(r, 3.0)*p);}

inline double drdh(double p, double e, double m, double r){return -(r*(r-2*m))/(m+4*M_PI*pow(r, 3.0)*p);}

double dydh(double p, double e, double m, double r, double y, double Gamma){
    double dydh1 = ((r-2*m)*(y+1)*y)/(m+4*M_PI*pow(r, 3.0)*p)+y+((m-4*M_PI*pow(r, 3.0)*e)*y+4*M_PI*pow(r, 3.0)*(5*e+9*p)-6*r)/(m+4*M_PI*pow(r, 3.0)*p);
    double dydh2 = (4*M_PI*pow(r, 3.0)*pow(e+p, 2))/((m+4*M_PI*pow(r, 3.0)*p)*p*Gamma)-4*(m+4*M_PI*pow(r, 3.0)*p)/(r-2*m); 
    return (dydh1+dydh2);
}
//inline double dbmdh(double p, double rho, double m, double r){return -(4*M_PI*rho*pow(r, 3.0)*(r-2*m))/(m+4*M_PI*pow(r, 3.0)*p);}
//inline double dIdh(double p, double rho, double m, double r){return -(8./3.*M_PI*rho*pow(r, 5.0)*(r-2*m))/(m+4*M_PI*pow(r, 3.0)*p);}


//----------------------------------------------------------------------------------------------
// EOS integrate function
//----------------------------------------------------------------------------------------------


// Calculate p, e, rho
void cal_eos(const state_type &x, state_type &dxdt, double h){
    double gamma = EOS->gamma(h, x[0], x[1]);
    dxdt[0] = dpdh(x[0], x[1]);
    dxdt[1] = dedh(x[0], x[1], gamma);
    dxdt[2] = dxdt[1]/exp(h); // previous expression: x[2]/(gamma*x[0])*dxdt[0]
    //cout<<h<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<dxdt[0]<<" "<<dxdt[1]<<" "<<dxdt[2]<<" "<<gamma<<endl;
}


//----------------------------------------------------------------------------------------------
//second layer utilities to integrate
//----------------------------------------------------------------------------------------------


void mrl_integration_system(const state_type &x, state_type &dxdt, double h){
    double p_h = EOS->ph(h), e_h = EOS->eh(h);
    double gamma = EOS->gamma(h, p_h, e_h);
    dxdt[0] = dmdh(p_h, e_h, x[0], x[1]);
    dxdt[1] = drdh(p_h, e_h, x[0], x[1]);
    dxdt[2] = dydh(p_h, e_h, x[0], x[1], x[2], gamma);
    //cout<<h<<"   "<<x[0]<<"   "<<x[1]*r_trans<<"   "<<x[2]<<"   "<<1./dxdt[1]<<"   "<<1./dxdt[1]/h<<endl;
}


//----------------------------------------------------------------------------------------------
//third layer utilities to integrate the whole star
//----------------------------------------------------------------------------------------------


/** @brief Integrated portable function to integrate ODE. */
state_type integrate_func(state_type x0, void func(const state_type &, state_type &, double), double start_t, double end_t, bool reverse){
    vector<state_type> x_o;//x for output
    state_type h_o;//h for output
    x_o.reserve(300); h_o.reserve(300);
    double dt;
    //lower bound must smaller than upper bound, so special care should be taken
    try{
        if (consid_const_inter_step){
            if (reverse) dt = -sg_const_step; else dt = sg_const_step;
            integrate_const(error_stepper_type(), func, x0, start_t, end_t, dt, push_back_state_and_time(x_o, h_o));
        }
        else{
            if (reverse) dt = -sg_step; else dt = sg_step;
            integrate_adaptive(controlled_stepper, func, x0, start_t, end_t, dt, push_back_state_and_time(x_o, h_o));
        }
        //calculate neutron star structure and store in an interpolation function list
        if (cal_internal_structure){
            state_type internal_structure_h, internal_structure_m, internal_structure_r;
            structure_function_h_base.clear();
            int intg_steps = h_o.size();
            for (int i=1; i<=intg_steps; i++){
                internal_structure_h.push_back(h_o[intg_steps-i]);
                internal_structure_m.push_back(x_o[intg_steps-i][0]);
                internal_structure_r.push_back(x_o[intg_steps-i][1]*r_trans);
            }
            state_type h_bk1(internal_structure_h), h_bk2(internal_structure_h);
            auto function_h_m = pchip(std::move(h_bk1), std::move(internal_structure_m));
            auto function_h_r = pchip(std::move(h_bk2), std::move(internal_structure_r));
            structure_function_h_base.push_back(function_h_m);
            structure_function_h_base.push_back(function_h_r);
        }
    }
    catch (exception & except){
        cout<<"Error: "<<except.what()<<endl;
        throw;
    }
    //debug
    if (vverbose){
        cout<<"iter steps:"<<h_o.size()<<endl<<endl;
        cout<<"h_i"<<"\t\t"<<"m"<<"\t\t"<<"r"<<"\t\t"<<"y"<<endl;
        cout<<std::setprecision(10);
        for(int j = 0; j<h_o.size(); j++){
            cout<<h_o[j]<<"\t\t"; 
            for (int i = 0; i<x_o[j].size(); i++) cout<<x_o[j][i]<<"\t\t";
            cout<<endl;
        }
    }
    return x_o.back();
}


/** @brief Integrate the most central core part of the compact star to avoid singular point. */
state_type initiate_core(double hc, double hig=h_ig){
    double ec, pc, Gamma_c;
    state_type init_x(3);
    try{
        pc = EOS->ph(hc), ec = EOS->eh(hc);
        Gamma_c = EOS->gamma(hc, pc, ec);
    }
    catch (exception & except){
        if (verbose) cout<<"Invalid value encountered in initiate_core: "<<except.what()<<endl;
        throw;
    }
    double r1 = pow(3.0/(2*M_PI*(ec+3*pc)), 0.5);
    double r3 = -r1/(4*(ec+3*pc))*(ec-3*pc-3*pow(ec+pc, 2.0)/(5*pc*Gamma_c));
    double m3 = 4*M_PI/3.0*ec*pow(r1, 3.0);
    double m5 = 4*M_PI*pow(r1, 3.0)*(r3*ec/r1-pow(ec+pc, 2.0)/(5*pc*Gamma_c));
    double y2 = -6*(ec/3+11*pc+pow(ec+pc, 2.0)/(pc*Gamma_c))/(7*(ec+3*pc));
    init_x[0] = m3*pow(hig, 1.5)+m5*pow(hig, 2.5);
    init_x[1] = r1*pow(hig, 0.5)+r3*pow(hig, 1.5);
    init_x[2] = 2+y2*(hig);
    if (vverbose) cout<<"M_init:"<<init_x[0]<<"    R_init:"<<init_x[1]<<endl;
    return init_x;
}


/** @brief Integrate the whole structure of the star.
  * @details Give h_c, get m, r and y.
  * @note Radius r in this function should multiply half Schwarzschild radius of the sun.
  * @param hc enthalpy at the center of the compact star.
*/
state_type int_whole_star(double hc, double h_surf=h_surface){
    state_type  x_result(3);
    double hc_lowerbound = h_0+0.01;
    if ((param_method!=MITBAG_QUARK_STAR) and (hc<hc_lowerbound)){
        cout<<"Invalid value encountered in int_whole_star, hc="<<hc<<" < "<<hc_lowerbound<<" not valid"<<endl;
        throw exception();
    }
    try{
        state_type init_x = initiate_core(hc);
        if (param_method==MITBAG_QUARK_STAR){//deal with the energy jump
            double h_surf_qks = EOS->params_eos[3];
            double p_surf_qks = EOS->ph(h_surf_qks);
            if (p_surf_qks<1e-300){ //avoid divided by zero in calculating lambda
                EOS->params_eos[3] += 1e-8;
                h_surf_qks = EOS->params_eos[3];
                p_surf_qks = EOS->ph(h_surf_qks);
            }
            double e_surf_qks = EOS->eh(h_surf_qks);
            if (h_surf!=h_surface) { h_surf_qks=h_surf; cout<<"Don't believe in the tidal parameters now!"<<endl; }
            x_result = integrate_func(init_x, mrl_integration_system, hc-h_ig, h_surf_qks);
            double average_e = x_result[0]/(4.*M_PI*pow(x_result[1],3)/3.);
            x_result[2] -= e_surf_qks/(average_e/3.+p_surf_qks); //see arXiv: 1004.5098v1(Eq.14), and the correction in 2007.01139v2(Eq.11)
            //cout<<e_surf_qks<<"  "<<average_e/3.<<"   "<<p_surf_qks<<endl;
        }
        else {
            x_result = integrate_func(init_x, mrl_integration_system, hc-h_ig, h_surf);
        }
    }
    catch (exception &){
        if (verbose) cout<<"Invalid value encountered in int_whole_star("<<hc<<"), set results to 0"<<endl;
        x_result[0]=0;x_result[1]=0;x_result[2]=0;
        throw;
    }
    return x_result;
}


#endif


//debuger
// cout<<"let's go"<<endl;
// cout<<"I am touched!"<<endl;
// cout<<"bye"<<endl;
// cout<<"I am far from ok!"<<endl;
// cout<<"I am pretty far from ok!"<<endl;
// cout<<"I am pretty fucking far from ok!"<<endl;
