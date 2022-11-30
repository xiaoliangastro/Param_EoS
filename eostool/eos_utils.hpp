#ifndef EOS_UTILS_HPP
#define EOS_UTILS_HPP


#include<iostream>
#include<exception>
#include<fstream>
#include"integrate_core.hpp"
#include<boost/math/interpolators/makima.hpp>


using namespace std;
using boost::math::interpolators::makima;
typedef std::vector< double > state_type;


struct control_params{
    double min_tov_mass; ///< Minimum mass allowed astropysically, used to contrain the EoS.
    double max_tov_mass; ///< Maximum mass allowed astropysically, used to contrain the EoS.
    bool check_causal; ///< Some EoS do not causal in the maximum mass place, then find the hc where cs=1, the maximum mass in this case is defined by M(hc|_{cs=1}).
    int param_method; ///< Parameterization method.
    int verbose_level; ///< Verbose level, 0-3, more details will be printed if verbose_level is larger.
    double const_inter_step; ///< If consider constant integration step, then please give the wanted step size.
    bool consid_const_inter_step; ///< Consider constant integration step or adaptive step size.
    bool cal_internal_structure; ///< whether you want to find the internal structure of the neutron star given central density
};


extern "C"{
    double eos_messenger(int ask_for_information_type);
    void init_control_params(control_params parameters);
    void create_low_density_eos_with_ep_table(double ee[], double pp[], const int len_eos);
    bool update_high_density_eos_parameters(double ppar[], int len_par, int len_border, double h_max, const char init_function_type[]);
    double* find_eos_properties(double known_aim, int find_type);
    void make_eos_table(double h_start, double h_end, double dh, int precision, char *f_name, const char *unit="cgs", int get_type=1);
    double interp_pe_likelihood(double *e_org, double *p_org, int size_org);
}


//----------------------------------------------------------------------------------------------
//                                    EoS calculators
//----------------------------------------------------------------------------------------------


double eos_messenger(int ask_for_information_type){
    return EOS->eos_messenger(ask_for_information_type);
}


bool integrate_eos(double h_start, double h_end, double dh, state_type *h_tb, state_type *p_tb, state_type *e_tb, \
                   state_type *rho_tb, state_type *gamma_tb, state_type *v_sq_tb, int get_type){
    bool success = true;
    double h = h_end, p, e, gamma, rho, v_square;
    if ((param_method==TABULATED_ONLY) or (param_method==MITBAG_QUARK_STAR)){
        for (int i = 0; i<floor((h_end-h_start)/dh); i++){
            p = EOS->ph(h), e = EOS->eh(h), rho = (p+e)/exp(h);
            gamma = EOS->gamma(h, p, e);
            v_square = (p*gamma)/(e+p);
            h_tb->push_back(h); p_tb->push_back(p); e_tb->push_back(e);
            rho_tb->push_back(rho); gamma_tb->push_back(gamma); v_sq_tb->push_back(v_square); 
            h -= dh;
        }
    }
    else {
        vector <state_type> per_o;
        state_type per_c, h_o;
        state_type x0(3); x0[0] = p_0; x0[1] = e_0; x0[2] = rho_0;
        try{
            per_c = integrate_func(x0, cal_eos, h_0, h_end, false);
            for (int i = 0; i<floor((h_end-h_start)/dh); i++){
                if (h>h_0 and get_type!=1){
                    integrate_adaptive(controlled_stepper, cal_eos, per_c, h, h-dh, -1e-6, push_back_state_and_time(per_o, h_o));
                    p = per_o.back()[0], e = per_o.back()[1], rho = per_o.back()[2];
                }
                else {
                    p = EOS->ph(h), e = EOS->eh(h), rho = (p+e)/exp(h);
                }
                gamma = EOS->gamma(h, p, e);
                v_square = (p*gamma)/(e+p);
                h_tb->push_back(h); p_tb->push_back(p); e_tb->push_back(e);
                rho_tb->push_back(rho); gamma_tb->push_back(gamma); v_sq_tb->push_back(v_square);
                h -= dh;
            }
        }
        catch (exception & except){
            cout<<"failed in integrate_eos, "<<except.what()<<endl;
            success = false;
        }
    }
    return success;
}


/** @brief Python interface function for finding eos properties from h(find_type=1), 
           e(find_type=2) or rho(find_type=3).
  * @params find_type must equal to init_function_type-1.
*/
double* find_eos_properties(double known_aim, int find_type){
    // h, p, e, rho, gamma, cssq; please also update eos with the same find_type
    double h, p, e, rho;
    #define fail_process_find_eos_prop() for (int i=0; i<6; i++) {eos_props[i] = 0.;} return eos_props;
    try{
        if (find_type==1){
            h = known_aim;
            p = EOS->ph(h);
            e = EOS->eh(h);
            rho = EOS->rhoh(h, p, e);
        }
        else if (find_type==2){
            e = known_aim;
            h = EOS->eos_table_function_e_base[0](e);
            p = EOS->eos_table_function_e_base[1](e);
            rho = EOS->eos_table_function_e_base[2](e);
        }
        else if (find_type==3){
            rho = known_aim;
            h = EOS->eos_table_function_rho_base[0](rho);
            p = EOS->eos_table_function_rho_base[1](rho);
            e = EOS->eos_table_function_rho_base[2](rho);
        }
        else {
            cout<<"Unknown find_type in find_eos_properties, return zero"<<endl;
            fail_process_find_eos_prop();
        }
    }
    catch (exception & except){
        cout<<"failed in find_eos_properties: "<<except.what()<<endl;
        fail_process_find_eos_prop();
    }
    double gamma = EOS->gamma(h, p, e), v_square = (p*gamma)/(e+p);
    eos_props[0] = h, eos_props[1] = p, eos_props[2] = e;
    eos_props[3] = rho, eos_props[4] = gamma, eos_props[5] = v_square;
    return eos_props;
}


void make_eos_table(double h_start, double h_end, double dh, int precision, char *f_name, const char *unit, int get_type){
//check whether e, p have been correctly interpolated, cgs: p, e, rho
    state_type h_tb, p_tb, e_tb, rho_tb, gamma_tb, v_sq_tb;
    if (verbose) cout<<"Output EoS table to file: "<<f_name<<endl;
    ofstream fp(f_name);
    fp.precision(precision);
    integrate_eos(h_start, h_end, dh, &h_tb, &p_tb, &e_tb, &rho_tb, &gamma_tb, &v_sq_tb, get_type);
    if (strcmp(unit, "cgs")==0) fp<<"h"<<"\t\t"<<"p/dyn*cm-2"<<"\t\t"<<"e/erg*cm-3"<<"\t\t"<<"(be careful)rho/g*cm-3"<<"\t\t"<<"gamma"<<"\t\t"<<"v_square"<<endl;
    else fp<<"h"<<"\t\t\t"<<"p"<<"\t\t\t"<<"e"<<"\t\t\t"<<"rho"<<"\t\t\t"<<"gamma"<<"\t\t\t"<<"v_square"<<endl;
    for (int i=0; i<h_tb.size(); i++){
        if (strcmp(unit, "cgs")==0) {fp<<h_tb[i]<<"\t\t"<<p_tb[i]*p_trans<<"\t\t"<<e_tb[i]*rho_trans<<"\t\t"<<rho_tb[i]*rho_trans<<"\t\t"<<gamma_tb[i]<<"\t\t"<<v_sq_tb[i]<<endl;}
        else fp<<h_tb[i]<<"\t\t"<<p_tb[i]<<"\t\t"<<e_tb[i]<<"\t\t"<<rho_tb[i]<<"\t\t"<<gamma_tb[i]<<"\t\t"<<v_sq_tb[i]<<endl;
    }
    fp.close();
}


//----------------------------------------------------------------------------------------------
//                                    EoS utils
//----------------------------------------------------------------------------------------------


bool cool_eos(double *hc){
//check whether gamma(h)<10
    bool valid = true;
    double h_lowest = EOS->suggested_hc_star;
    double h = *hc, dh=0.001, gamma, p, e, v_square;
    try{
        for (int i = 0; i<floor((*hc-h_lowest)/dh); i++){
            p = EOS->ph(h), e = EOS->eh(h); gamma = EOS->gamma(h, p, e);
            v_square = (p*gamma)/(e+p);
            valid = (v_square<1.) and (gamma<10); // must be this way, double checked
            if (verbose) cout<<"h: "<<h<<"   p: "<<p<<"   e: "<<e<<"    gamma: "<<gamma<<"   v^2: "<<v_square<<endl;
            if (valid) break;
            h -= dh;
        }
    }
    catch (exception & except){
        cout<<"Invalid value encountered in cool eos: "<<except.what()<<endl;
        valid = false;
    }
    if (valid) *hc = h;
    return valid;
}


/** @brief calculate the segement point of the energy denisty in the piecewise approx (see Appendix B of 10.1103/PhysRevD.86.084003) */
state_type cal_eborders(state_type params, state_type rho_borders){
    state_type e_borders, gamma, pb, rb;
    if ((param_method==PIECEWISE_PRESSURE) or (param_method==ADAPT_PIECEWISE) or (param_method==PIECE_SPEC_PHTR_CSS) or (param_method==PHYSICAL_SPEC_PHTR_CSS)){
        for(int i=0; i<params.size(); i++){
            if (i==0) gamma.push_back(log(params[0]/p_0)/log(rho_borders[0]/rho_0));
            else gamma.push_back(log(params[i]/params[i-1])/log(rho_borders[i]/rho_borders[i-1]));
        }
    }
    else if (param_method==PIECEWISE_GAMMA) {
        for(int i=0; i<params.size(); i++) gamma.push_back(params[i]);
    }
    else { cout<<"Ask yourself: Why do you want to calculate eborders?!"<<endl; exit(1); }
    int n_pieces = gamma.size();
    double ai_plus1;
    e_borders.push_back(e_0);
    rb.push_back(rho_0);
    pb.push_back(p_0);
    for (int i=0; i<n_pieces; i++){
        rb.push_back(rho_borders[i]);
        pb.push_back(pb[i]*pow(rb[i+1]/rb[i],gamma[i]));
        if (gamma[i]==1){
            ai_plus1 = e_borders[i]/rb[i]-pb[i]*log(rb[i])/rb[i];
            e_borders.push_back(rb[i+1]*ai_plus1+pb[i+1]*log(rb[i+1]));
        }
        else {
            ai_plus1 = e_borders[i]/rb[i]-pb[i]/(rb[i]*(gamma[i]-1));
            e_borders.push_back(rb[i+1]*ai_plus1+pb[i+1]/(gamma[i]-1));
        }
    }
    if(verbose){
        cout<<"gamma: ";
        for(int i=0; i<gamma.size(); i++) cout<<gamma[i]<<"  ";
        cout<<endl;
        cout<<"e_borders: ";
        for(int i=0; i<e_borders.size(); i++) cout<<e_borders[i]<<"  ";
        cout<<endl;
    }
    return e_borders;
}


//----------------------------------------------------------------------------------------------
//                                  Likelihood of EoS fitting
//----------------------------------------------------------------------------------------------


double interp_pe_likelihood(double *e_org, double *p_org, int size_org){
    double likelihood = 0., increase = 0.;
    try {
        for (int i=0; i<size_org; i++) {
            increase = pow(log(p_org[i]/EOS->eos_table_function_e_base[1](e_org[i])), 2.);
            if (not std::isnan(increase)) likelihood += increase;
        }
    }
    catch (exception & except) {
        cout<<"interpolate error: "<<except.what()<<endl;
        likelihood = +1.e100;
    }
    return likelihood;
}

double interp_pe_likelihood_precise(double *e_org, double *p_org, int size_org, double h_start, double h_end, double dh){
    double likelihood = 0., increase = 0.;
    state_type h_tb, p_tb, e_tb, rho_tb, gamma_tb, v_sq_tb;
    state_type e_buff, p_buff, eorg_buff;
    for (int j=0; j<size_org; j++) eorg_buff.push_back(e_org[j]);
    bool int_eos_success = integrate_eos(h_start, h_end, dh, &h_tb, &p_tb, &e_tb, &rho_tb, &gamma_tb, &v_sq_tb, 2);
    if (not int_eos_success) {
        cout<<"integrate eos error at interp_pe_likelihood."<<endl;
        return +1.e100;
    }
    state_type::iterator e_tb_ub = std::upper_bound(e_tb.begin(), e_tb.end(), e_org[0], [](double val, double element){ return val>element;});
    if (e_tb_ub != e_tb.end()) e_tb_ub += 1;
    int distl = std::distance(e_tb.begin(), e_tb_ub);
    std::reverse_copy(e_tb.begin(), e_tb_ub, std::back_inserter(e_buff));
    std::reverse_copy(p_tb.begin(), p_tb.begin()+distl, std::back_inserter(p_buff));
    state_type::iterator e_tb_lb = std::lower_bound(eorg_buff.begin(), eorg_buff.end(), e_buff[e_buff.size()-1]);
    int distu = std::distance(eorg_buff.begin(), e_tb_lb);
    try {
        auto f_interp = makima(std::move(e_buff), std::move(p_buff)); //'move' will clear the table
        for (int i=0; i<distu; i++) { 
            increase = pow(log(p_org[i]/f_interp(e_org[i])), 2.);
            if (not std::isnan(increase)) likelihood += increase;
        }
    }
    catch (exception & except) {
        cout<<"interpolate error: "<<except.what()<<endl;
        likelihood = +1.e100;
    }
    return likelihood;
}


//----------------------------------------------------------------------------------------------
//                                    Python control utilities
//----------------------------------------------------------------------------------------------


/** @brief initiate global control parameters */
void init_control_params(control_params pars){
    int verbose_level = pars.verbose_level;
    param_method = pars.param_method;
    minm_tov=pars.min_tov_mass, maxm_tov = pars.max_tov_mass, check_causal=pars.check_causal;// maximum mass constraint
    sg_const_step = pars.const_inter_step; consid_const_inter_step = pars.consid_const_inter_step;// constant integrate step
    cal_internal_structure = pars.cal_internal_structure;
    //cout<<"TOV: "<<minm_tov<<"  "<<maxm_tov<<" "<<check_causal<<" "<<param_method<<" "<<verbose_level<<" "<<sg_const_step;
    // verbose level
    if (verbose_level==0) {}
    else if (verbose_level==1) {verbose=1; vverbose=0; vvverbose=0;}
    else if (verbose_level==2) {verbose=1; vverbose=1; vvverbose=0;}
    else if (verbose_level==3){verbose=1; vverbose=1; vvverbose=1;}
    else {cout<<"Unknown verbose level: "<<verbose_level<<endl; exit(1);}
}


void create_low_density_eos_with_ep_table(double ee[], double pp[], const int len_eos){
    //parameterization method
    switch (param_method){
        case TABULATED_ONLY:
            {EOS = new EoS_tabulated_only(ee, pp, len_eos); break;}
        case PIECEWISE_GAMMA:
            {EOS = new EoS_piecewise_gamma(ee, pp, len_eos); break;}
        case PIECEWISE_PRESSURE:
            {EOS = new EoS_piecewise_pressure(ee, pp, len_eos); break;}
        case SPECTRAL_HBASE:
            {EOS = new EoS_spectral_hbase(ee, pp, len_eos, false); break;}
        case SPECTRAL_HBASE_CAUSAL:
            {EOS = new EoS_spectral_hbase(ee, pp, len_eos, true); break;}
        case PWG_PHASE_TRANS:
            {EOS = new EoS_pwg_pt(ee, pp, len_eos); break;}
        case MITBAG_QUARK_STAR:
            {EOS = new EoS_MIT_bag(); break;}
        case ADAPT_PIECEWISE:
            {EOS = new EoS_adp_pw(ee, pp, len_eos); break;}
        case CONS_CS:
            {EOS = new EoS_cons_cs(ee, pp, len_eos); break;}
        case PIECE_SPEC_PHTR_CSS:
            {EOS = new EoS_pw_sp_pt_css(ee, pp, len_eos); break;}
        case PHYSICAL_SPEC_PHTR_CSS:
            {EOS = new EoS_ph_sp_pt_css(ee, pp, len_eos); break;}
        case PARAM_MU_CS:
            {EOS = new EoS_param_mu_cs(ee, pp, len_eos); break;}
        case PARAM_MU_CS_PT:
            {EOS = new EoS_param_mu_cs_PT(ee, pp, len_eos); break;}
        default:
            {cout<<"Unknown parameterization method: "<<param_method<<endl; exit(1);}
    }
}


/** @brief accept important parameters from python script */
bool update_high_density_eos_parameters(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    bool ret_val = true;
    try{
        ret_val = EOS->update_eos(ppar, len_par, extra_par, h_max, init_function_type);
    }
    catch (exception &except){
        cout<<"Failed to change parameters in change_pars, "<<except.what()<<endl;
        ret_val = false;
    }
    return ret_val;
}


#endif