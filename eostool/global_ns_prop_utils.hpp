#ifndef PARAM_UTILS_HPP
#define PARAM_UTILS_HPP


#include"eos_utils.hpp"
#include<cstdlib>


// python interfases
extern "C"{
    double* get_mrl(double hc);
    double* calculate_internal_structure(double h);
    double* get_mr_with_specific_hsurf(double hc, double h_surf);
    bool check_mmax_gd(double *hc, double *M_max, double h_start=EOS->suggested_hc_star+0.08, bool compare_ftype_is_less=true, double jump_size=0.1, bool check_ok=true);
    bool check_mmax_pt_two_branch(double* h1, double* h2, double* h3, double* h_max, double* m1, double* m2, double* m3, double* m_max, double rho_tr, double drho);
    bool find_closest_global_property(double m_aim, double h_i, double h_max, double *h_closest, double *lambda, bool use_user_start_point=false);
    bool find_closest_global_property_with_maxm_known(double known_aim, double h_i, double h_max, double *h_closest, double *unknown, int get_type, bool use_user_start_point=false);
    bool get_unknowns_from_knowns(double known1, double known2, double *unknown1, double *unknown2, double *h_max, int get_type);
}


//----------------------------------------------------------------------------------------------
//       Mass-Radius-Lambda Tools
//----------------------------------------------------------------------------------------------


/** @brief Python interface function for integrating the whole structure of the star.
* @details Give h_c, get M/M_sun, R/km and L, which are stored in global variable double *mrl_result.
* @attention Use this for python control only, because double* is convenient for python, and global variable is needed,  
          which may cause confuse problem if used in c++ main, so please use int_whole_star instead for c++, 
          because that do not need to return global variable.
* @param hc enthalpy at the center of the compact star.
* @retval double* A double array of [mass, radius, tidal-deformability].
*/
double* get_mrl(double hc){
    state_type x_result(3);
    try {
        x_result = int_whole_star(hc);
        mrl_result[0] = x_result[0], mrl_result[1] = x_result[1]*r_trans;
        mrl_result[2] = cal_lambda(x_result[0]/x_result[1], x_result[2]);
    }
    catch (exception &){
        if (verbose) cout<<"Invalid value encountered in get_mrl("<<hc<<"), set results to 0"<<endl;
        mrl_result[0] = 0, mrl_result[1] = 0, mrl_result[2] = 0;
    }
    return mrl_result;
}

// get the m-r structure of neutron star: do the whole TOV integration again and again
double* get_mr_with_specific_hsurf(double hc, double h_surf){
    state_type x_result(3);
    try {
        x_result = int_whole_star(hc, h_surf);
        structure_mr[0] = x_result[0], structure_mr[1] = x_result[1]*r_trans;
    }
    catch (exception &){
        if (verbose) cout<<"Invalid value encountered in get_mr_with_specific_hsurf("<<hc<<", "<<h_surf<<"), set results to 0"<<endl;
        structure_mr[0] = 0, structure_mr[1] = 0;
    }
    return structure_mr;
}

// get the m-r structure of neutron star: do the whole TOV integration once in get_mrl(hc), then calculate m, r using the interpolated integration steps in get_mrl
double* calculate_internal_structure(double h_surf){
    try {
        structure_mr[0] = structure_function_h_base[0](h_surf);
        structure_mr[1] = structure_function_h_base[1](h_surf);
    }
    catch (exception &){
        if (verbose) cout<<"Invalid value encountered in calculate_internal_structure("<<h_surf<<"), set results to 0"<<endl;
        structure_mr[0] = 0, structure_mr[1] = 0;
    }
    return structure_mr;
}

//make table of [h, m, r, l]
void make_hmrl_tool_table(double minh, double maxh, double dh){
    double h = minh;
    state_type result_x(3);
    string tool_fname = "test_hmrl.txt";
    ofstream ofp(tool_fname);
    for (int i = 0; i<floor((maxh-minh)/dh); i++){
        try{
            result_x = int_whole_star(h);
            double tidal_deform = cal_lambda(result_x[0]/result_x[1], result_x[2]);
            ofp<<h<<"\t\t"<<result_x[0]<<"\t\t"<<result_x[1]*r_trans<<"\t\t"<<tidal_deform<<endl;
            //tidal_deform*pow(result_x[0]*r_trans, 5), tidal_deform, 3./2.*tidal_deform*pow(result_x[0]/result_x[1], 5)
        }
        catch (exception &){
            cout<<"Invalid value encountered in make_tool_table, h="<<h<<endl;
        }
        h += dh;
    }
    ofp.close();
}


//----------------------------------------------------------------------------------------------
//       Find maximum mass related
//----------------------------------------------------------------------------------------------


/**
* @brief Find non-rotating maximum supported mass of an EoS.
* @details
* The challenge is the "noise" or none monotonicity property, the solution is:
* 1. A bigger guard to check again, 2. Judge whether a step before is larger.
* error: ~0.005 M_sun
* @param hc  Pointer to the central enthalpy, return -1 if finding error.
* @param M_max Pointer to the maximum mass, return 0 if finding error.
* @param h_start Where to start finding the maximum mass.
* @param compare_ftype_is_less Find mimimum or maximum, set to true if maximum is wanted.
* @param jump_size Set the initial jump size of h to find the maximum mass.
* @param check_ok Whether to set the (mtov_min, mtov_max) constraint or not.
* @retval bool whether finding process correctly worked.
*/
bool check_mmax_old(double *hc, double *M_max, double h_start, bool compare_ftype_is_less, double jump_size, bool check_ok){
    auto cp_func = [=](double cp1, double cp2){
        if (compare_ftype_is_less) return cp1<cp2;
        else return cp1>cp2;
    };
    bool jump_back=false, force_jump_back=false;
    state_type x_result(3), x_guard(3), x_check_again(3);
    double h_max=EOS->eos_table_max_h, h_guard=4e-3, h_check_again=0.01;
    double h = h_start;
    double local_minimum_allowed_h = h_start-1e-4;
    double M=0., M_guard=0., M_check_again=0., M_temp=0.;
    #define fail_process_check_mmax(message) \
        *hc = -1; *M_max=0.; \
        if (verbose) { \
            cout<<message<<", with h="<<h<<", h_guard="<<h+h_guard<<", h_check_again="<<h+h_check_again; \
            cout<<", M="<<M<<", M_guard="<<M_guard<<", M_check_again="<<M_check_again<<endl;\
        } \
        return false;
    //find m_max
    if(verbose) cout<<endl<<"find maximum mass, with h_min: "<<local_minimum_allowed_h<<", h_max: "<<h_max<<endl;
    while(h>local_minimum_allowed_h and h<h_max and jump_size>=5e-4){
        if (jump_back) {
            if (force_jump_back) h -= jump_size;
            h -= jump_size;
            jump_size /= 5.;
            h += jump_size;
        }
        else h += jump_size;
        try{
            x_result = int_whole_star(h);
            M = x_result[0];
            x_guard = int_whole_star(h+h_guard);
            M_guard = x_guard[0];
        }
        catch (exception &){
            if (verbose) cout<<"Invalid value encountered in check_mmax, with h="<<h<<", or h="<<h+h_guard<<endl;
            M=0; M_guard=0;
        }
        if (verbose) cout<<"h: "<<h<<"     M:"<<M<<"     M_guard:"<<M_guard<<"     "<<endl;
        jump_back = (M==0 or M_guard==0 or std::isnan(M_guard) or std::isnan(M));
        force_jump_back = (M==0 or M_guard==0);
        if ((not jump_back) and (cp_func(M_guard, M) or M_guard==M)){
            try {
                x_check_again = int_whole_star(h+h_check_again);
                M_check_again = x_check_again[0];
                if (verbose) cout<<"M_check_again: "<<M_check_again<<endl;
                if (cp_func(M_check_again, M)) {jump_back = true; if (cp_func(std::max({M, M_guard}, cp_func), M_temp)) force_jump_back = true;} 
            }
            catch (exception &) { if (verbose) cout<<"Invalid value encountered in check_mmax: h="<<h+0.01<<endl;}
        }
        M_temp = std::max({M, M_guard}, cp_func);
    }
    // take the extreme of M and h
    M = std::max({M, M_guard, M_check_again}, cp_func);
    if (M==M_guard) h = h+h_guard;
    else if (M==M_check_again) h = h+h_check_again;
    else {}// M = M, h = h
    if (verbose) cout<<"find h: "<<h<<",    with M:"<<M<<endl<<endl;
    if (std::isnan(M) or M==0 or h<=local_minimum_allowed_h) {fail_process_check_mmax("check_mmax: find max mass error(got M=(NaN, 0) or h too small)");}
    if (check_ok){
        if (M<minm_tov) {if (verbose) cout<<"M_max: "<<M<<" <min_mtov: "<<minm_tov<<endl; fail_process_check_mmax("check_mmax: find max mass error(M<min_m_tov)");}
        if (check_causal){
            if (verbose) cout<<endl<<"cooling: "<<endl;
            if (cool_eos(&h)){
                try {M = int_whole_star(h)[0];}
                catch (exception&){ fail_process_check_mmax("Invalid value encountered in integrate after cooling");}
                if (M>maxm_tov or M<minm_tov) {fail_process_check_mmax("check_mmax: find max mass error(maximum mass invalid after cooling)");}
                if (verbose) cout<<"Max mass after cooling: "<<M<<endl<<endl;
            }
            else{fail_process_check_mmax("check_mmax: find max mass error(cool_eos not success)");}
        }
        else{
            if (M>maxm_tov) {if (verbose) cout<<"M_max: "<<M<<" >max_mtov: "<<maxm_tov<<endl; fail_process_check_mmax("check_mmax: find max mass error without cooling(M>max_m_tov)");}      
        }
    }
    *hc = h; *M_max=M;
    return true;
}


// Find maximum mass using information of gradient, all the parameters have the same meaning with check_mmax.
bool check_mmax_gd(double *hc, double *M_max, double h_start, bool compare_ftype_is_less, double jump_size, bool check_ok){
    auto cp_func = [=](double cp1, double cp2){
        if (compare_ftype_is_less) return cp1<cp2;
        else return cp1>cp2;
    };
    int step_i = 0, sum_sign = 0;
    bool jump_back=false, force_jump_back=false;
    state_type x_result(3), x_guard(3), x_check_again(3);
    double h_max=EOS->eos_table_max_h, h_guard=1e-3, h_check_again=0.01;
    double h = h_start;
    double local_minimum_allowed_h = h_start-1e-10;
    double M=0., M_guard=0., gradient=0.0;
    #define fail_process_check_mmax_gd(message) \
        *hc = -1; *M_max=0.; \
        if (verbose) { \
            cout<<message<<", with h="<<h<<", h_guard="<<h+h_guard; \
            cout<<", M="<<M<<", M_guard="<<M_guard<<endl;\
        } \
        return false;
    //find m_max
    if(verbose) cout<<endl<<"find maximum mass, with h_min: "<<local_minimum_allowed_h<<", h_max: "<<h_max<<endl;
    while(h>local_minimum_allowed_h and h<h_max){
        h += gradient*jump_size;//+rand()%10/1.e4; // avoid dead loop
        try{
            x_result = int_whole_star(h);
            M = x_result[0];
            x_guard = int_whole_star(h+h_guard);
            M_guard = x_guard[0];
            if (compare_ftype_is_less) gradient = (M_guard-M)/h_guard;
            else gradient = (M-M_guard)/h_guard;
        }
        catch (exception &){
            if (verbose) cout<<"Invalid value encountered in check_mmax, with h="<<h<<", or h="<<h+h_guard<<endl;
            M=0; M_guard=0;
        }
        step_i += 1;
        if (verbose) cout<<"step: "<<step_i<<"    h: "<<h<<"     M:"<<M<<"     M_guard:"<<M_guard<<"     gradient:"<<gradient<<"    jump_size:"<<jump_size<<"     dh:"<<gradient*jump_size<<endl;
        if (std::signbit(gradient)) sum_sign -= 1;
        else sum_sign += 1;
        if (abs(gradient*jump_size)<1e-5 or abs(M-M_guard)<1e-5 or step_i>200) break;
        if (step_i%8==0){
            if (sum_sign==0) jump_size /= 2.;
            else if (sum_sign==8) jump_size *= 5.;
            else {}
            sum_sign = 0.;
        }
    }
    // take the extreme of M and h
    M = std::max({M, M_guard}, cp_func);
    if (M==M_guard) h = h+h_guard;
    else {}// M = M, h = h
    if (verbose) cout<<"find h: "<<h<<",    with M:"<<M<<endl<<endl;
    if (std::isnan(M) or M==0 or h<=local_minimum_allowed_h) {fail_process_check_mmax_gd("check_mmax: find max mass error(got M=(NaN, 0) or h too small)");}
    if (check_ok){
        if (M<minm_tov) {fail_process_check_mmax_gd("check_mmax: find max mass error(M<min_m_tov)");}
        if (check_causal){
            if (verbose) cout<<endl<<"cooling: "<<endl;
            if (cool_eos(&h)){
                try {M = int_whole_star(h)[0];}
                catch (exception&){ fail_process_check_mmax_gd("Invalid value encountered in integrate after cooling");}
                if (M>maxm_tov or M<minm_tov) {fail_process_check_mmax_gd("check_mmax: find max mass error(maximum mass invalid after cooling)");}
                if (verbose) cout<<"Max mass after cooling: "<<M<<endl<<endl;
            }
            else{fail_process_check_mmax_gd("check_mmax: find max mass error(cool_eos not success)");}
        }
        else{
            if (M>maxm_tov) {fail_process_check_mmax_gd("check_mmax: find max mass error without cooling(M>max_m_tov)");}      
        }
    }
    *hc = h; *M_max=M;
    return true;
}


// 
/**
* @brief Find maximum mass in the case that there are two branches of M-R (phase transition of EoS).
* @details
* Need the transition density and strenth to find three critical structure
* and output the corresponding hc and mass.
* @param h1  Pointer to the central enthalpy of first peak in M-R graph, return -1 if finding error.
* @param h2  Pointer to the central enthalpy of first valley in M-R graph, return -1 if finding error.
* @param h3  Pointer to the central enthalpy of second peak in M-R graph, return -1 if finding error.
* @param h_max  Pointer to the central enthalpy of maximum mass in M-R graph, return -1 if finding error.
* @param m1 Pointer to the mass corresponding to h1, return 0 if finding error.
* @param m2 Pointer to the mass corresponding to h2, return 0 if finding error.
* @param m3 Pointer to the mass corresponding to h3, return 0 if finding error.
* @param m_max Pointer to the mass corresponding to h_max, return 0 if finding error.
* @param rho_tr The starting mass density of phase transition, to boost the finding process.
* @param drho The length of mass density of phase transition, to boost the finding process.
* @retval bool whether finding process correctly worked.
*/
bool check_mmax_pt_two_branch(double* h1, double* h2, double* h3, double* h_max, double* m1, double* m2, double* m3, double* m_max, double rho_tr, double drho){
    double h_tr_start, h_tr_end, jump_size;
    bool ck_tbranch_l=false, ck_tbranch_u=false;
    *h1 = -1, *h2 = -1, *h3 = -1, *h_max = -1, *m1 = 0., *m2 = 0., *m3 = 0., *m_max=0.;
    #define fail_process_check_mmax_pt_two_branch(message) \
        *h1 = -1, *h2 = -1, *h3 = -1, *h_max = -1, *m1 = 0., *m2 = 0., *m3 = 0., *m_max=0.; \
        if (verbose) cout<<message<<endl; \
        return false;
    try {
        h_tr_start = EOS->eos_table_function_rho_base[0](rho_tr);
        h_tr_end = EOS->eos_table_function_rho_base[0](rho_tr+drho);
    }
    catch (exception&) {h_tr_start=0.; h_tr_end=0.;}
    if (check_mmax_gd(h1, m1, 0.05, true, 0.0003, false)){
        if (verbose) cout<<"h_tr_start="<<h_tr_start<<", h_tr_end="<<h_tr_end<<", h_find="<<*h1<<", delta_h="<<h_tr_start-*h1<<endl;
        if (abs(h_tr_start-*h1)<0.03 or abs(h_tr_end-*h1)<0.03){
            ck_tbranch_l = check_mmax_gd(h2, m2, *h1+0.005, false, 0.002, false);
            if (ck_tbranch_l){
                cout<<"two branches!"<<endl;
                ck_tbranch_u = check_mmax_gd(h3, m3, *h2+0.01, true, 0.2, false);
                if (ck_tbranch_u) {*h_max = *h3; *m_max = *m3;}
                else {fail_process_check_mmax_pt_two_branch("check second branch error");}
            }
        }
        else {
            *h_max = *h1; *m_max = *m1;
        }
        if ((not ck_tbranch_l) or (ck_tbranch_u and *m1>*m3)) {
            *h_max = *h1; *m_max = *m1;
        }
        if (*m_max<minm_tov) {fail_process_check_mmax_pt_two_branch("check_mmax: find max mass error(M<min_m_tov)");}
        if (check_causal){
            if (verbose) cout<<endl<<"cooling: "<<endl;
            if (cool_eos(h_max)){
                try {*m_max = int_whole_star(*h_max)[0];}
                catch (exception&){ fail_process_check_mmax_pt_two_branch("Invalid value encountered in integrate after cooling");}
                if (*m_max>maxm_tov or *m_max<minm_tov) {fail_process_check_mmax_pt_two_branch("check_mmax: find max mass error(maximum mass invalid after cooling)");}
                if (verbose) cout<<"Max mass after cooling: "<<*m_max<<endl<<endl;
            }
            else{fail_process_check_mmax_pt_two_branch("check_mmax: find max mass error(cool_eos not success)");}
        }
        else{
            if (*m_max>maxm_tov) {fail_process_check_mmax_pt_two_branch("check_mmax: find max mass error without cooling(M>max_m_tov)");}      
        }
        return true;
    }
    *h1 = -1, *h2 = -1, *h3 = -1, *h_max = -1, *m1 = 0., *m2 = 0., *m3 = 0., *m_max=0.;
    return false;
}


//----------------------------------------------------------------------------------------------
//       Find global properties using another global properties
//----------------------------------------------------------------------------------------------


/**
* @brief Find a global property a with another global property b.
* @param known_aim Global_b.
* @param h_i Initial guess.
* @param h_closest Pointer to the closest enthalpy, return -1 if finding error.
* @param unknown Pointer to the global_a to be found, return 0 if finding error.
* @param get_type 1. using mass to find lambda; 2. using mass to find radius; 3. using lambda to find mass.
* @param use_user_start_point Whether to use user specified initial guess h_i or use system default value.
* @retval bool whether finding process correctly worked.
*/
bool find_closest_global_property(double m_aim, double h_i, double h_max, double *h_closest, double *lambda, bool use_user_start_point){
//find nearest mass, return lambda;
    bool jump_back=false, force_jump_back=false;
    state_type x_result(3), x_guard(3), x_check_again(3);
    double h_guard=4e-3;
    double h, h_init, h_save, jump_size = 0.1;
    double M=0., M_guard=0., M_check_again=0., M_temp=0.;
    double minimum_allowed_h = EOS->minimum_allowed_h;
    if (use_user_start_point) h_init = h_i;
    else h_init = EOS->suggested_hc_star;
    if(verbose){
        cout<<endl<<"find closest aim (known maximum mass): "<<m_aim<<", user specified start point: "<<h_i<<", suggested start point: "<<EOS->suggested_hc_star;
        cout<<", minimum_allowed_h: "<<minimum_allowed_h<<", use user start point: "<<use_user_start_point<<", h start: "<<h_init<<", h end: "<<h_max<<endl;
    }
    h = h_init, h_save = h_init;
    while(h>minimum_allowed_h and h<h_max and jump_size>=1e-3){
        try{
            x_result = int_whole_star(h);
            M = x_result[0];
            x_guard = int_whole_star(h+h_guard);
            M_guard = x_guard[0];
        }
        catch (exception &){
            cout<<"Invalid value encountered in find nearest mass: h="<<h<<endl;
            M=0; M_guard=0;
        }
        if (verbose) cout<<h<<"    "<<M_guard<<"     "<<M<<endl;
        jump_back = (M==0 or M_guard==0 or std::isnan(M_guard) or std::isnan(M) or max(M, M_guard)>m_aim);
        force_jump_back = false;
        if ((not jump_back) and M_guard<=M){
            try {
                x_check_again = int_whole_star(h+0.01);
                M_check_again = x_check_again[0];
                if (verbose) cout<<"M_check_again: "<<M_check_again<<endl;
                if (M_check_again<M) {jump_back = true; if (max(M, M_guard)<M_temp) force_jump_back=true;} 
            }
            catch (exception &) { if (verbose) cout<<"Invalid value encountered in find maximum mass: h="<<h+0.01<<endl;}
        }
        h_save = h;
        if (jump_back) {
            if (force_jump_back) h -= jump_size;
            h -= jump_size;
            jump_size /= 5.;
            h += jump_size;
        }
        else h += jump_size;
        M_temp = max(M, M_guard);
    }
    if (verbose) cout<<"find h: "<<h_save<<",    with M:"<<M<<endl<<endl;
    if (std::isnan(M) or h_save<=minimum_allowed_h or force_jump_back) {*h_closest = -1; *lambda=0; return false;}
    else {
        *h_closest = h_save;
        if (x_result[0]>0. and x_result[1]>0.) *lambda = cal_lambda(x_result[0]/x_result[1], x_result[2]);
        else *lambda = 0.;
        return true;
    }
}


/**
* @brief Find a global property a with another global property b after knowing the maximum central specific enthalpy h_max.
* @param known_aim Global_b.
* @param h_i Initial guess.
* @param h_max Known maximum enthalpy.
* @param h_closest Pointer to the closest enthalpy, return -1 if finding error.
* @param unknown Pointer to the global_a to be found, return 0 if finding error.
* @param get_type 1. using mass to find lambda; 2. using mass to find radius; 3. using lambda to find mass.
* @param use_user_start_point Whether to use user specified initial guess h_i or use system default value.
* @retval bool whether finding process correctly worked.
*/
bool find_closest_global_property_with_maxm_known(double known_aim, double h_i, double h_max, double *h_closest, double *unknown, int get_type, bool use_user_start_point){
    int iter = 0;
    bool jump_back = false;
    state_type x_result(3);
    double h_guard=4e-3;
    double h, h_init, h_save, jump_size=0.1, test_known=0.;
    double minimum_allowed_h = EOS->minimum_allowed_h;
    if (use_user_start_point) h_init = h_i;
    else h_init = EOS->suggested_hc_star;
    #define fail_process_find_global_prop() *h_closest = -1; *unknown=0; return false;
    if(get_type!=1 and get_type!=2 and get_type!=3) {cout<<"not known parameter value of get_type(in func find_closest_with_mmax_known): "<<get_type<<endl; fail_process_find_global_prop();}
    if(verbose){
        cout<<endl<<"find closest aim (known maximum mass): "<<known_aim<<", user specified start point: "<<h_i<<", suggested start point: "<<EOS->suggested_hc_star;
        cout<<", minimum_allowed_h: "<<minimum_allowed_h<<", use user start point: "<<use_user_start_point<<", h start: "<<h_init<<", h end: "<<h_max<<endl;
    }
    h = h_init, h_save = h_init;
    while(h>minimum_allowed_h and h<h_max and jump_size>=1e-4){
        try{
            x_result = int_whole_star(h);
            if(get_type==1 or get_type==2) test_known = x_result[0];
            else if(get_type==3) test_known = cal_lambda(x_result[0]/x_result[1], x_result[2]);
            else {}
        }
        catch (exception &){
            cout<<"Invalid value encountered in find nearest mass: h="<<h<<endl;
            test_known=0;
        }
        if (verbose) cout<<"step:"<<iter<<"    ,h: "<<h<<"     ,M: "<<test_known<<endl;
        if(get_type==3) jump_back = (test_known==0 or std::isnan(test_known) or test_known<known_aim);
        else jump_back = (test_known==0 or std::isnan(test_known) or test_known>known_aim);
        if (h!=h_init) jump_back |= (h+jump_size>h_max);
        else { if (jump_back and (iter==0)) jump_size = h-minimum_allowed_h;}
        h_save = h;
        if (jump_back) {
            h -= jump_size;
            jump_size /= 5.;
            h += jump_size;
        }
        else h += jump_size;
        iter += 1;
    }
    if (verbose) cout<<"aim: "<<known_aim<<", h: "<<h_save<<", used steps: "<<iter<<", test_known:"<<test_known<<endl;
    if (std::isnan(test_known) or h_save<=minimum_allowed_h) {fail_process_find_global_prop();}
    else {
        *h_closest = h_save;
        if (x_result[0]>0. and x_result[1]>0.) {
            if(get_type==1) {*unknown = cal_lambda(x_result[0]/x_result[1], x_result[2]);}
            else if(get_type==2) {*unknown = x_result[1]*r_trans;}
            else if(get_type==3) {*unknown = x_result[0];}
            else {}
        }
        else *unknown = 0.;
        return true;
    }
}


bool find_closest_global_property_with_maxm_known_old(double known_aim, double h_i, double h_max, double *h_closest, double *unknown, int get_type, bool use_user_start_point){
    int iter = 0;
    bool jump_back = false;
    state_type x_result(3);
    double h_guard=4e-3;
    double h, h_init, h_save, jump_size=0.1, test_known=0.;
    double minimum_allowed_h = EOS->minimum_allowed_h;
    if (use_user_start_point) h_init = h_i;
    else h_init = EOS->suggested_hc_star;
    #define fail_process_find_global_prop() *h_closest = -1; *unknown=0; return false;
    if(get_type!=1 and get_type!=2 and get_type!=3) {cout<<"not known parameter value of get_type(in func find_closest_with_mmax_known): "<<get_type<<endl; fail_process_find_global_prop();}
    if(verbose){
        cout<<endl<<"find closest aim (known maximum mass): "<<known_aim<<", user specified start point: "<<h_i<<", suggested start point: "<<EOS->suggested_hc_star;
        cout<<", minimum_allowed_h: "<<minimum_allowed_h<<", use user start point: "<<use_user_start_point<<", h start: "<<h_init<<", h end: "<<h_max<<endl;
    }
    h = h_init, h_save = h_init;
    while(h>minimum_allowed_h and h<h_max and jump_size>=1e-4){
        try{
            x_result = int_whole_star(h);
            if(get_type==1 or get_type==2) test_known = x_result[0];
            else if(get_type==3) test_known = cal_lambda(x_result[0]/x_result[1], x_result[2]);
            else {}
        }
        catch (exception &){
            cout<<"Invalid value encountered in find nearest mass: h="<<h<<endl;
            test_known=0;
        }
        if (verbose) cout<<"step:"<<iter<<"    ,h: "<<h<<"     ,M: "<<test_known<<endl;
        if(get_type==3) jump_back = (test_known==0 or std::isnan(test_known) or test_known<known_aim);
        else jump_back = (test_known==0 or std::isnan(test_known) or test_known>known_aim);
        if (h!=h_init) jump_back |= (h+jump_size>h_max);
        else { if (jump_back and (iter==0)) jump_size = h-minimum_allowed_h;}
        h_save = h;
        if (jump_back) {
            h -= jump_size;
            jump_size /= 5.;
            h += jump_size;
        }
        else h += jump_size;
        iter += 1;
    }
    if (verbose) cout<<"aim: "<<known_aim<<", h: "<<h_save<<", used steps: "<<iter<<", test_known:"<<test_known<<endl;
    if (std::isnan(test_known) or h_save<=minimum_allowed_h) {fail_process_find_global_prop();}
    else {
        *h_closest = h_save;
        if (x_result[0]>0. and x_result[1]>0.) {
            if(get_type==1) {*unknown = cal_lambda(x_result[0]/x_result[1], x_result[2]);}
            else if(get_type==2) {*unknown = x_result[1]*r_trans;}
            else if(get_type==3) {*unknown = x_result[0];}
            else {}
        }
        else *unknown = 0.;
        return true;
    }
}


/**
* @brief In a BNS system: Find nearest known global properties, return corresponding unknown global ones.
* @param known1 Known global property 1 to find with.
* @param known2 Known global property 2 to find with.
* @param unknown1 Unknown global property 1 to be found.
* @param unknown2 Unknown global property 2 to be found.
* @param h_max Pointer to stare maximum enthalpy, this should be the output also, this means this function also try to find the maximum mass.
* @param get_type 1. using mass to find lambda; 2. using mass to find radius; 3. using lambda to find mass.
* @retval bool Whether finding process correctly worked.
*/
bool get_unknowns_from_knowns(double known1, double known2, double *unknown1, double *unknown2, double *h_max, int get_type){
    double h1=0., h2=0., M_max=0., L_max=0.;
    *unknown1 = 0.; *unknown2 = 0.; *h_max=0.;
    if(not check_mmax_gd(h_max, &M_max)) return false;
    if(get_type==1 or get_type==2){
        if ((known1>M_max) or (known2>M_max)) return false;
    }
    else if (get_type==3){
        get_mrl(*h_max);
        L_max = mrl_result[2];
        if ((known1<L_max) or (known2<L_max)) return false;
    }
    else{
        cout<<"not known parameter value of get_type(in func get_unkowns_from_knowns): "<<get_type<<endl;
        return false;
    }
    if (known1<known2){
        if (find_closest_global_property_with_maxm_known(known1, EOS->suggested_hc_star, *h_max, &h1, unknown1, get_type)) find_closest_global_property_with_maxm_known(known2, h1, *h_max, &h2, unknown2, get_type);
        else find_closest_global_property_with_maxm_known(known2, EOS->suggested_hc_star, *h_max, &h2, unknown2, get_type);
    }
    else if(known1>known2){
        if (find_closest_global_property_with_maxm_known(known2, EOS->suggested_hc_star, *h_max, &h2, unknown2, get_type)) find_closest_global_property_with_maxm_known(known1, h2, *h_max, &h1, unknown1, get_type);
        else find_closest_global_property_with_maxm_known(known1, EOS->suggested_hc_star, *h_max, &h1, unknown1, get_type);
    }
    else{
        find_closest_global_property_with_maxm_known(known1, EOS->suggested_hc_star, *h_max, &h1, unknown1, get_type);
        *unknown2 = *unknown1;
    }
    if (*unknown1>0. and *unknown2>0.) return true;
    else return false;
}


#endif
