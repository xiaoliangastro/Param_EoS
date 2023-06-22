#ifndef EOS_HPP
#define EOS_HPP


#include<map>
#include<cmath>
#include<string>
#include<iostream>
#include<exception>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_roots.h>
#include<boost/math/interpolators/pchip.hpp>
#include<boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include<boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#include<boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include<boost/numeric/odeint/integrate/integrate_const.hpp>
#include"global_variable_constants.hpp"


using boost::math::interpolators::pchip;
using namespace boost::numeric::odeint;
using namespace std;
typedef std::vector< double > state_type;



state_type cal_eborders(state_type params, state_type rho_borders);
void cal_eos(const state_type &x, state_type &dxdt, double h);
typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
controlled_stepper_type controlled_stepper_cal_eos(default_error_checker< double, range_algebra, default_operations >( 1e-15, 1e-9, int_a_x, int_a_dxdt ) );//controlled_stepper_type


////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                            //
//                                   Virtural Base EoS Classes                                //
//                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////



class EoS{
    public:
    	double minimum_allowed_h = 0.0316;
    	double eos_table_max_h=2;
    	double suggested_hc_star = minimum_allowed_h+0.1;///< suggested start hc to find neutron star for specific properties.
    	state_type params_eos;
    	vector< pchip<state_type> > eos_table_function_h_base;           // reachable only after calling update_eos()
        vector< pchip<state_type> > eos_table_function_rho_base;         // reachable only after calling update_eos()
        vector< pchip<state_type> > eos_table_function_e_base;           // reachable only after calling update_eos()
        state_type eos_table_h, eos_table_p, eos_table_e, eos_table_rho; // reachable only after calling update_eos()
    	EoS(){}
    	virtual bool eos_valid() { return false; }
    	virtual bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]) { return false; }
    	virtual double eos_messenger(int ask_for_information_type){ return 0; }
    	virtual double ph(double h)=0;
    	virtual double eh(double h)=0;
    	virtual double rhoh(double h, double p, double e)=0;
    	virtual double gamma(double h, double p, double e)=0;
    	virtual double cs_sq(double p, double e, double Gamma) { return (p*Gamma)/(e+p); }
    	virtual ~EoS(){}
};


class EoS_tablulated: public EoS{
    public:
    	map<double, double[3]>::iterator coeftb_idx;///< index iterator-get the corresponding interval of h to use corresponding interpolation coefficients.
        map<double, double[3]> coeftb;///< look up table for the coefficient of p, e and c, taking h as key.
        void check_coeftb_idx(double h){
        	if (coeftb_idx==coeftb.end()){
        		coeftb_idx = --coeftb.upper_bound(h);
        	}
        	else{
        	    if (not ((coeftb_idx->first>h) and ((--coeftb_idx)->first<=h))) coeftb_idx = --coeftb.upper_bound(h);
        	}
        }
    	double eh_tabulated(double h){
    		check_coeftb_idx(h);
            double co_p = coeftb_idx->second[0], co_e = coeftb_idx->second[1], co_c = coeftb_idx->second[2], co_h = coeftb_idx->first;
            if (coeftb_idx!= coeftb.begin()) return co_e*pow((co_e+co_p)/co_p*exp((co_c-1.0)/co_c*(h-co_h))-co_e/co_p, 1.0/(co_c-1.0));
            else return co_e*pow(co_e/co_p*(exp(2.0*h/5.0)-1.0), 1.5);
        }
        double ph_tabulated(double h){
        	check_coeftb_idx(h);
        	return coeftb_idx->second[0]*pow(eh_tabulated(h)/coeftb_idx->second[1], coeftb_idx->second[2]);
        }
        double rhoh_tabulated(double h, double p, double e){ return (p+e)/exp(h); }
        double gamma_tabulated(double h, double p, double e){
        	check_coeftb_idx(h);
        	return coeftb_idx->second[2]*(p/e+1.0);
        }
};


class EoS_hybrid: public EoS_tablulated{
    public:
    	double h_border=h_0; ///< the border that seperate apart the interpolated part and the tabulated part.
        EoS_hybrid(double ee[], double pp[], const int len_eos){
        	init_interp_coeff_through_ep_array(ee, pp, len_eos);
        	update_h0_related_variables();
        	h_border=h_0;
        }
		/** @brief Calculate pressure using h only. */
		double ph(double h){
			if (h>h_border) return ph_interp(h);
			else return ph_tabulated(h);
		}
		/** @brief Calculate energy density using h only. */
		double eh(double h){
			if (h>h_border) return eh_interp(h);
			else return eh_tabulated(h);
		}
		/** @brief Calculate mass density using h only. */
		double rhoh(double h, double p, double e){
			if (h>h_border) return rhoh_interp(h);
			else return rhoh_tabulated(h, p, e);
		}
		/** @brief Calculate gamma. */
		double gamma(double h, double p, double e){
			if (h>h_border) return gamma_interp(h, p, e);
			else return gamma_tabulated(h, p, e);
		}
        double ph_interp(double h){ return eos_table_function_h_base[0](h); }
        double eh_interp(double h){ return eos_table_function_h_base[1](h); }
        double rhoh_interp(double h){ return eos_table_function_h_base[2](h); }
        virtual double gamma_interp(double h, double p, double e) {return 0;}

		/** @brief update global variables h_0, p_0, e_0, and rho_0 with the last term of the coefficient table */
		void update_h0_related_variables(){
		//must be called before function update_eos
		    auto iter = coeftb.end();
		    iter--;
		    h_0 = iter->first; p_0 = iter->second[0];
		    e_0 = iter->second[1]; rho_0 = (e_0+p_0)/exp(h_0);
		    if (verbose) cout<<"h_0: "<<h_0<<", p_0: "<<p_0<<", e_0: "<<e_0<<", rho_0: "<<rho_0<<endl;
		}

		/** @brief calculate eos tables as list of interpolating functions and store them */
		void cal_eos_table(double h_max, const char init_function_type[]){
		    state_type eos_tb{p_0, e_0, rho_0};
		    vector<state_type> x_o;
		    eos_table_h.clear(); eos_table_p.clear(); eos_table_e.clear(); eos_table_rho.clear();
		    try {
			    if (consid_const_inter_step) integrate_const(error_stepper_type(), cal_eos, eos_tb, h_0, h_max, sg_const_step, push_back_state_and_time(x_o, eos_table_h));
		        else integrate_adaptive(controlled_stepper_cal_eos, cal_eos, eos_tb, h_0, h_max, sg_step, push_back_state_and_time(x_o, eos_table_h));
		        for (int i=0; i<eos_table_h.size(); i++){ eos_table_p.push_back(x_o[i][0]), eos_table_e.push_back(x_o[i][1]), eos_table_rho.push_back(x_o[i][2]); }
		        if (init_function_type[0]=='1'){
		            eos_table_function_h_base.clear();
		            state_type h_bk1(eos_table_h), h_bk2(eos_table_h), h_bk3(eos_table_h);
		            state_type p_bk(eos_table_p), e_bk(eos_table_e), rho_bk(eos_table_rho);
		            auto function_h_p = pchip(std::move(h_bk1), std::move(p_bk));
		            auto function_h_e = pchip(std::move(h_bk2), std::move(e_bk));
		            auto function_h_rho = pchip(std::move(h_bk3), std::move(rho_bk));
		            eos_table_function_h_base.push_back(function_h_p);
		            eos_table_function_h_base.push_back(function_h_e);
		            eos_table_function_h_base.push_back(function_h_rho);
		        }
		        if (init_function_type[1]=='1'){
		            eos_table_function_e_base.clear();
		            state_type e_bk1(eos_table_e), e_bk2(eos_table_e), e_bk3(eos_table_e);
		            state_type h_bk(eos_table_h), p_bk(eos_table_p), rho_bk(eos_table_rho);
		            auto function_e_h = pchip(std::move(e_bk1), std::move(h_bk));
		            auto function_e_p = pchip(std::move(e_bk2), std::move(p_bk));
		            auto function_e_rho = pchip(std::move(e_bk3), std::move(rho_bk));
		            eos_table_function_e_base.push_back(function_e_h);
		            eos_table_function_e_base.push_back(function_e_p);
		            eos_table_function_e_base.push_back(function_e_rho);
		        }
		        if (init_function_type[2]=='1'){
		            eos_table_function_rho_base.clear();
		            state_type rho_bk1(eos_table_rho), rho_bk2(eos_table_rho), rho_bk3(eos_table_rho);
		            state_type h_bk(eos_table_h), p_bk(eos_table_p), e_bk(eos_table_e);
		            auto function_rho_h = pchip(std::move(rho_bk1), std::move(h_bk));
		            auto function_rho_p = pchip(std::move(rho_bk2), std::move(p_bk));
		            auto function_rho_e = pchip(std::move(rho_bk3), std::move(e_bk));
		            eos_table_function_rho_base.push_back(function_rho_h);
		            eos_table_function_rho_base.push_back(function_rho_p);
		            eos_table_function_rho_base.push_back(function_rho_e);
		        }
		        if (vverbose) {cout<<"initiate type of interpolation function : "<<string(init_function_type)<<endl;}
		    }
		    catch (exception & except) {
		        cout<<"error encountered in cal_eos_table: "<<except.what()<<endl;
		        throw;
		    }
		}

        /** @brief transform ep table to get interpolation coefficients (see Appendix B of 10.1103/PhysRevD.86.084003) */
		void init_interp_coeff_through_ep_array(double ee[], double pp[], const int len_eos){
			int lm1 = len_eos-1;
		    double cc[len_eos], hh[len_eos];
		    coeftb.clear();
		    coeftb[0.0][0] = pp[0], coeftb[0.0][1] = ee[0], coeftb[0.0][2] = 5.0/3.0;
		    hh[0] = (5.0/2.0)*log((ee[0]+pp[0])/ee[0]);
		    for (int i=0; i<lm1; i++){
		        cc[i] = log(pp[i+1]/pp[i])/log(ee[i+1]/ee[i]);
		        hh[i+1] = hh[i]+(cc[i]/(cc[i]-1))*log((ee[i]*(ee[i+1]+pp[i+1]))/(ee[i+1]*(ee[i]+pp[i])));
		        coeftb[hh[i]][0] = pp[i], coeftb[hh[i]][1] = ee[i], coeftb[hh[i]][2] = cc[i];
		    }
		    cc[lm1] = cc[lm1-1];
		    coeftb[hh[lm1]][0] = pp[lm1], coeftb[hh[lm1]][1] = ee[lm1], coeftb[hh[lm1]][2] = cc[lm1];
		    coeftb_idx=coeftb.end(); //initialize coeftb_idx
		    //map<double, double[3]>::iterator id;
			//cout<<"fuck here!"<<endl;
		    //for (id = coeftb.begin(); id!= coeftb.end(); id++) cout<<"indexes:"<<id->first<<"\t\t"<<id->second[0]<<"\t\t"<<id->second[1]<<"\t\t"<<id->second[2]<<endl;
		}
};



////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                            //
//                                   Eos of various kinds                                     //
//                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////



//----------------------------------------------------------------------------------------------
//                                   EoS of Tabulated Only
//----------------------------------------------------------------------------------------------


class EoS_tabulated_only: public EoS_tablulated{
    public:
    	EoS_tabulated_only(double ee[], double pp[], const int len_eos){
    		init_interp_coeff_through_ep_array(ee, pp, len_eos, false);
    	}
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		// expected parameters: ppar: [e, p], extra_par: length of e or p
    		double ee[extra_par], pp[extra_par];
    		bool interpolation = ((init_function_type[0]=='1') or (init_function_type[1]=='1') or (init_function_type[2]=='1'));
    		for (int i=0; i<extra_par; i++) {ee[i]=ppar[i]; pp[i]=ppar[extra_par+i];}
    		try{
    			if (interpolation){
		    		init_interp_coeff_through_ep_array(ee, pp, extra_par, true);
		    		cal_eos_table(h_max, init_function_type);
		    	}
		    	else init_interp_coeff_through_ep_array(ee, pp, extra_par, false);
	    	}
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	return false;
    	    }
    	    return true;
    	}
    	double eh(double h){ return eh_tabulated(h);}
    	double ph(double h){ return ph_tabulated(h);}
		double rhoh(double h, double p, double e){ return rhoh_tabulated(h, p, e); }
    	double gamma(double h, double p, double e){ return gamma_tabulated(h, p, e);}
		/** @brief calculate eos tables as list of interpolating functions and store them */
		void cal_eos_table(double h_max, const char init_function_type[]){
		    try {
		        if (init_function_type[0]=='1'){
		            eos_table_function_h_base.clear();
		            state_type h_bk1(eos_table_h), h_bk2(eos_table_h), h_bk3(eos_table_h);
		            state_type p_bk(eos_table_p), e_bk(eos_table_e), rho_bk(eos_table_rho);
		            auto function_h_p = pchip(std::move(h_bk1), std::move(p_bk));
		            auto function_h_e = pchip(std::move(h_bk2), std::move(e_bk));
		            auto function_h_rho = pchip(std::move(h_bk3), std::move(rho_bk));
		            eos_table_function_h_base.push_back(function_h_p);
		            eos_table_function_h_base.push_back(function_h_e);
		            eos_table_function_h_base.push_back(function_h_rho);
		        }
		        if (init_function_type[1]=='1'){
		            eos_table_function_e_base.clear();
		            state_type e_bk1(eos_table_e), e_bk2(eos_table_e), e_bk3(eos_table_e);
		            state_type h_bk(eos_table_h), p_bk(eos_table_p), rho_bk(eos_table_rho);
		            auto function_e_h = pchip(std::move(e_bk1), std::move(h_bk));
		            auto function_e_p = pchip(std::move(e_bk2), std::move(p_bk));
		            auto function_e_rho = pchip(std::move(e_bk3), std::move(rho_bk));
		            eos_table_function_e_base.push_back(function_e_h);
		            eos_table_function_e_base.push_back(function_e_p);
		            eos_table_function_e_base.push_back(function_e_rho);
		        }
		        if (init_function_type[2]=='1'){
		            eos_table_function_rho_base.clear();
		            state_type rho_bk1(eos_table_rho), rho_bk2(eos_table_rho), rho_bk3(eos_table_rho);
		            state_type h_bk(eos_table_h), p_bk(eos_table_p), e_bk(eos_table_e);
		            auto function_rho_h = pchip(std::move(rho_bk1), std::move(h_bk));
		            auto function_rho_p = pchip(std::move(rho_bk2), std::move(p_bk));
		            auto function_rho_e = pchip(std::move(rho_bk3), std::move(e_bk));
		            eos_table_function_rho_base.push_back(function_rho_h);
		            eos_table_function_rho_base.push_back(function_rho_p);
		            eos_table_function_rho_base.push_back(function_rho_e);
		        }
		        if (vverbose) {cout<<"initiate type of interpolation function : "<<string(init_function_type)<<endl;}
		    }
		    catch (exception & except) {
		        cout<<"error encountered in cal_eos_table: "<<except.what()<<endl;
		        throw;
		    }
		}
		
        /** @brief transform ep table to get interpolation coefficients (see Appendix B of 10.1103/PhysRevD.86.084003) */
		void init_interp_coeff_through_ep_array(double ee[], double pp[], const int len_eos, bool interpolation){
			int lm1 = len_eos-1;
		    double cc[len_eos], hh[len_eos], rr[len_eos];
		    coeftb.clear();
		    coeftb[0.0][0] = pp[0], coeftb[0.0][1] = ee[0], coeftb[0.0][2] = 5.0/3.0;
		    hh[0] = (5.0/2.0)*log((ee[0]+pp[0])/ee[0]);
		    if (interpolation) eos_table_h.clear(); eos_table_p.clear(); eos_table_e.clear(); eos_table_rho.clear();
		    for (int i=0; i<lm1; i++){
		        cc[i] = log(pp[i+1]/pp[i])/log(ee[i+1]/ee[i]);
		        hh[i+1] = hh[i]+(cc[i]/(cc[i]-1))*log((ee[i]*(ee[i+1]+pp[i+1]))/(ee[i+1]*(ee[i]+pp[i])));
		        coeftb[hh[i]][0] = pp[i], coeftb[hh[i]][1] = ee[i], coeftb[hh[i]][2] = cc[i];
		        if (interpolation){
		        	rr[i] = (ee[i]+pp[i])/exp(hh[i]);
		    	    eos_table_h.push_back(hh[i]);
		    	    eos_table_e.push_back(ee[i]);
		    	    eos_table_p.push_back(pp[i]);
		    	    eos_table_rho.push_back(rr[i]);
		        }
		    }
		    cc[lm1] = cc[lm1-1];
		    coeftb[hh[lm1]][0] = pp[lm1], coeftb[hh[lm1]][1] = ee[lm1], coeftb[hh[lm1]][2] = cc[lm1];
		    coeftb_idx=coeftb.end(); //initialize coeftb_idx
		    //map<double, double[3]>::iterator id;
		    //for (id = coeftb.begin(); id!= coeftb.end(); id++) cout<<"indexes:"<<id->first<<"\t\t"<<id->second[0]<<"\t\t"<<id->second[1]<<"\t\t"<<id->second[2]<<endl;
		}
};


//----------------------------------------------------------------------------------------------
//                                   Gamma Piecewise Model
//----------------------------------------------------------------------------------------------


class EoS_piecewise_gamma: public EoS_hybrid{
    public:
    	state_type e_borders;
    	state_type rho_borders{1.*rho_sat, 1.85*rho_sat, 3.7*rho_sat, 7.4*rho_sat};
    	using EoS_hybrid::EoS_hybrid;
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		params_eos.clear(); e_borders.clear();
    		eos_table_max_h = h_max;
    		for (int i=0; i<len_par; i++) params_eos.push_back(ppar[i]);
    		try{
    		    e_borders = cal_eborders(params_eos, rho_borders);
    	        cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	return false;
    	    }
    	    return true;
    	}
		double gamma_interp(double h, double p, double e){
		    auto upper = std::lower_bound(e_borders.begin()+1, e_borders.end(), e);
		    int idx = std::distance(e_borders.begin(), upper)-1;
		    if (upper==e_borders.end()) {idx -= 1;}
		    return params_eos[idx];
		}
};


//----------------------------------------------------------------------------------------------
//                                 Pressure Piecewise Model
//----------------------------------------------------------------------------------------------


class EoS_piecewise_pressure: public EoS_hybrid{
    public:
    	state_type e_borders;
    	state_type rho_borders{1.*rho_sat, 1.85*rho_sat, 3.7*rho_sat, 7.4*rho_sat};
    	using EoS_hybrid::EoS_hybrid;
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		// expected ppar: pressures at rho_borders, thus if you want more pieces, just change rho_borders and transfer the corresponding pressures
    		params_eos.clear(); e_borders.clear();
    		eos_table_max_h = h_max;
    		for (int i=0; i<len_par; i++) params_eos.push_back(ppar[i]);
    		try{
    		    e_borders = cal_eborders(params_eos, rho_borders);
    	        cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	return false;
    	    }
    	    return true;
    	}
		/** @brief Calculate gamma in the pressure piecewise model. */
		double gamma_interp(double h, double p, double e){
		    auto upper = std::lower_bound(e_borders.begin()+1, e_borders.end(), e);
		    int idx = std::distance(e_borders.begin(), upper)-1;
		    if (upper==e_borders.end()) {idx -= 1;}
		    if(idx==0) return log(params_eos[0]/p_0)/log(rho_borders[0]/rho_0);
		    else return log(params_eos[idx]/params_eos[idx-1])/log(rho_borders[idx]/rho_borders[idx-1]);
		}
};


//----------------------------------------------------------------------------------------------
//                                H-based Spectral(causal) Model
//----------------------------------------------------------------------------------------------


class EoS_spectral_hbase: public EoS_hybrid{
    public:
    	bool spectral_causal;
    	double h_base;
        EoS_spectral_hbase(double ee[], double pp[], const int len_eos, bool causal): EoS_hybrid(ee, pp, len_eos){
        	spectral_causal = causal;
        }
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		// expected ppar: g0, g1, g2, g3, h_base
    		params_eos.clear();
    		eos_table_max_h = h_max;
	        for (int i=0; i<extra_par; i++) params_eos.push_back(ppar[i]);
	        h_base = ppar[extra_par];
    		try{
    	        cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	return false;
    	    }
    	    return true;
    	}
		/** @brief Calculate gamma in the spectral expansion model. */
		double gamma_interp(double h, double p, double e){
		    double loggamma = 0.0, ret;
		    for(int j = 0; j<params_eos.size(); j++) loggamma += params_eos[j]*pow(log(h/h_base), j);
		    ret = exp(loggamma);
		    if (spectral_causal) ret = (e+p)/((ret+1.)*p);
		    return ret;
		}
};


//----------------------------------------------------------------------------------------------
//                             Piecewise_gamma + Phase Transition Model
//----------------------------------------------------------------------------------------------


class EoS_pwg_pt: public EoS_hybrid{
    public:
    	state_type e_borders;
    	using EoS_hybrid::EoS_hybrid;
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		// expected ppar: eb1, eb2, eb3, gamma1, gamma2, gamma3
    		params_eos.clear(); e_borders.clear();
	        int len_border = extra_par;
    		eos_table_max_h = h_max;
	        for (int i=0; i<len_border; i++) e_borders.push_back(ppar[i]);
	        for (int i=len_border; i<len_par; i++) params_eos.push_back(ppar[i]);
    		try{
    	        cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	return false;
    	    }
	        return true;
	    }
		/** @brief Calculate gamma in the phase transition model. */
		double gamma_interp(double h, double p, double e){
		    auto upper = std::lower_bound(e_borders.begin()+1, e_borders.end(), e);
		    int idx = std::distance(e_borders.begin(), upper)-1;
		    if(idx<=2) return params_eos[idx];
		    else if(idx==3) return e/p+1.;
		    else return (e/p+1.)*(1./3.);
		}
};


//----------------------------------------------------------------------------------------------
//                                   MIT bag model
//----------------------------------------------------------------------------------------------


class EoS_MIT_bag: public EoS{
    public:
    	static constexpr double m_qk = 100.; ///< mass of strange quark in MeV
        static constexpr double m_br = 930.; ///< rest baryon mass on stellar surface in MeV
        struct quark_params{double a4, Beff, egap;};
        EoS_MIT_bag(){}
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		// expected parameters: a4, Beff, egap
    		params_eos.clear();
    		for (int i=0; i<len_par; i++) params_eos.push_back(ppar[i]);
    		try{
    		    params_eos.push_back(solve_surface_h(-1.0, 0.7));///< push back where the pressure(h)=0
    		}
    		catch (exception &except){
                cout<<"failed to initialize MIT bag model: "<<except.what()<<endl;
                return false;
            }
            minimum_allowed_h = params_eos[3];
            suggested_hc_star = minimum_allowed_h+0.1;
            if (verbose) cout<<"surf_h: "<<params_eos[3]<<", suggested h: "<<suggested_hc_star<<endl;
            return true;
    	}
        double ph(double h){return ph_qks_org(h)*MeV4_to_dmls;}
        double eh(double h){return eh_qks_org(h, ph_qks_org(h))*MeV4_to_dmls;}
        double rhoh(double h, double p, double e){return 0.;} // not sure it is n_CFL*(m_u+m_d+m_s)/3 or not!
        double gamma(double h, double p, double e){
            //As soon as the unit of p and e are the same, it will be OK, no matter it is MeV^4 or dimensionless, because they will cancel each other.
            double mu = exp(h)*m_br/3., n = n_CFL(mu, params_eos[0], params_eos[1], params_eos[2]);
            double mu_dot_n_mu =  -(3.*pow(m_qk,2.)*mu/2./M_PI/M_PI-9.*pow(mu,3.)/M_PI/M_PI \
            +3.*pow(m_qk,4.)/8./M_PI/M_PI/mu+(1.-params_eos[0])*9.*pow(mu,3.)/M_PI/M_PI-6.*pow(params_eos[2],2.)*mu/M_PI/M_PI)/3.;
            return (e+p)*n/(mu_dot_n_mu*p);
        }
        double ph_qks_org(double h){return -Omega_CFL(exp(h)*m_br/3., params_eos[0], params_eos[1], params_eos[2]);}
        double eh_qks_org(double h, double p){return -p+3.*n_CFL(exp(h)*m_br/3., params_eos[0], params_eos[1], params_eos[2])*(exp(h)*m_br/3.);}
    	double Omega_CFL(double mu, double a4, double Beff, double egap){return 3.*pow(m_qk,2.)*pow(mu,2.)/4./M_PI/M_PI- \
	                            3.*pow(mu,4.)/4./M_PI/M_PI-(1.-12.*log(m_qk/2./mu))/32./M_PI/M_PI*pow(m_qk,4.)+(1.-a4)*3.*pow(mu,4.)/4./M_PI/M_PI-3.*pow(egap,2.)*pow(mu,2.)/M_PI/M_PI+Beff;}
        double n_CFL(double mu, double a4, double Beff, double egap){return -(3.*pow(m_qk,2.)*mu/2./M_PI/M_PI-3.*pow(mu,3.)/M_PI/M_PI \
	                            -3.*pow(m_qk,4.)/8./M_PI/M_PI/mu+(1.-a4)*3.*pow(mu,3.)/M_PI/M_PI-6.*pow(egap,2.)*mu/M_PI/M_PI)/3.;}
        static double minfunc_qks(double h, void *params){
			struct quark_params *p = static_cast<struct quark_params *>(params);
			double mu=exp(h)*m_br/3., a4=p->a4, Beff=p->Beff, egap=p->egap;
			return 3.*pow(m_qk,2.)*pow(mu,2.)/4./M_PI/M_PI-3.*pow(mu,4.)/4./M_PI/M_PI-(1.-12.*log(m_qk/2./mu))/32./M_PI/M_PI*pow(m_qk,4.)+\
			             (1.-a4)*3.*pow(mu,4.)/4./M_PI/M_PI-3.*pow(egap,2.)*pow(mu,2.)/M_PI/M_PI+Beff; //Omega_CFL
		}
		double solve_surface_h(double lower_h, double upper_h, int max_iter=100){
		    //initiate
		    int status;
		    int iter = 0;
		    double guess = (lower_h+upper_h)/2., h_low, h_up, h_tol = 1e-8;
		    const gsl_root_fsolver_type *ST;
		    gsl_root_fsolver *solver;
		    gsl_function Func;
		    struct quark_params params = {params_eos[0], params_eos[1], params_eos[2]};
		    Func.function = &minfunc_qks;
		    Func.params = &params;
		    ST = gsl_root_fsolver_brent;
		    solver = gsl_root_fsolver_alloc(ST);
		    //check root exist
		    double f_lower = minfunc_qks(lower_h, &params), f_upper = minfunc_qks(upper_h, &params);
		    if ((f_lower < 0.0 && f_upper < 0.0) || (f_lower > 0.0 && f_upper > 0.0)){
		        cout<<"f_lower: "<<f_lower<<", f_upper: "<<f_upper<<endl;
		        throw std::runtime_error("Do not have root in given region!(in function solve_surface_h)");
		    }
		    else{gsl_root_fsolver_set(solver, &Func, lower_h, upper_h);}
		    //find root
		    do{
		    	iter++;
		    	status = gsl_root_fsolver_iterate(solver);
		    	guess = gsl_root_fsolver_root(solver);
		    	h_low = gsl_root_fsolver_x_lower(solver);
		    	h_up = gsl_root_fsolver_x_upper(solver);
		    	status = gsl_root_test_interval(h_low, h_up, 0., h_tol);
		        //cout<<h_low<<"\t"<<guess<<"\t"<<h_up<<"\t"<<minfunc_qks(guess, &params)<<endl;
		    }
		    while (status==GSL_CONTINUE && iter<max_iter);
		    gsl_root_fsolver_free(solver);
		    return guess;
		}
        void make_eos_table(double h_begin, double h_end, int n_table, string unit="MeV"){
            double e, p, n;
            double h = h_begin;
            double h_step = (h_end-h_begin)/n_table;
            if (unit=="cgs"){cout<<"h"<<"\t\t"<<"n/fm-3"<<"\t\t"<<"e/(g*cm-3)"<<"\t\t"<<"p/(erg*cm-3)"<<endl;}
            else {cout<<"h"<<"\t\t"<<"n/fm-3"<<"\t\t"<<"e/(MeV*fm-3)"<<"\t\t"<<"p/(MeV*fm-3)"<<endl;}
            for (int i=0; i<n_table+1; i++){
                n = n_CFL(exp(h)*m_br/3., params_eos[0], params_eos[1], params_eos[2]);
                p = ph_qks_org(h);
                e = eh_qks_org(h, p);
                if (unit=="cgs"){cout<<h<<"\t\t"<<n*MeV3_to_ifm3<<"\t\t"<<e*MeV3_to_ifm3*MeVEifm3_to_gEicm3<<"\t\t"<<p*MeV3_to_ifm3*MeVEifm3_to_ergEicm3<<endl;}
                else {cout<<h<<"\t\t"<<n*MeV3_to_ifm3<<"\t\t"<<e*MeV3_to_ifm3<<"\t\t"<<p*MeV3_to_ifm3<<endl;}
                h += h_step;
            }
        }
};


//----------------------------------------------------------------------------------------------
//                                     Adaptive Piecewise Model
//----------------------------------------------------------------------------------------------


class EoS_adp_pw: public EoS_hybrid{
    public:
    	state_type e_borders, rho_borders;
    	using EoS_hybrid::EoS_hybrid;
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
	        //transfered ppar: rho_borders_ratios*N+p_ratios*N
    		params_eos.clear(); e_borders.clear(); rho_borders.clear();
    		eos_table_max_h = h_max;
	        int len_border = extra_par;
	        for (int i=0; i<len_border; i++){
	            if (i==0){
	                rho_borders.push_back(rho_0*ppar[i]);
	                params_eos.push_back(p_0*ppar[len_border+i]);
	            }
	            else{
	                rho_borders.push_back(rho_borders[i-1]*ppar[i]);
	                params_eos.push_back(params_eos[i-1]*ppar[len_border+i]);
	            }
	        }
    		try{
    		    e_borders = cal_eborders(params_eos, rho_borders);
    	        cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	return false;
    	    }
    	    return true;
    	}
		/** @brief Calculate gamma in the adaptive piecewise model. */
		double gamma_interp(double h, double p, double e){
		    auto upper = std::lower_bound(e_borders.begin()+1, e_borders.end(), e);
		    int idx = std::distance(e_borders.begin(), upper)-1;
		    if (upper==e_borders.end()) {idx -= 1;}
		    if(idx==0) return log(params_eos[0]/p_0)/log(rho_borders[0]/rho_0);
		    else return log(params_eos[idx]/params_eos[idx-1])/log(rho_borders[idx]/rho_borders[idx-1]);
		}
};


//----------------------------------------------------------------------------------------------
//                                     Constant CS Model
//----------------------------------------------------------------------------------------------

// it seems mathematically unstable to join the low density eos with css, which causes a sharp jump in gamma
// physically they are different EOSs if starting from different points
class EoS_cons_cs: public EoS_hybrid{
    public:
    	using EoS_hybrid::EoS_hybrid;
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		//expected ppar: cs
    		params_eos.clear(); 
    		eos_table_max_h = h_max;
	        params_eos.push_back(ppar[0]);
    		try{
    	        cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	return false;
    	    }
	        return true;
	    }
		/** @brief Calculate gamma in the constant speed of sound model.
		  * @attention Define v^2=dp/drho, and p=K*rho^{gamma}.*/
		double gamma_interp(double h, double p, double e){
		    return (e+p)*pow(params_eos[0], 2)/p;
		}
};


//----------------------------------------------------------------------------------------------
//                                  piece_spec_phtr_css model
//----------------------------------------------------------------------------------------------


class EoS_pw_sp_pt_css: public EoS_hybrid{
    public:
	    state_type e_borders, rho_borders;
    	using EoS_hybrid::EoS_hybrid;
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
	        //transfered ppar: p1, g0, g1, g2, rho_tr, deltar_tr, gamma_pt, cs
	        //expected params_eos: rhob1, gamma1, rhob2, h0, g0, g1, g2, rhob3, gamma_pt, cs
	        params_eos.clear(); e_borders.clear(); rho_borders.clear();
	        double rhob_1 = rho_sat, pb_1 = ppar[0];
	        double eb_1, h0, gamma1;
	        eos_table_max_h = h_max;
	        gamma1 = log(pb_1/p_0)/log(rhob_1/rho_0);
	        params_eos.push_back(rhob_1); params_eos.push_back(gamma1); params_eos.push_back(ppar[4]);
	        state_type p_borders{pb_1};
	        rho_borders.push_back(rhob_1);
	        e_borders = cal_eborders(p_borders, rho_borders);
	        eb_1 = e_borders[1];
	        h0 = log((eb_1+pb_1)/rhob_1);
	        params_eos.push_back(h0);
	        params_eos.push_back(ppar[1]); params_eos.push_back(ppar[2]); params_eos.push_back(ppar[3]);
	        params_eos.push_back(ppar[4]+ppar[5]); params_eos.push_back(ppar[6]); params_eos.push_back(ppar[7]);
	        if (verbose){
	        	cout<<"h_0: "<<h_0<<", p_0: "<<p_0<<", rho_0: "<<rho_0<<", e_0: "<<e_0<<endl;
	            cout<<"EoS param:"<<endl<<endl;
	            for (int j=0; j<params_eos.size(); j++) cout<<params_eos[j]<<"  ";
	            cout<<endl<<endl;
	        }
    		try{
    	        cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	return false;
    	    }
	        // check eos constraints
	        if (extra_par==1){
                if (not check_p_at_rho185()){
                    if (verbose) cout<<"Error in update eos: p="<<eos_table_function_rho_base[1](1.85*rho_sat)*p_trans<<" at rho=1.85*rho_sat not allowed!"<<endl;
                    return false;
                }
            }
            if (extra_par==2){
                double start_h = eos_table_function_rho_base[0]((1.+1e-10)*rho_sat);
                if (not check_gamma(start_h, 5.)){
                    if (verbose) cout<<"Error in update eos: gamma is not allowed in the 1-5 rho_sat range!"<<endl;
                    return false;
                }
            }
            return true;
    	}
		/** @brief Calculate gamma in the piece_spec_phtr_css model. */
		double gamma_interp(double h, double p, double e){
		//expected params_eos: rhob1, gamma1, rhob2, h0, g0, g1, g2, rhob3, gamma_pt, cs
		    double rho = (e+p)/exp(h), Gamma, loggamma = 0.;
		    if (rho<params_eos[0]) Gamma = params_eos[1];
		    else if (rho<params_eos[2]) {
		        for (int j = 0; j<3; j++){ loggamma += params_eos[4+j]*pow(log(h/params_eos[3]), j); }
		        Gamma = (e+p)/((exp(loggamma)+1.)*p);
		    }
		    else if (rho<params_eos[7]) Gamma = params_eos[8];
		    else Gamma = (e+p)*pow(params_eos[9], 2)/p;
		    return Gamma;
		}
		/** @brief check whether the pressure at the rho=1.85\rho_sat satisfies constraint from 10.1103/PhysRevD.80.103003 */
		bool check_p_at_rho185(){ return eos_table_function_rho_base[1](1.85*rho_sat)*p_trans>1.21e34;}
		/** @brief special need of the Gamma in the piece_spec_phtr_css method */
		bool check_gamma(double start_h=h_0, double max_check_rho=5.){
		    double h = start_h;
		    double rho, p, e, Gamma;
		    double p_start = eos_table_function_rho_base[1](1.*rho_sat), e_start = eos_table_function_rho_base[2](1.*rho_sat); 
		    double gamma_start = (e_start+p_start)/(p_start*(1.+exp(params_eos[4])));
		    bool valid = ((gamma_start>1.4) and (gamma_start<10.));
		    if (verbose and (not valid)) cout<<"Error in change_pars: Gamma is not allowed at the 1rho_sat point!"<<endl;
		    while(valid){
		        p = eos_table_function_h_base[0](h);
		        e = eos_table_function_h_base[1](h);
		        rho = eos_table_function_h_base[2](h);
		        if (rho>=max_check_rho*rho_sat) break;
		        else if (rho<=rho_sat) {}
		        else {
		            Gamma = gamma(h, p, e);
		            valid &= ((Gamma>1.4) and (Gamma<10.));
		        }
		        h += 0.01;
		    }
		    return valid;
	    }
};


//----------------------------------------------------------------------------------------------
//                               Physical_spec_phtr_css model
//----------------------------------------------------------------------------------------------


class EoS_ph_sp_pt_css: public EoS_hybrid{
    public:
	    state_type e_borders, rho_borders;
    	using EoS_hybrid::EoS_hybrid;
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		//transfered ppar: p1, href, g0, g1, g2, g3, rho_tr, deltar_tr, gamma_pt, cs
	        //expected params_eos: rhob1, gamma1, rhob2, href, g0, g1, g2, g3, rhob3, gamma_pt, cs
	        params_eos.clear(); e_borders.clear(); rho_borders.clear();
	        int p_href = 1;//place of href in parameters
	        double rhob_1 = rho_sat, pb_1 = ppar[0];
	        double gamma1;
	        eos_table_max_h = h_max;
	        gamma1 = log(pb_1/p_0)/log(rhob_1/rho_0);
	        params_eos.push_back(rhob_1); params_eos.push_back(gamma1); params_eos.push_back(ppar[p_href+5]);
	        state_type p_borders{pb_1};
	        rho_borders.push_back(rhob_1);
	        e_borders = cal_eborders(p_borders, rho_borders);
	        params_eos.push_back(ppar[p_href]);//href
	        params_eos.push_back(ppar[p_href+1]); params_eos.push_back(ppar[p_href+2]); params_eos.push_back(ppar[p_href+3]); params_eos.push_back(ppar[p_href+4]);//gamma
	        params_eos.push_back(ppar[p_href+5]+ppar[p_href+6]); params_eos.push_back(ppar[p_href+7]); params_eos.push_back(ppar[p_href+8]);//rhob3, gamma_pt, cs
	        if (verbose){
	            cout<<"eos param:"<<endl<<endl;
	            for (int j=0; j<params_eos.size(); j++) cout<<params_eos[j]<<"  ";
	            cout<<endl<<endl;
            }
    		try{
    	        cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	return false;
    	    }
            return true;
    	}
		/** @brief Calculate gamma in the physical_spec_phtr_css model. */
		double gamma_interp(double h, double p, double e){
		//expected params_eos: rhob1, gamma1, rhob2, href, g0, g1, g2, g3, rhob3, gamma_pt, cs
		    double rho = (e+p)/exp(h), Gamma, loggamma = 0.;
		    if (rho<params_eos[0]) Gamma = params_eos[1];
		    else if (rho<params_eos[2]) {
		        for (int j = 0; j<4; j++){ loggamma += params_eos[4+j]*pow(log(h/params_eos[3]), j); }
		        Gamma = (e+p)/((exp(loggamma)+1.)*p);
		    }
		    else if (rho<params_eos[8]) Gamma = params_eos[9];
		    else Gamma = (e+p)*pow(params_eos[10], 2)/p;
		    return Gamma;
		}
};


//----------------------------------------------------------------------------------------------
//                                   Mu-CS Parameterization Model
//----------------------------------------------------------------------------------------------


class EoS_param_mu_cs_deprecated: public EoS_hybrid{
    public:
    	using EoS_hybrid::EoS_hybrid;
    	state_type h_borders;
    	double eos_jump_tolerance; // extra_par/100, it means extra_par%%
    	double p_pqcd_N, e_pqcd_N, rho_pqcd_N, cssq_pqcd_N; // pressure, energy density and mass density at where pQCD start
    	double gammap, qcd_x; // gamma of piecewise, x of pqcd
    	double h1, cssq1; // properties of the first piecewise segment
    	int n_seg; // N segments of the speed of sound parameterization
    	vector< pchip<state_type> > eos_table_function_h_base_near_pqcd; // reachable only after calling update_eos()
    	state_type eos_table_h_near_pqcd, eos_table_p_near_pqcd, eos_table_e_near_pqcd, eos_table_rho_near_pqcd; // reachable only after calling update_eos()
    	
    	/** @brief Calculate pressure using h only. */
		double ph(double h){
			if (h>h_borders[n_seg-2]) return eos_table_function_h_base_near_pqcd[0](h);
			else if (h>h_border) return eos_table_function_h_base[0](h);
			else return ph_tabulated(h);
		}

		/** @brief Calculate energy density using h only. */
		double eh(double h){
			if (h>h_borders[n_seg-2]) return eos_table_function_h_base_near_pqcd[1](h);
			else if (h>h_border) return eos_table_function_h_base[1](h);
			else return eh_tabulated(h);
		}

		/** @brief Calculate mass density using h only. */
		double rhoh(double h, double p, double e){
			if (h>h_borders[n_seg-2]) return eos_table_function_h_base_near_pqcd[2](h);
			else if (h>h_border) return eos_table_function_h_base[2](h);
			else return rhoh_tabulated(h, p, e);
		}

    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		//transfered ppar: [N*h_borders, (N-1)*cssq, gammap, x]; len_par: N; extra_par: percentage delta_difference; an extra function will be generated if h_max>h_{N-1}.
	        //expected params_eos: (N-1)*cssq
	        bool ret_val = true;
	        params_eos.clear(); h_borders.clear();
	        for (int i=0; i<len_par; i++) h_borders.push_back(ppar[i]);
	        for (int i=0; i<len_par-1; i++) params_eos.push_back(ppar[i+len_par]);
	        gammap = ppar[2*len_par-1]; qcd_x = ppar[2*len_par];
	        eos_jump_tolerance = extra_par/100.; n_seg = len_par;
	        double rho1 = 1.1*rho_sat;
	        double p1 = p_0*pow((rho1/rho_0), gammap);
	        double e1 = rho1/rho_0*e_0+(p1-rho1/rho_0*p_0)/(gammap-1.);
	        double c1 = log(p1/p_0)/log(e1/e_0);
	        double mu_N = exp(h_borders[len_par-1])*m_neutron/1000; // GeV
	        cssq1 = p1*gammap/(e1+p1);
	        h1 = h_0+(c1/(c1-1))*log((e_0*(e1+p1))/(e1*(e_0+p_0)));
            //cout<<h1<<endl;
	        eos_table_max_h = h_max;
	        // update properties at mu_N
	        p_pqcd_N = pressure_renorm(mu_N, qcd_x)*(MeV4_to_dmls*1.e12);
	        e_pqcd_N = energy_density_renorm(mu_N, qcd_x)*(MeV4_to_dmls*1.e12);
	        rho_pqcd_N = number_density_renorm(mu_N, qcd_x)*(m_neutron/1000.)*(MeV4_to_dmls*1.e12);
	        cssq_pqcd_N = cssq_renorm(mu_N, qcd_x);
    		try{
    			if (h_max<h_borders[len_par-2]) cal_eos_table(h_max, init_function_type);
                else {
                	cal_eos_table(h_borders[len_par-2], init_function_type);
                	cal_eos_table_near_pqcd(h_max, init_function_type);// only accept init_function_type=='100'
                	//if (not eos_valid()) {
                	//	if (verbose) cout<<"Failed in update eos: pQCD constraint not fulfilled!"<<endl;
                	//	ret_val = false;}
                }
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	ret_val = false;
    	    }
	        if (verbose) {
	        	cout<<"Border of the piecewise and parameterize speed of sound, h_1: "<<h1<<", cssq1: "<<cssq1<<endl;
	        	cout<<"mu_N/GeV: "<<mu_N<<", p_pqcd_N/(dmls): "<<p_pqcd_N<<", e_pqcd_N/(dmls): "<<e_pqcd_N<<", rho_pqcd_N/(dmls)"<<rho_pqcd_N<<endl;
	        }
	        return ret_val;
	    }

	    double get_loglike_consistency_with_pqcd(){
	    	double delta_p, delta_e, delta_rho;
	    	delta_p = (eos_table_p.back()-eos_table_p_near_pqcd[0])/eos_table_p.back()*100;
	    	delta_e = (eos_table_e.back()-eos_table_e_near_pqcd[0])/eos_table_e.back()*100;
	    	delta_rho = (eos_table_rho.back()-eos_table_rho_near_pqcd[0])/eos_table_rho.back()*100;
	    	if (verbose) cout<<"delta_p: "<<delta_p<<"%, "<<"delta_e: "<<delta_e<<"%, "<<"delta_rho: "<<delta_rho<<"%, standard: "<<eos_jump_tolerance<<"%."<<endl;
	    	//cout<<delta_p<<"  "<<delta_e<<"  "<<delta_rho<<endl;
	    	return -(pow(delta_p, 2)+pow(delta_e, 2)+pow(delta_rho, 2))/(2*pow(eos_jump_tolerance, 2));
	    }

	    double eos_messenger(int ask_for_information_type){ return get_loglike_consistency_with_pqcd(); }

		/** @brief Calculate gamma in the speed of sound model  */
		double gamma_interp(double h, double p, double e){
		//expected params_eos: gamma1, h_i(N_seg), cssq_i(N_seg); expexted auxiliary_variables: h and cssq at the end of the piecewise
		    double Gamma;
		    if (h<=h1) Gamma = gammap;
		    else Gamma = (e+p)/p*interpolate_cssq(h, h_borders, params_eos);
		    return Gamma;
		}

	    double interpolate_cssq(double h, state_type hb, state_type pars){
		    double cssq;
		    auto upper = std::lower_bound(hb.begin(), hb.end(), h);
		    int idx = std::distance(hb.begin(), upper);
		    if (idx==n_seg) cssq = cssq_renorm(exp(h)*m_neutron/1000., qcd_x);
		    else if (idx==(n_seg-1)) cssq = ((exp(hb[n_seg-1])-exp(h))*pars[n_seg-2]+(exp(h)-exp(hb[n_seg-2]))*cssq_pqcd_N)/(exp(hb[n_seg-1])-exp(hb[n_seg-2]));
		    else if (idx==0) cssq = ((exp(hb[0])-exp(h))*cssq1+(exp(h)-exp(h1))*pars[0])/(exp(hb[0])-exp(h1));
		    else cssq = ((exp(hb[idx])-exp(h))*pars[idx-1]+(exp(h)-exp(hb[idx-1]))*pars[idx])/(exp(hb[idx])-exp(hb[idx-1]));
		    return cssq;
		}

		double pressure_renorm(double mu, double x){ // see DOI: 10.1088/2041-8205/781/2/L25
			return pow(mu, 4)/(108*pow(M_PI, 2))*(0.9008-(0.5034*pow(x, -0.3553))/(mu-1.452*pow(x, -0.9101)));
		}

		double number_density_renorm(double mu, double x){
			return 4*pressure_renorm(mu, x)/mu+((pow(mu, 4)/(108*pow(M_PI, 2)))*(0.5034*pow(x, -0.3553)))/pow(mu-1.452*pow(x, -0.9101), 2);
		}

		double energy_density_renorm(double mu, double x){
			return 3*pressure_renorm(mu, x)+mu*((pow(mu, 4)/(108*pow(M_PI, 2)))*(0.5034*pow(x, -0.3553)))/pow(mu-1.452*pow(x, -0.9101), 2);
		}

		double cssq_renorm(double mu, double x){
			return 1./(8-20*pressure_renorm(mu, x)/(mu*number_density_renorm(mu, x))-(2*mu/number_density_renorm(mu, x))*((pow(mu, 4)/(108*pow(M_PI, 2)))*(0.5034*pow(x, -0.3553)))/pow(mu-1.452*pow(x, -0.9101), 3));
		}

	    /** @brief calculate eos tables as list of interpolating functions and store them */
		void cal_eos_table_near_pqcd(double h_max, const char init_function_type[]){
			int intg_steps; // integration steps
		    state_type eos_tb{p_pqcd_N, e_pqcd_N, rho_pqcd_N};
		    state_type eos_table_h_reverse;
		    vector<state_type> x_o;
		    eos_table_h_near_pqcd.clear(); eos_table_p_near_pqcd.clear(); eos_table_e_near_pqcd.clear(); eos_table_rho_near_pqcd.clear();
		    try {
		        integrate_adaptive(controlled_stepper_cal_eos, cal_eos, eos_tb, h_borders[n_seg-1], h_borders[n_seg-2], -sg_step, push_back_state_and_time(x_o, eos_table_h_reverse));
		        intg_steps = eos_table_h_reverse.size();
		        for (int i=1; i<=intg_steps; i++){
		        	eos_table_p_near_pqcd.push_back(x_o[intg_steps-i][0]);
		        	eos_table_e_near_pqcd.push_back(x_o[intg_steps-i][1]);
		        	eos_table_rho_near_pqcd.push_back(x_o[intg_steps-i][2]);
		        	eos_table_h_near_pqcd.push_back(eos_table_h_reverse[intg_steps-i]);
		        }
		        if (h_max>h_borders[n_seg-1]){
					const int segment_pqcd = 30;
		        	double increase_h = (h_max-h_borders[n_seg-1])/segment_pqcd, mu;
		            for (int i=1; i<=segment_pqcd; i++) {
                        mu = exp(h_borders[n_seg-1]+increase_h*i)*m_neutron/1000.;
		            	eos_table_h_near_pqcd.push_back(h_borders[n_seg-1]+increase_h*i);
		         		eos_table_p_near_pqcd.push_back(pressure_renorm(mu, qcd_x)*(MeV4_to_dmls*1.e12));
		        	    eos_table_e_near_pqcd.push_back(energy_density_renorm(mu, qcd_x)*(MeV4_to_dmls*1.e12));
		        	    eos_table_rho_near_pqcd.push_back(number_density_renorm(mu, qcd_x)*(m_neutron/1000.)*(MeV4_to_dmls*1.e12));
		            }
		        }
		        if (init_function_type[0]=='1'){
		            eos_table_function_h_base_near_pqcd.clear();
		            state_type h_bk1(eos_table_h_near_pqcd), h_bk2(eos_table_h_near_pqcd), h_bk3(eos_table_h_near_pqcd);
		            state_type p_bk(eos_table_p_near_pqcd), e_bk(eos_table_e_near_pqcd), rho_bk(eos_table_rho_near_pqcd);
		            auto function_h_p = pchip(std::move(h_bk1), std::move(p_bk));
		            auto function_h_e = pchip(std::move(h_bk2), std::move(e_bk));
		            auto function_h_rho = pchip(std::move(h_bk3), std::move(rho_bk));
		            eos_table_function_h_base_near_pqcd.push_back(function_h_p);
		            eos_table_function_h_base_near_pqcd.push_back(function_h_e);
		            eos_table_function_h_base_near_pqcd.push_back(function_h_rho);
		        }
		        else{
		        	cout<<"For h based only!"<<endl;
		        	if (h_max>h_borders[n_seg-2]) cout<<"Error occured, this kind of param method do not support rho/e based interpolation functions at density near pQCD!"<<endl;
		        	// eos higher than h_borders[n_seg-2] must be calculated this way (and h_based only)!
		        }
		        if (vverbose) {cout<<"initiate type of interpolation function : "<<string(init_function_type)<<endl;}
		    }
		    catch (exception & except) {
		        cout<<"error encountered in cal_eos_table: "<<except.what()<<endl;
		        throw;
		    }
		}
};


class EoS_param_mu_cs: public EoS_hybrid{
    public:
    	using EoS_hybrid::EoS_hybrid;
    	state_type h_borders;
    	double cssq_pqcd_N; // speed of sound where pQCD start
    	double gammap, qcd_x; // gamma of piecewise, x of pqcd
    	double h1, cssq1; // properties of the first piecewise segment
    	int n_seg; // N segments of the speed of sound parameterization
    	
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		//transfered ppar: [(N+1)*h_borders, N*cssq, gammap, x]; len_par: N;
	        //expected params_eos: N*cssq
	        bool ret_val = true;
	        params_eos.clear(); h_borders.clear();
	        for (int i=0; i<len_par+1; i++) h_borders.push_back(ppar[i]);
	        for (int i=0; i<len_par; i++) params_eos.push_back(ppar[i+1+len_par]);
	        gammap = ppar[2*len_par+1]; qcd_x = ppar[2*len_par+2]; n_seg = len_par;
	        double rho1 = 1.1*rho_sat;
	        double p1 = p_0*pow((rho1/rho_0), gammap);
	        double e1 = rho1/rho_0*e_0+(p1-rho1/rho_0*p_0)/(gammap-1.);
	        double mu_N = exp(h_borders[len_par])*m_neutron/1000; // GeV
	        cssq1 = p1*gammap/(e1+p1);
	        h1 = log((p1+e1)/rho1);
	        eos_table_max_h = h_max;
	        cssq_pqcd_N = cssq_renorm(mu_N, qcd_x);
    		try{
    			cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	ret_val = false;
    	    }
	        if (verbose) cout<<"Border of the piecewise and parameterize speed of sound, h_1: "<<h1<<", cssq1: "<<cssq1<<endl;
	        return ret_val;
	    }

		/** @brief Calculate gamma in the speed of sound model  */
		double gamma_interp(double h, double p, double e){
		//expected params_eos: gamma1, h_i(N_seg), cssq_i(N_seg); expexted auxiliary_variables: h and cssq at the end of the piecewise
		    double Gamma;
		    if (h<=h1) Gamma = gammap;
		    else Gamma = (e+p)/p*interpolate_cssq(h, h_borders, params_eos);
		    return Gamma;
		}

	    double interpolate_cssq(double h, state_type hb, state_type pars){
		    double cssq;
		    auto upper = std::lower_bound(hb.begin(), hb.end(), h);
		    int idx = std::distance(hb.begin(), upper);
		    if (idx==n_seg+1) cssq = cssq_renorm(exp(h)*m_neutron/1000., qcd_x);
		    else if (idx==(n_seg)) cssq = ((exp(hb[n_seg])-exp(h))*pars[n_seg-1]+(exp(h)-exp(hb[n_seg-1]))*cssq_pqcd_N)/(exp(hb[n_seg])-exp(hb[n_seg-1]));
		    else if (idx==0) cssq = ((exp(hb[0])-exp(h))*cssq1+(exp(h)-exp(h1))*pars[0])/(exp(hb[0])-exp(h1));
		    else cssq = ((exp(hb[idx])-exp(h))*pars[idx-1]+(exp(h)-exp(hb[idx-1]))*pars[idx])/(exp(hb[idx])-exp(hb[idx-1]));
		    return cssq;
		}

		double pressure_renorm(double mu, double x){ // see DOI: 10.1088/2041-8205/781/2/L25
			return pow(mu, 4)/(108*pow(M_PI, 2))*(0.9008-(0.5034*pow(x, -0.3553))/(mu-1.452*pow(x, -0.9101)));
		}

		double number_density_renorm(double mu, double x){
			return 4*pressure_renorm(mu, x)/mu+((pow(mu, 4)/(108*pow(M_PI, 2)))*(0.5034*pow(x, -0.3553)))/pow(mu-1.452*pow(x, -0.9101), 2);
		}

		double energy_density_renorm(double mu, double x){
			return 3*pressure_renorm(mu, x)+mu*((pow(mu, 4)/(108*pow(M_PI, 2)))*(0.5034*pow(x, -0.3553)))/pow(mu-1.452*pow(x, -0.9101), 2);
		}

		double cssq_renorm(double mu, double x){
			return 1./(8-20*pressure_renorm(mu, x)/(mu*number_density_renorm(mu, x))-(2*mu/number_density_renorm(mu, x))*((pow(mu, 4)/(108*pow(M_PI, 2)))*(0.5034*pow(x, -0.3553)))/pow(mu-1.452*pow(x, -0.9101), 3));
		}

};

class EoS_param_mu_cs_PT: public EoS_param_mu_cs{
    public:
    	using EoS_param_mu_cs::EoS_param_mu_cs;
    	double h_pt, h_pt_end=0, cssq_pt, delta_rho_pt;
    	
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		//transfered ppar: [(N+1)*h_borders, N*cssq, gammap, x]; len_par: N;
	        //expected params_eos: N*cssq
	        bool ret_val = true;
	        params_eos.clear(); h_borders.clear();
	        for (int i=0; i<len_par+1; i++) h_borders.push_back(ppar[i]);
	        for (int i=0; i<len_par; i++) params_eos.push_back(ppar[i+1+len_par]);
	        gammap = ppar[2*len_par+1]; qcd_x = ppar[2*len_par+2]; n_seg = len_par;
	        h_pt = log((ppar[2*len_par+3])/(m_neutron/1000.)); delta_rho_pt = ppar[2*len_par+4]*rho_sat; cssq_pt = ppar[2*len_par+5];
	        double rho1 = 1.1*rho_sat;
	        double p1 = p_0*pow((rho1/rho_0), gammap);
	        double e1 = rho1/rho_0*e_0+(p1-rho1/rho_0*p_0)/(gammap-1.);
	        double mu_N = exp(h_borders[len_par])*m_neutron/1000; // GeV
	        double mu_pt = exp(h_pt)*m_neutron/1000; // GeV
	        cssq1 = p1*gammap/(e1+p1);
	        h1 = log((p1+e1)/rho1);
	        eos_table_max_h = h_max;
	        cssq_pqcd_N = cssq_renorm(mu_N, qcd_x);
    		try{
    			cal_eos_table(h_pt, init_function_type);
    			double rho_pt = eos_table_function_h_base[2](h_pt);
    			double delta_mu_pt = (pow(1+delta_rho_pt/rho_pt, cssq_pt)-1)*mu_pt;
    			h_pt_end = log((delta_mu_pt+mu_pt)/(m_neutron/1000.));
    			cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	ret_val = false;
    	    }
	        if (verbose) cout<<"Border of the piecewise and parameterize speed of sound, h_1: "<<h1<<", cssq1: "<<cssq1<<endl;
	        return ret_val;
	    }

        double eos_messenger(int ask_for_information_type){ return h_pt_end; }

		/** @brief Calculate gamma in the speed of sound model  */
		double gamma_interp(double h, double p, double e){
		//expected params_eos: gamma1, h_i(N_seg), cssq_i(N_seg); expexted auxiliary_variables: h and cssq at the end of the piecewise
		    double Gamma;
		    if (h<=h1) Gamma = gammap;
		    else if (h>h_pt and h<=h_pt_end) Gamma = (e+p)/p*cssq_pt;
		    else Gamma = (e+p)/p*interpolate_cssq(h, h_borders, params_eos);
		    return Gamma;
		}

		/** @brief calculate eos tables as list of interpolating functions and store them */
		void cal_eos_table(double h_max, const char init_function_type[]){
		    state_type eos_tb{p_0, e_0, rho_0};
		    vector<state_type> x_o;
		    eos_table_h.clear(); eos_table_p.clear(); eos_table_e.clear(); eos_table_rho.clear();
		    try {
		    	if (h_max<=h_pt){
		    		integrate_adaptive(controlled_stepper_cal_eos, cal_eos, eos_tb, h_0, h_max, sg_step, push_back_state_and_time(x_o, eos_table_h));
		    	}
		    	else if(h_max<=h_pt_end){ 
		    		integrate_adaptive(controlled_stepper_cal_eos, cal_eos, eos_tb, h_0, h_pt, sg_step, push_back_state_and_time(x_o, eos_table_h));
		    		state_type eos_tb1{x_o.back()[0], x_o.back()[1], x_o.back()[2]};
		    		x_o.pop_back(); eos_table_h.pop_back();
		    		integrate_adaptive(controlled_stepper_cal_eos, cal_eos, eos_tb1, h_pt, h_max, sg_step, push_back_state_and_time(x_o, eos_table_h));
		    	}
		    	else {
		    		integrate_adaptive(controlled_stepper_cal_eos, cal_eos, eos_tb, h_0, h_pt, sg_step, push_back_state_and_time(x_o, eos_table_h));
		    		state_type eos_tb1{x_o.back()[0], x_o.back()[1], x_o.back()[2]};
		    		x_o.pop_back(); eos_table_h.pop_back();
		    		integrate_adaptive(controlled_stepper_cal_eos, cal_eos, eos_tb1, h_pt, h_pt_end, sg_step, push_back_state_and_time(x_o, eos_table_h));
		    		state_type eos_tb2{x_o.back()[0], x_o.back()[1], x_o.back()[2]};
		    		x_o.pop_back(); eos_table_h.pop_back();
		    		integrate_adaptive(controlled_stepper_cal_eos, cal_eos, eos_tb2, h_pt_end, h_max, sg_step, push_back_state_and_time(x_o, eos_table_h));
		    	}
		        for (int i=0; i<eos_table_h.size(); i++){ eos_table_p.push_back(x_o[i][0]), eos_table_e.push_back(x_o[i][1]), eos_table_rho.push_back(x_o[i][2]); }
		        if (init_function_type[0]=='1'){
		            eos_table_function_h_base.clear();
		            state_type h_bk1(eos_table_h), h_bk2(eos_table_h), h_bk3(eos_table_h);
		            state_type p_bk(eos_table_p), e_bk(eos_table_e), rho_bk(eos_table_rho);
		            auto function_h_p = pchip(std::move(h_bk1), std::move(p_bk));
		            auto function_h_e = pchip(std::move(h_bk2), std::move(e_bk));
		            auto function_h_rho = pchip(std::move(h_bk3), std::move(rho_bk));
		            eos_table_function_h_base.push_back(function_h_p);
		            eos_table_function_h_base.push_back(function_h_e);
		            eos_table_function_h_base.push_back(function_h_rho);
		        }
		        if (init_function_type[1]=='1'){
		            eos_table_function_e_base.clear();
		            state_type e_bk1(eos_table_e), e_bk2(eos_table_e), e_bk3(eos_table_e);
		            state_type h_bk(eos_table_h), p_bk(eos_table_p), rho_bk(eos_table_rho);
		            auto function_e_h = pchip(std::move(e_bk1), std::move(h_bk));
		            auto function_e_p = pchip(std::move(e_bk2), std::move(p_bk));
		            auto function_e_rho = pchip(std::move(e_bk3), std::move(rho_bk));
		            eos_table_function_e_base.push_back(function_e_h);
		            eos_table_function_e_base.push_back(function_e_p);
		            eos_table_function_e_base.push_back(function_e_rho);
		        }
		        if (init_function_type[2]=='1'){
		            eos_table_function_rho_base.clear();
		            state_type rho_bk1(eos_table_rho), rho_bk2(eos_table_rho), rho_bk3(eos_table_rho);
		            state_type h_bk(eos_table_h), p_bk(eos_table_p), e_bk(eos_table_e);
		            auto function_rho_h = pchip(std::move(rho_bk1), std::move(h_bk));
		            auto function_rho_p = pchip(std::move(rho_bk2), std::move(p_bk));
		            auto function_rho_e = pchip(std::move(rho_bk3), std::move(e_bk));
		            eos_table_function_rho_base.push_back(function_rho_h);
		            eos_table_function_rho_base.push_back(function_rho_p);
		            eos_table_function_rho_base.push_back(function_rho_e);
		        }
		        if (vverbose) {cout<<"initiate type of interpolation function : "<<string(init_function_type)<<endl;}
		    }
		    catch (exception & except) {
		        cout<<"error encountered in cal_eos_table: "<<except.what()<<endl;
		        throw;
		    }
		}
};


class EoS_param_mu_cs_modified: public EoS_hybrid{
    public:
    	using EoS_hybrid::EoS_hybrid;
    	state_type h_borders;
    	double cssq_pqcd_N; // speed of sound where pQCD start
    	double gammap; // gamma of piecewise, x of pqcd
    	double h1, cssq1; // properties of the first piecewise segment
    	int n_seg; // N segments of the speed of sound parameterization
    	
    	bool update_eos(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    		//transfered ppar: [(N+1)*h_borders, N*cssq, gammap, cssq_pqcd_N]; len_par: N;
	        //expected params_eos: N*cssq
	        bool ret_val = true;
	        params_eos.clear(); h_borders.clear();
	        for (int i=0; i<len_par+1; i++) h_borders.push_back(ppar[i]);
	        for (int i=0; i<len_par; i++) params_eos.push_back(ppar[i+1+len_par]);
	        gammap = ppar[2*len_par+1]; n_seg = len_par;
	        double rho1 = 1.1*rho_sat;
	        double p1 = p_0*pow((rho1/rho_0), gammap);
	        double e1 = rho1/rho_0*e_0+(p1-rho1/rho_0*p_0)/(gammap-1.);
	        cssq1 = p1*gammap/(e1+p1);
	        h1 = log((p1+e1)/rho1);
	        eos_table_max_h = h_max;
	        cssq_pqcd_N = ppar[2*len_par+2];
    		try{
    			cal_eos_table(h_max, init_function_type);
    	    }
    	    catch (exception &except){
    	    	cout<<"Failed in update eos: "<<except.what()<<endl;
    	    	ret_val = false;
    	    }
	        if (verbose) cout<<"Border of the piecewise and parameterize speed of sound, h_1: "<<h1<<", cssq1: "<<cssq1<<endl;
	        return ret_val;
	    }

		/** @brief Calculate gamma in the speed of sound model  */
		double gamma_interp(double h, double p, double e){
		//expected params_eos: gamma1, h_i(N_seg), cssq_i(N_seg); expexted auxiliary_variables: h and cssq at the end of the piecewise
		    double Gamma;
		    if (h<=h1) Gamma = gammap;
		    else Gamma = (e+p)/p*interpolate_cssq(h, h_borders, params_eos);
		    return Gamma;
		}

	    double interpolate_cssq(double h, state_type hb, state_type pars){
		    double cssq;
		    auto upper = std::lower_bound(hb.begin(), hb.end(), h);
		    int idx = std::distance(hb.begin(), upper);
		    if (idx==n_seg+1) cssq = cssq_pqcd_N;
		    else if (idx==(n_seg)) cssq = ((exp(hb[n_seg])-exp(h))*pars[n_seg-1]+(exp(h)-exp(hb[n_seg-1]))*cssq_pqcd_N)/(exp(hb[n_seg])-exp(hb[n_seg-1]));
		    else if (idx==0) cssq = ((exp(hb[0])-exp(h))*cssq1+(exp(h)-exp(h1))*pars[0])/(exp(hb[0])-exp(h1));
		    else cssq = ((exp(hb[idx])-exp(h))*pars[idx-1]+(exp(h)-exp(hb[idx-1]))*pars[idx])/(exp(hb[idx])-exp(hb[idx-1]));
		    return cssq;
		}

};


#endif
