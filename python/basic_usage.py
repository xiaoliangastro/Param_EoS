import numpy as np
import matplotlib.pyplot as plt
import load_cpp_func as cfunc
import ctypes


# function of the script
cal_single_ns = 0
find_maximum_mass = 0
cal_eos_properties = 1


# constant & unit
C = 29979245800.                   # cgs
p_trans = 5.55174079257738e38      # Cactus to cgs
rho_trans = 6.17714470405638e+17   # Cactus to cgs
rho_sat = 2.7e+14                  # cgs
m_neutron = 939.56542052           # mass of neutron in MeV
e_trans_mevifm3_cactus = 2.88630841e-06    # MeV/fm^3 to cactus


# init control parameters
param_method = 1
# interp_only=1, piecewise=2, piecewise_p=3, spectral=4, spectral_causal=5, causal_sp+phase_trans=6
# quark_star=7, adapt_piecewise=8, const_cs=9, piece_spec_pt_css=10, phy_spec_pt_css=11, mu_cs=12
cfunc.my_cpp_struct.min_tov_mass = 0.     # used in EoS constraint only
cfunc.my_cpp_struct.max_tov_mass = 3.      # used in EoS constraint only
cfunc.my_cpp_struct.check_causal = True    # check whether the EoS is causal at the maximum mass, if not, cool down to a causal point and the maximum mass is re-defined there
cfunc.my_cpp_struct.param_method = param_method 
cfunc.my_cpp_struct.verbose_level = 1 # set to zero to print nothing
cfunc.my_cpp_struct.consid_const_inter_step = False # use constant integration step or not
cfunc.my_cpp_struct.const_inter_step = 0.#1e-5 # integration step
cfunc.init_control_params(cfunc.my_cpp_struct)


############################################################################
#          Initiate EoS in different parameterization methods              #
############################################################################


### Initiate EoS at low density
pp, ee, _, _ = np.loadtxt("../eos_tables/std_ebase_lowdense_eos.txt").T
array_eos = ctypes.c_double*len(pp)
cfunc.create_low_density_eos_with_ep_table(array_eos(*ee), array_eos(*pp), len(pp))

### update EoS at high density
# default parameters
max_eos_h = 1.05                    # specify the maximum enthalpy of the EoS
init_function_type = b'100'        # use '100' default to save time, if want to find_eos_properties, change it accroding to the rules in the README.md
if param_method==1:
    extra_par = 0
    ee, pp = np.loadtxt("/Users/jiangjinliang/Software/code-compose-master/extract_ep/result/ep_HS_DD2.txt", skiprows=1).T*e_trans_mevifm3_cactus
    #da = np.loadtxt("/Users/jiangjinliang/Work/Param_cs/data/eos/DD2+HRG+VQCD_interm_beta.lorene", skiprows=9).T[-2:, :]
    #ee, pp = da[0]/rho_trans, da[1]/p_trans
    #_, pp, ee = np.loadtxt("../eos_tables/EOS_TB_Ozel/H4.dat").T/rho_trans
    ee, unique_index = np.unique(ee, return_index=True)
    pp = pp[unique_index]
    params = np.append(ee, pp)
    len_pars = len(params)
    array_wanted = ctypes.c_double*len_pars
    params = array_wanted(*params)
    cfunc.update_high_density_eos_parameters(params, len_pars, int(len_pars/2), max_eos_h, init_function_type) # this will cover low density one initialized before
elif param_method==2:
    extra_par = 0
    len_pars = 4
    array_wanted = ctypes.c_double*len_pars
    gm1, gm2, gm3, gm4 = 2, 3, 4, 3 # gamma among [ rho_0, 1, 1.85, 3.7, 7.4]*rho_sat, respectively
    params = array_wanted(gm1, gm2, gm3, gm4)
    cfunc.update_high_density_eos_parameters(params, len_pars, extra_par, max_eos_h, init_function_type)
elif param_method==3:
    extra_par = 0
    len_pars = 4
    p_bases = np.array([1.e33, 1.e34, 1.e35, 1.e36])/p_transl
    array_wanted = ctypes.c_double*len_pars
    p1, p2, p3, p4 = p_bases*np.array([6, 3, 4, 3]) # pressures at [1, 1.85, 3.7, 7.4]*rho_sat, respectively
    params = array_wanted(p1, p2, p3, p4)
    cfunc.update_high_density_eos_parameters(params, len_pars, extra_par, max_eos_h, init_function_type)
elif (param_method==4 or param_method==5):
    extra_par = 4 # where to find the h_base
    len_pars = 5
    h_base = 0.03
    array_wanted = ctypes.c_double*len_pars
    g1, g2, g3, g4 = 1.213, -0.15, 0.3, -0.0015 # expantion parameters gamma
    params = array_wanted(g1, g2, g3, g4, h_base)
    cfunc.update_high_density_eos_parameters(params, len_pars, extra_par, max_eos_h, init_function_type)
elif param_method==6:
    len_pars = 8
    extra_par = 5
    array_wanted = ctypes.c_double*len_pars
    eb1, eb2, eb3, eb4, eb5 = np.array([1, 1.2, 1.35, 1.35+0.25, 2.3])*rho_transl
    gm1, gm2, gm3 = 3, 4, 1.07
    params = array_wanted(eb1, eb2, eb3, eb4, eb5, gm1, gm2, gm3)
    cfunc.update_high_density_eos_parameters(params, len_pars, extra_par, max_eos_h, init_function_type)
elif param_method==7:
    len_pars = 3
    extra_par = 0
    array_wanted = ctypes.c_double*len_pars
    a4, beff, egap = 0.5, 130**4, 50
    params = array_wanted(a4, beff, egap)
    cfunc.update_high_density_eos_parameters(params, len_pars, extra_par, max_eos_h, init_function_type)
elif param_method==8:
    extra_par = 5 # segements of pieces
    len_pars = 2*extra_par
    array_wanted = ctypes.c_double*len_pars
    rho_ratios = [1.1, 1.2, 1.05, 1.6, 1.2]
    p_ratios = [1.1, 1.7, 1.3, 2., 1.6]
    params = array_wanted(*(rho_ratios+p_ratios))
    cfunc.update_high_density_eos_parameters(params, len_pars, extra_par, max_eos_h, init_function_type)
elif param_method==9:
    extra_par = 0
    len_pars = 1
    array_wanted = ctypes.c_double*len_pars
    params = array_wanted(np.sqrt(1./3)) # const cs with a user specified curst
    cfunc.update_high_density_eos_parameters(params, len_pars, extra_par, max_eos_h, init_function_type)
elif param_method==10:
    extra_par = 0  # whether to check the eos constraint
    len_pars = 8
    array_wanted = ctypes.c_double*len_pars
    p1, g0, g1, g2, rhotr, drho, gammatr, cs2 = [4.1, 0.53, 1.90, -0.57, 2.33, 0.32, 1.18, 0.60]
    params = array_wanted(p1*p_bases[0], g0, g1, g2, rhotr*rho_transl, drho*rho_transl, gammatr, np.sqrt(cs2))
    cfunc.update_high_density_eos_parameters(params, len_pars, extra_par, max_eos_h, init_function_type)
elif param_method==11:
    extra_par = 0  # whether to check the eos constraint and in which way
    len_pars = 10
    array_wanted = ctypes.c_double*len_pars
    p1, href, g0, g1, g2, g3, rhotr, drho, gammatr, cs2 = [4.1, 0.03, 0.53, 1.90, -0.57, 0.034, 2.33, 0.32, 1.18, 0.60]
    params = array_wanted(p1*p_bases[0], href, g0, g1, g2, g3, rhotr*rho_transl, drho*rho_transl, gammatr, np.sqrt(cs2))
    cfunc.update_high_density_eos_parameters(params, len_pars, extra_par, max_eos_h, init_function_type)
elif param_method==12:
    extra_par = 50  # delta miss-match
    N_seg = 7
    h_borders = np.log(np.linspace(np.exp(0.0816780310142671), np.exp(1.01785), N_seg)).tolist()
    array_wanted = ctypes.c_double*(2*N_seg+1)
    params = array_wanted(*(h_borders+[0.3, 0.4, 0.6, 0.55, 0.35, 0.4]+[2.3, 3.]))
    aa = cfunc.update_high_density_eos_parameters(params, N_seg, extra_par, max_eos_h, init_function_type)
elif param_method==13:
    extra_par = 0
    N_seg = 10
    h_pt = 0.27
    h_borders = np.linspace(0.0816780310142671, 1.0178492742434604, N_seg+2)[1:].tolist() # mu0=1.0GeV
    array_wanted = ctypes.c_double*(2*N_seg+6)
    params = array_wanted(*(h_borders+[0.3, 0.4, 0.6, 0.5, 0.7, 0.4, 0.5, 0.4, 0.3, 0.5]+[2.3, 3.]+[np.exp(h_pt)*m_neutron/1000, 1, 0.01]))
    aa = cfunc.update_high_density_eos_parameters(params, N_seg, extra_par, max_eos_h, init_function_type)
    h_pt_end = cfunc.eos_messenger(0)
    print(h_pt_end-h_pt)
else:
    print("Unknown parameterization method: {}".format(param_method))


############################################################################
#                     Usage after EoS initialization                       #
############################################################################


if find_maximum_mass:
    import time
    model_phase_transition = 0
    if model_phase_transition:
        cfunc.check_mmax_pt_two_branch(cfunc.phc_1, cfunc.phc_2, cfunc.phc_3, cfunc.phc_max, cfunc.pm_1, cfunc.pm_2, cfunc.pm_3, cfunc.pm_max, h_pt, h_pt_end)
        hc_max = cfunc.get_fp_val(cfunc.phc_max)
    else:
        cfunc.check_mmax_gd(cfunc.phc_max, cfunc.pm_max, 0.3, True, 0.1, True)
        hc_max, m_max = cfunc.get_fp_val(cfunc.phc_max), cfunc.get_fp_val(cfunc.pm_max)
        print("Find maximum mass: {}, correspond central enthalpy: {}".format(m_max, hc_max))
        mass_wanted = 1.4
        cfunc.find_closest_global_property_with_maxm_known(mass_wanted, 0., hc_max, cfunc.phc_closest, cfunc.punknown_closest, 1, False)
        print("Find tidal: {} of mass: {}, correspond to h: {}".format(cfunc.get_fp_val(cfunc.punknown_closest), mass_wanted, cfunc.get_fp_val(cfunc.phc_closest)))
        mass_wanted1, mass_wanted2 = 1.1, 1.4
        cfunc.get_unknowns_from_knowns(mass_wanted1, mass_wanted2, cfunc.unknown_find1, cfunc.unknown_find2, cfunc.phc_max, 1)
        print("Find tidal1: {} of mass1: {}; tidal2: {} of mass2: {}".format(\
            cfunc.get_fp_val(cfunc.unknown_find1), mass_wanted1, cfunc.get_fp_val(cfunc.unknown_find2), mass_wanted2))
    # plot the M-R relation of this EoS
    num_points = 120
    hcs = np.linspace(0.1, hc_max+0.5, num_points)
    ms, rs, ls = [], [], []
    t_start = time.time()
    for hc in hcs:
        m,r,l = cfunc.get_mrl(hc).contents[:]
        ms.append(m); rs.append(r); ls.append(l)
    print("Single TOV integration used time: {} s.".format((time.time()-t_start)/num_points))
    plt.scatter(rs, ms, s=1.3)
    plt.xlabel(r"$R/\rm km$")
    plt.ylabel(r"$M/M_{\odot}$")
    plt.show()


if cal_single_ns:
    h_center = 0.20055999999999996
    # get the internal structure of neutron star
    cfunc.my_cpp_struct.cal_internal_structure = True
    cfunc.init_control_params(cfunc.my_cpp_struct)
    mrl = cfunc.get_mrl(h_center).contents[:]
    print("M: {}/M_sun, R: {}/km, Lambda: {}".format(*mrl))
    hss = np.linspace(8.e-3, h_center, 100)
    # method1: interpolation using the integration steps in get_mrl, TOV integration once
    ms_interp, rs_interp = [], []
    for h_surf in hss:
        # function calculate_internal_structure must be called after get_mrl
        mr_interp = cfunc.calculate_internal_structure(h_surf).contents[:]
        ms_interp.append(mr_interp[0]); rs_interp.append(mr_interp[1]);
    # method2: do the whole TOV integration again and again
    ms_interg, rs_interg = [], []
    for h_surf in hss:
        mr_interg = cfunc.get_mr_with_specific_hsurf(h_center, h_surf).contents[:]
        ms_interg.append(mr_interg[0]); rs_interg.append(mr_interg[1]);
    plt.scatter(rs_interg, ms_interg, label='method intergration', s=2, color='r')
    plt.plot(rs_interp, ms_interp, label='method interpolation', lw=1.5)
    plt.xlabel(r"$R/\rm km$")
    plt.ylabel(r"$M/M_{\odot}$")
    plt.legend()
    plt.show()


if cal_eos_properties:
    # find eos properties (h, p, e, rho, gamma, vs) using (h / e / rho)
    find_type = 1 # 1->h, 2->e, 3->rho, please also update the EoS using the corresponding init_function_type
    if find_type==3:
        # corespondance of find_type and init_function_type: 1->'100', 2->'010', 3->'001'
        rhos = np.linspace(0.5, 7.4, 100)*rho_sat/rho_trans
        print("Use property  \t\t find: h \t\t p/(dyn/cm^2) \t\t e/(erg/cm^3) \t\t rho/(g/cm^3) \t\t gamma \t\t\t cs^2/C^2")
        for rho in rhos:
            # please input the cactus unit to function find_eos_properties, the output should also be cactus unit
            properties = cfunc.find_eos_properties(rho, find_type).contents[:]
            h = properties[0]; p = properties[1]*p_trans
            e = properties[2]*rho_trans*C**2; rho = properties[3]*rho_trans
            gamma = properties[4]; cs_square = properties[5]
            print("\t".join([str(it) for it in [rho, h, p, e, rho, gamma, cs_square]]))
    if find_type==1:
        hs = np.linspace(8.e-3, max_eos_h, 50)
        show_x, show_y = [], []
        for h in hs:
            properties = cfunc.find_eos_properties(h, find_type).contents[:]
            show_y.append(properties[2]*rho_trans/rho_sat)
            #show_y.append(properties[-1])
            show_x.append(properties[0])
        properties_max = cfunc.find_eos_properties(0.72056, find_type).contents[:]
        print("{:e}".format(properties_max[3]*rho_trans))
        #properties_pt = cfunc.find_eos_properties(h_pt, find_type).contents[:]
        #properties_pt_end = cfunc.find_eos_properties(h_pt_end, find_type).contents[:]
        #print((properties_pt_end[3]-properties_pt[3])*rho_trans/rho_sat)
        plt.plot(show_x, show_y)
        plt.show()