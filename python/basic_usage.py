import numpy as np
import matplotlib.pyplot as plt
import load_cpp_func as cfunc
import ctypes


# function of the script
cal_single_ns = 0
find_maximum_mass = 1
cal_eos_properties = 0


# constant & unit
C = 29979245800.                   # cgs
p_trans = 5.55174079257738e38      # Cactus to cgs
rho_trans = 6.17714470405638e+17   # Cactus to cgs
rho_sat = 2.7e+14                  # cgs


# init control parameters
param_method = 8
# interp_only=1, piecewise=2, piecewise_p=3, spectral=4, spectral_causal=5, causal_sp+phase_trans=6
# quark_star=7, adapt_piecewise=8, const_cs=9, piece_spec_pt_css=10, phy_spec_pt_css=11, mu_cs=12
cfunc.my_cpp_struct.min_tov_mass = 1.4     # used in EoS constraint only
cfunc.my_cpp_struct.max_tov_mass = 3.      # used in EoS constraint only
cfunc.my_cpp_struct.check_causal = True    # check whether the EoS is causal at the maximum mass, if not, cool down to a causal point and the maximum mass is re-defined there
cfunc.my_cpp_struct.param_method = param_method 
cfunc.my_cpp_struct.verbose_level = 1 # set to zero to print nothing
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
max_eos_h = 1.2                    # specify the maximum enthalpy of the EoS
init_function_type = b'111'        # use '100' default to save time, if want to find_eos_properties, change it accroding to the rules in the README.md
if param_method==1:
    extra_par = 0
    _, pp, ee = np.loadtxt("../eos_tables/EOS_TB_Ozel/H4.dat").T/rho_trans
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
else:
    print("Unknown parameterization method: {}".format(param_method))


############################################################################
#                     Usage after EoS initialization                       #
############################################################################


if find_maximum_mass:
    import time
    cfunc.check_mmax_gd(cfunc.phc_max, cfunc.pm_max, 0.3, True, 0.1, True)
    hc_max, m_max = cfunc.get_fp_val(cfunc.phc_max), cfunc.get_fp_val(cfunc.pm_max)
    print("Find maximum mass: {}, correspond central enthalpy: {}".format(hc_max, m_max))
    mass_wanted = 2.1
    cfunc.find_closest_global_property_with_maxm_known(1.2, 0., hc_max, cfunc.phc_closest, cfunc.punknown_closest, 1, False)
    print("Find tidal: {} of mass: {}".format(cfunc.get_fp_val(cfunc.phc_closest), mass_wanted))
    mass_wanted1, mass_wanted2 = 1.1, 1.4
    cfunc.get_unknowns_from_knowns(mass_wanted1, mass_wanted2, cfunc.unknown_find1, cfunc.unknown_find2, cfunc.phc_max, 1)
    print("Find tidal1: {} of mass1: {}; tidal2: {} of mass2: {}".format(\
        cfunc.get_fp_val(cfunc.unknown_find1), mass_wanted1, cfunc.get_fp_val(cfunc.unknown_find2), mass_wanted2))
    # plot the M-R relation of this EoS
    num_points = 50
    hcs = np.linspace(0.1, hc_max, num_points)
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
    h_center = 0.3001
    # get the internal structure of neutron star
    cfunc.my_cpp_struct.cal_internal_structure = True
    cfunc.init_control_params(cfunc.my_cpp_struct)
    print("M: {}/M_sun, R: {}/km, Lambda: {}".format(*cfunc.get_mrl(h_center).contents[:]))
    hss = np.linspace(0.01, 0.3, 100)
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
    plt.plot(rs_interp, ms_interp, label='method interpolation', lw=1.5)
    plt.scatter(rs_interg, ms_interg, label='method intergration', s=2, color='r')
    plt.xlabel(r"$R/\rm km$")
    plt.ylabel(r"$M/M_{\odot}$")
    plt.legend()
    plt.show()


if cal_eos_properties:
    # find eos properties (h, p, e, rho, gamma, vs) using (h / e / rho)
    find_type = 3 # 1->h, 2->e, 3->rho, please also update the EoS using the corresponding init_function_type
    # corespondance of find_type and init_function_type: 1->'100', 2->'010', 3->'001'
    rhos = np.linspace(0.5, 7.4, 100)*rho_sat/rho_trans
    print("Use property  \t\t find: h \t\t p/(dyn/cm^2) \t\t e/(erg/cm^3) \t\t rho/(g/cm^3) \t\t gamma \t\t\t cs^2/C^2")
    for rho in rhos:
        # please input the cactus unit to function find_eos_properties, the output should also be cactus unit
        properties = cfunc.find_eos_properties(rho, 3).contents[:]
        h = properties[0]; p = properties[1]*p_trans
        e = properties[2]*rho_trans*C**2; rho = properties[3]*rho_trans
        gamma = properties[4]; cs_square = properties[5]
        print("\t".join([str(it) for it in [rho, h, p, e, rho, gamma, cs_square]]))
