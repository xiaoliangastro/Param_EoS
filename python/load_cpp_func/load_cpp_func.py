import ctypes
import os


fthis_file = os.path.realpath(__file__)
fbase = fthis_file.split('python')[0]
fdll = fbase+"eostool/param_eos.so"
eos_name = fbase+"eos_tables/std_ebase_lowdense_eos.txt"
lib = ctypes.CDLL(fdll)
cb = ctypes.c_bool
cf = ctypes.c_double
ci = ctypes.c_int
ccp = ctypes.c_char_p
#cstr = ctypes.c_wchar_p
cc = ctypes.c_char
cfp = ctypes.POINTER(cf)
class My_cpp_structure(ctypes.Structure):
    _fields_ = [('min_tov_mass', cf), ('max_tov_mass', cf), ('check_causal', cb), ('param_method', ci), ('verbose_level', ci), \
                ('const_inter_step', cf), ('consid_const_inter_step', cb), ('cal_internal_structure', cb)] # the order should be exactly the same as which in the cpp file


# eos related C function type (defined in eos_utils.hpp)
init_control_params = lib.init_control_params
init_control_params.restype = None
init_control_params.argtypes = (My_cpp_structure, )
create_low_density_eos_with_ep_table = lib.create_low_density_eos_with_ep_table
create_low_density_eos_with_ep_table.restype = None
create_low_density_eos_with_ep_table.argtypes = (cfp, cfp, ci)
update_high_density_eos_parameters = lib.update_high_density_eos_parameters
update_high_density_eos_parameters.restype = cb
update_high_density_eos_parameters.argtypes = (cfp, ci, ci, cf, ccp)
eos_messenger = lib.eos_messenger
eos_messenger.restype = cf
eos_messenger.argtypes = (ci, )
find_eos_properties = lib.find_eos_properties
find_eos_properties.restype = ctypes.POINTER(cf*6)
find_eos_properties.argtypes = (cf, ci)
make_eos_table = lib.make_eos_table
make_eos_table.restype = None
make_eos_table.argtypes = (cf, cf, cf, ci, ccp, ccp, ci)
interp_pe_likelihood = lib.interp_pe_likelihood
interp_pe_likelihood.restype = cf
interp_pe_likelihood.argtypes = (cfp, cfp, ci)


# global property related C function type (defined in global_ns_prop_utils.hpp)
get_mrl = lib.get_mrl
get_mrl.restype = ctypes.POINTER(cf*3)
get_mrl.argtypes = (cf, )
calculate_internal_structure = lib.calculate_internal_structure
calculate_internal_structure.restype = ctypes.POINTER(cf*2)
calculate_internal_structure.argtypes = (cf, )
get_mr_with_specific_hsurf = lib.get_mr_with_specific_hsurf
get_mr_with_specific_hsurf.restype = ctypes.POINTER(cf*2)
get_mr_with_specific_hsurf.argtypes = (cf, cf)
check_mmax_gd = lib.check_mmax_gd
check_mmax_gd.restype = cb
check_mmax_gd.argtypes = (cfp, cfp, cf, cb, cf, cb)
check_mmax = lib.check_mmax
check_mmax.restype = cb
check_mmax.argtypes = (cfp, cfp, cf, cb, cf, cb)
check_mmax_pt_two_branch = lib.check_mmax_pt_two_branch
check_mmax_pt_two_branch.restype = cb
check_mmax_pt_two_branch.argtypes = (cfp, cfp, cfp, cfp, cfp, cfp, cfp, cfp, cf, cf)
find_closest_global_property_with_maxm_known = lib.find_closest_global_property_with_maxm_known
find_closest_global_property_with_maxm_known.restype = cb
find_closest_global_property_with_maxm_known.argtypes = (cf, cf, cf, cfp, cfp, ci, cb)
get_unknowns_from_knowns = lib.get_unknowns_from_knowns
get_unknowns_from_knowns.restype = cb
get_unknowns_from_knowns.argtypes = (cf, cf, cfp, cfp, cfp, ci)


# storage variables use frequently
unknown_find1 = cfp(cf(0.0))
unknown_find2 = cfp(cf(0.0))
phc_closest = cfp(cf(0.0))
punknown_closest = cfp(cf(0.0))
pm_1 = cfp(cf(0.0))
pm_2 = cfp(cf(0.0))
pm_3 = cfp(cf(0.0))
phc_1 = cfp(cf(0.0))
phc_2 = cfp(cf(0.0))
phc_3 = cfp(cf(0.0))
pm_max = cfp(cf(0.0))
phc_max = cfp(cf(0.0))
get_fp_val = lambda x: x.contents.value
# initiate control parameter struct
my_cpp_struct = My_cpp_structure()
my_cpp_struct.consid_const_inter_step = False   # Consider constant integration step or adaptive step size.
my_cpp_struct.const_inter_step = 0.             # If consider constant integration step, then please give the wanted step size.
my_cpp_struct.min_tov_mass = 0.0                # Minimum mass allowed astropysically, useful in check_mmax.
my_cpp_struct.max_tov_mass = 100.               # Maximum mass allowed astropysically, useful in check_mmax.
my_cpp_struct.verbose_level = 0                 # Verbose level, 0-3, more details will be printed if verbose_level is larger.
my_cpp_struct.check_causal = False              # Some EoS do not causal in the maximum mass place, then find the hc where cs=1, the maximum mass in this case is defined by M(hc|_{cs=1}).
my_cpp_struct.cal_internal_structure = False    # Whether to calculate the internal structure of the neutron star given hc.
