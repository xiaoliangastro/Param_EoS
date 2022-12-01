# 1. Introduction

## 1.1. Code 
- Aim: light weight and fast tool to integrate the TOV function;
- Motivation: in Bayesian analysis of EoS parameterization, single likelihood evaluation time is vital;
- Method: adaptive step [runge_kutta_fehlberg78](https://www.boost.org/doc/libs/1_73_0/libs/numeric/odeint/doc/html/index.html) method implemented in boost;
- Benefit: 1. Fast: the integration and much of the computational expensive works are done in c++, and be embedded in a [dynamic link library](./eostool/param_eos.so)(DLL). 2. Convenient: a [python API](./python/load_cpp_func/load_cpp_func.py) is also provided for the cpp DLL.

## 1.2 Implement
- EoS at low density: most of the EoS are fixed in low density using a specific EoS nowadays and be interpolated in a way that is introduced in Appendix B of [this paper](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.86.084003);
- EoS at high density: the high density part are described by a few parameters that you want to infer, the EOS is integrated independent of the TOV system, and then stored in a series of interpolation functions;
- TOV: integration system is introduced in [this paper](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.89.064003);
- Unit system: to avoid steep integrations, we use Cactus Unit.

# 2. Dependency

### 2.1. C++
- A g++ compiler accept std c++17;
- Boost 1.73.0 or newer
- GSL 2.7 or newer
### 2.2. python
- ctypes

# 3. Install
```bash 
cd eostool
bash compile.sh dynamic
# if failed, specify where are the boost and the gsl yourself in the file compile.sh 
# please also append 'Your_Dir/python/load_cpp_func' to your envirenment variable PYTHONPATH
```



# 4. Usage & Important API
You should: a. initialize the global control parameters that control the behavior of  the whole script; b. initialize the low density fixed part of the EoS; c. initialize the high density parameterized EoS of a specific type; d. do the calculations below. We also provide a basic [example](./python/basic_usage.py) of usage. More detailed descriptions can also be found  in [documentation](./documentation/html/index.html).

#### 1. get_mrl
```cpp
double* get_mrl(double hc)
```
- **Aim**: get mass, radius and tidal deformability with a given central enthalpy *hc*

#### 2. get_mr_with_specific_hsurf
```cpp
double* get_mr_with_specific_hsurf(double hc, double h_surf)
```
- **Aim**: get the m-r relation of a neutron star with central enthalpy hc
- return: mass (in unit of $M_{\odot}$) and radius (in unit of $\rm km$)

#### 3. check_mmax_gd
```cpp
bool check_mmax_gd(double *hc, double *M_max, double h_start, bool compare_ftype_is_less, double jump_size, bool check_ok)
```
- **Aim**: find TOV mass *M_max* and central enthalpy *hc* related using a gradient method, both of them are set to 0 if failed
- h_start: where to start the finding process
- compare_ftype_is_less: set to true if you want to find the maximum mass, set to false for the minimum (useful only in the case of phase transition where in some cases increase h causes a decrease in mass and you want to know the h in a critical point)
- jump_size: step size of hc, it‘s a good choice to set to 0.1
- check_ok: whether to check the TOV mass allowed in the range [minm_tov, maxm_tov], return false if not

#### 4. find_closest_global_property_with_maxm_known
```cpp
 bool find_closest_global_property_with_maxm_known(double known_aim, double h_i, double h_max, double *h_closest, double *unknown, int get_type, bool use_user_start_point)
```
- **Aim**: find the closest global property with another one after the maximum mass is known
- known_aim: global property already known
- h_i: where to start the finding process
- h_max: central enthalpy related to TOV mass (you can get it using *check_mmax_gd*)
- h_cloest: hc that give global property that very close to known_aim
- unknown: global property you want to find
- get_type: 1-->using mass to find lambda; 2-->using mass to find radius; 3-->using lambda to find mass
- use_user_start_point: set to false if you don't know how to set start h_i

#### 5. get_unknowns_from_knowns
```cpp
bool get_unknowns_from_knowns(double known1, double known2, double *unknown1, double *unknown2, double *h_max, int get_type)
```
- **Aim**:  in a binary system, find corresponding global properties from properties that already known, e.g., find Mass1 and Mass2 using Lambda1 and Lambda2, this is a higher level combination function of *check_mmax_gd* and *find_closest_global_property_with_maxm_known*
- get_type: 1-->using mass to find lambda; 2-->using mass to find radius; 3-->using lambda to find mass

#### 6. find_eos_properties
```cpp
double* find_eos_properties(double known_aim, int find_type)
```
- **Aim**: find eos property y using x; y include (h, p, e, rho, gamma, cs^2), x include (h, e, rho)
- find_type: 1-->x=h, 2-->x=e, 3-->x=rho
- Note: all the input and output are in Cactus unit, and the find_type should be correspond to **init_function_type** in function *update_high_density_eos_parameters* with the coresponding relation: 1-->'100', 2-->'010', 3-->'001'


# 5. Extend the tool yourself

Do not always complain the shit hills inherited from the senior brothers and keep watching&waiting at the same time. You should create your shit hill and bother other novices, too. This can be done in the following simple way.
### 5.1.  Implement your new parameterization method :
- Create a class that publicly inherent the class *EoS_hybrid*, by just re-implement two functions: update_eos and gamma_interp (An example is class *EoS_ph_sp_pt_css*);
- Add new name of your eos class in emum *EoS_type*;
- Create a new branch in switch structure of function *create_low_density_eos_with_ep_table*;
- Change `load_cpp_func.my_cpp_struct.param_method` to your method in the python script.

### 5.2. Add new function to the tool
- Write a function that makes your dream come true;
- Add a function declaration copy to the `extern "C"{}` structure;
- Specify the input and output argtype in the file `load_cpp_func.py`;
- Call it in your python script.

### 5.3. Write your new tools completely different
- Feel free to be inspired by skills of this simple tool and even copy some of the codes.
