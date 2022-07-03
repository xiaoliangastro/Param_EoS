#include"global_ns_prop_utils.hpp"


/* ===================================================================================
 *  Author: Jin-liang Jiang, Shao-peng Tang
 *  Aim: Using adaptive runge_kutta_fehlberg78 method implemented in boost to integrate DE system
 *  Unit system: C = 1, G = 1, M_sun = 1
 *  Code dependent: c++-17 boost-1.73.0, gsl-2.7 or newer version
 *  Note: 1. Init value-x0 will change after calling integrate_adaptive();
          2. std::move will clear the table;
          3. Please remember to check whether to clear a vector before a new function call of push_back();
 *        4. Always use the abs path please.
 * ===================================================================================
 */


int main(int argc, char *argv[])
{
    return 0;
}