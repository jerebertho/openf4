/* 
 * Copyright (C) 2015 Antoine Joux, Vanessa Vitse and Titouan Coladon
 * 
 * This file is part of openf4.
 * 
 * openf4 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * openf4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with openf4.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 *  \file benchmark-long.cpp
 *  \example benchmark-long.cpp
 *  \brief Benchmark with integer 64 bits coefficients.
 *  \ingroup benchmark
 *  \author Vanessa VITSE, Antoine JOUX, Titouan COLADON
 */

#include <iostream>
#include <openf4.h>

using namespace F4;
using namespace std;

// Global variable
int F4::VERBOSE=1;
#ifdef USE_OPENMP
int F4::NB_THREAD=min(1/*8*/, omp_get_num_procs());
#else
int F4::NB_THREAD=1;
#endif

// Init element-prime tools
typedef ElementPrime<int64_t> eltType;
int64_t modulo=1073741827LL;


int cp_d3_n6_p6F4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                      CP(3,6,6)                        #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    
    // Init monomial tools
    Monomial::initMonomial(6);
    
    // Create variable name array
    string * vars = new string[6];
    for(int i = 0; i < 6; i++)
      {
	vars[i]='x'+to_string(i+1);
      }
    Monomial::setVariable(vars);

    // Create polynomial array
    vector<Polynomial<eltType>> polCp_D3_N6_P6;
    
    // Fill the polynomial array
    polCp_D3_N6_P6.emplace_back("-7*x1^3+22*x1^2*x2-62*x1*x2^2+6*x2^3-55*x1^2*x3+97*x1*x2*x3+74*x2^2*x3+62*x1*x3^2+44*x2*x3^2+68*x3^3-94*x1^2*x4-73*x1*x2*x4+72*x2^2*x4-82*x1*x3*x4+29*x2*x3*x4-10*x3^2*x4-17*x1*x4^2-61*x2*x4^2+95*x3*x4^2-96*x4^3+87*x1^2*x5-4*x1*x2*x5+37*x2^2*x5+80*x1*x3*x5+98*x2*x3*x5+31*x3^2*x5-75*x1*x4*x5-8*x2*x4*x5+x3*x4*x5+72*x4^2*x5-40*x1*x5^2+11*x2*x5^2-28*x3*x5^2-90*x4*x5^2+53*x5^3-56*x1^2*x6-83*x1*x2*x6-23*x2^2*x6-44*x1*x3*x6-23*x2*x3*x6-51*x3^2*x6-10*x1*x4*x6-29*x2*x4*x6+x3*x4*x6-87*x4^2*x6+42*x1*x5*x6-49*x2*x5*x6+16*x3*x5*x6+43*x4*x5*x6-28*x5^2*x6+23*x1*x6^2+40*x2*x6^2-27*x3*x6^2-91*x4*x6^2+13*x5*x6^2+71*x6^3-10*x1*x2+87*x2^2+71*x1*x3+10*x2*x3+77*x3^2-7*x1*x4+95*x2*x4+55*x3*x4+47*x4^2-50*x1*x5-47*x2*x5+30*x3*x5+92*x4*x5+5*x5^2+75*x1*x6-81*x2*x6-15*x3*x6-88*x4*x6-10*x5*x6+16*x6^2-92*x1+91*x2-59*x3-48*x4-82*x5+83*x6+9");
    polCp_D3_N6_P6.emplace_back("-60*x1^3-83*x1^2*x2+5*x1*x2^2+33*x2^3+98*x1^2*x3+96*x1*x2*x3+10*x2^2*x3+98*x1*x3^2-25*x2*x3^2-89*x3^3-48*x1^2*x4-17*x1*x2*x4+7*x2^2*x4-64*x1*x3*x4-96*x2*x3*x4-77*x3^2*x4-34*x1*x4^2+7*x2*x4^2-42*x3*x4^2-69*x4^3-19*x1^2*x5+25*x1*x2*x5-89*x2^2*x5+64*x1*x3*x5+50*x2*x3*x5+69*x3^2*x5-13*x1*x4*x5-89*x2*x4*x5-33*x3*x4*x5-46*x4^2*x5+71*x1*x5^2-68*x2*x5^2+97*x3*x5^2-34*x4*x5^2-85*x5^3+62*x1^2*x6+91*x1*x2*x6+65*x2^2*x6-90*x1*x3*x6-60*x2*x3*x6+80*x3^2*x6+44*x1*x4*x6-70*x2*x4*x6+21*x3*x4*x6-33*x4^2*x6-47*x1*x5*x6-60*x2*x5*x6+30*x3*x5*x6+40*x4*x5*x6+54*x5^2*x6-53*x1*x6^2+52*x2*x6^2+89*x3*x6^2+x4*x6^2+18*x5*x6^2+91*x6^3+37*x1^2+12*x2^2-60*x1*x3-42*x2*x3+28*x3^2-2*x1*x4+34*x2*x4-35*x3*x4+87*x4^2-39*x1*x5+16*x2*x5-64*x3*x5+77*x4*x5-72*x1*x6-20*x2*x6-16*x3*x6-10*x4*x6+52*x5*x6-22*x6^2-97*x1-4*x2+59*x3-65*x4+36*x5+51*x6-27");
    polCp_D3_N6_P6.emplace_back("50*x1^3+60*x1^2*x2+25*x1*x2^2-67*x2^3-91*x1^2*x3+31*x1*x2*x3+58*x2^2*x3-6*x1*x3^2-57*x2*x3^2+36*x3^3-47*x1^2*x4-27*x1*x2*x4+29*x2^2*x4+80*x1*x3*x4+85*x2*x3*x4-8*x3^2*x4+31*x1*x4^2+27*x2*x4^2-91*x3*x4^2-85*x4^3-97*x1^2*x5+65*x1*x2*x5+37*x2^2*x5-84*x1*x3*x5+80*x2*x3*x5+30*x3^2*x5+73*x1*x4*x5+9*x2*x4*x5-70*x3*x4*x5-44*x4^2*x5-29*x1*x5^2+65*x2*x5^2-21*x3*x5^2+49*x4*x5^2-97*x5^3-2*x1^2*x6+88*x1*x2*x6+5*x2^2*x6+57*x1*x3*x6+90*x2*x3*x6-3*x3^2*x6+95*x1*x4*x6-91*x2*x4*x6+42*x3*x4*x6-31*x4^2*x6+5*x1*x5*x6-12*x2*x5*x6-27*x3*x5*x6-58*x4*x5*x6-14*x5^2*x6-51*x1*x6^2+5*x2*x6^2-22*x3*x6^2+x4*x6^2-96*x5*x6^2+62*x6^3-31*x1^2+10*x1*x2-36*x2^2-49*x1*x3+74*x2*x3-56*x3^2+68*x1*x4+81*x2*x4+9*x3*x4+45*x4^2-26*x1*x5+78*x2*x5-79*x3*x5+49*x4*x5+83*x5^2+88*x1*x6-63*x2*x6-51*x3*x6-95*x4*x6-8*x5*x6+96*x6^2+97*x1-5*x2+16*x3+86*x4-54*x5-51*x6+89");
    polCp_D3_N6_P6.emplace_back("14*x1^3-79*x1^2*x2+57*x1*x2^2-22*x2^3-58*x1^2*x3-35*x1*x2*x3+48*x2^2*x3-71*x1*x3^2+18*x2*x3^2-39*x3^3-95*x1^2*x4+57*x1*x2*x4-82*x2^2*x4-66*x1*x3*x4+16*x2*x3*x4-20*x3^2*x4-68*x1*x4^2-59*x3*x4^2-17*x4^3+61*x1^2*x5+28*x1*x2*x5+46*x2^2*x5-34*x1*x3*x5-22*x2*x3*x5+31*x3^2*x5-15*x1*x4*x5-87*x2*x4*x5-54*x3*x4*x5+83*x4^2*x5+87*x1*x5^2-63*x2*x5^2+63*x3*x5^2-5*x4*x5^2-19*x5^3-2*x1^2*x6+63*x1*x2*x6-49*x2^2*x6+72*x1*x3*x6+38*x2*x3*x6-52*x3^2*x6-32*x1*x4*x6-13*x2*x4*x6-22*x3*x4*x6-65*x4^2*x6-60*x1*x5*x6-31*x2*x5*x6+21*x3*x5*x6+27*x4*x5*x6-62*x5^2*x6-61*x1*x6^2+66*x2*x6^2-21*x3*x6^2+26*x4*x6^2+48*x5*x6^2+96*x6^3+86*x1^2+21*x1*x2-66*x2^2-40*x1*x3-48*x2*x3+37*x3^2+17*x1*x4-87*x2*x4+17*x3*x4-28*x4^2+7*x1*x5+90*x2*x5-12*x3*x5-12*x5^2+45*x1*x6+88*x2*x6-74*x3*x6+13*x4*x6-84*x5*x6+6*x6^2-50*x1+74*x2+55*x3+83*x4-83*x5+76*x6+95");
    polCp_D3_N6_P6.emplace_back("-14*x1^3+61*x1^2*x2-40*x1*x2^2-15*x2^3+83*x1^2*x3+45*x1*x2*x3-67*x2^2*x3+92*x1*x3^2+95*x2*x3^2-55*x3^3-83*x1^2*x4-86*x1*x2*x4-62*x2^2*x4+35*x1*x3*x4+65*x2*x3*x4-34*x3^2*x4+31*x1*x4^2-x2*x4^2+11*x3*x4^2+60*x4^3+53*x1^2*x5-88*x1*x2*x5+70*x2^2*x5-31*x1*x3*x5-10*x2*x3*x5+52*x3^2*x5-52*x1*x4*x5+49*x2*x4*x5-11*x3*x4*x5-90*x4^2*x5-37*x1*x5^2-63*x2*x5^2+6*x3*x5^2-53*x4*x5^2+2*x5^3-6*x1^2*x6-21*x1*x2*x6+93*x2^2*x6-17*x1*x3*x6+11*x2*x3*x6+x3^2*x6+44*x1*x4*x6-27*x2*x4*x6+88*x3*x4*x6+11*x4^2*x6-59*x1*x5*x6+61*x2*x5*x6+21*x3*x5*x6+37*x4*x5*x6+4*x5^2*x6+58*x1*x6^2-36*x2*x6^2+75*x3*x6^2-44*x4*x6^2-60*x5*x6^2+84*x6^3+52*x1^2+28*x1*x2+74*x2^2+17*x1*x3-56*x2*x3+87*x3^2-16*x1*x4+88*x2*x4+94*x3*x4-91*x4^2-47*x1*x5+46*x2*x5+18*x3*x5+76*x4*x5-84*x5^2-56*x1*x6-76*x2*x6+10*x3*x6+67*x4*x6-83*x5*x6+77*x6^2-45*x1-22*x2+18*x3-72*x4-46*x6+63");
    polCp_D3_N6_P6.emplace_back("27*x1^3-43*x1^2*x2-84*x1*x2^2-11*x2^3+72*x1^2*x3+78*x1*x2*x3+2*x2^2*x3+58*x1*x3^2+4*x2*x3^2-48*x3^3-64*x1^2*x4-12*x1*x2*x4-4*x2^2*x4-62*x1*x3*x4-89*x2*x3*x4+53*x3^2*x4+94*x1*x4^2-38*x2*x4^2+33*x3*x4^2-73*x4^3+8*x1^2*x5-2*x1*x2*x5+5*x2^2*x5-77*x1*x3*x5-54*x2*x3*x5-68*x3^2*x5+67*x1*x4*x5-50*x2*x4*x5-87*x3*x4*x5-88*x4^2*x5+45*x1*x5^2+80*x2*x5^2-8*x3*x5^2-86*x4*x5^2+94*x5^3+43*x1^2*x6-45*x1*x2*x6-40*x2^2*x6-45*x1*x3*x6-81*x2*x3*x6+27*x3^2*x6+95*x1*x4*x6-75*x2*x4*x6+24*x3*x4*x6-52*x4^2*x6-99*x1*x5*x6+14*x2*x5*x6+37*x3*x5*x6-80*x4*x5*x6+5*x5^2*x6+42*x1*x6^2+8*x2*x6^2-82*x3*x6^2-72*x4*x6^2-45*x5*x6^2-3*x6^3+58*x1^2-85*x1*x2+5*x2^2-34*x1*x3+22*x2*x3-5*x3^2-38*x1*x4+11*x2*x4-18*x3*x4+3*x4^2-80*x1*x5-56*x2*x5-96*x3*x5-96*x4*x5+65*x5^2+86*x1*x6+80*x2*x6+50*x3*x6-16*x4*x6-59*x5*x6-21*x6^2+7*x1-91*x2+88*x3-11*x4+83*x5-80*x6+98");

    // Create cp_d3_n6_p6 ideal;
    Ideal<eltType> cp_d3_n6_p6(polCp_D3_N6_P6, 6, 20000000);
    
    // Compute a reduced groebner basis;
    nbGen=cp_d3_n6_p6.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        cp_d3_n6_p6.printReducedGroebnerBasis("benchmark-cp_d3_n6_p6.res", modulo);
    }
    return nbGen;
}

int main (int argc, char **argv)
{
    // Time
    chrono::steady_clock::time_point start;
    typedef chrono::duration<int,milli> millisecs_t;
    
    // Number of thread(s)
    cout << NB_THREAD << " thread(s) used " << endl << endl;
    
    // Magma output
    bool magma = false;
    
    // Number of generator
    int nbGen;
    
    cout << "Benchmark for ideal with integer long type coefficient." << endl;
    
    start=chrono::steady_clock::now();
    nbGen=cp_d3_n6_p6F4(magma);

    cout << "cp(3,6,6): " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;
   
    return 0;
}


