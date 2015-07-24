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
 *  \file check-cyclic6-32bits.cpp
 *  \example check-cyclic6-32bits.cpp
 *  \brief Library check tests.
 *  \ingroup examples
 *  \author Vanessa VITSE, Antoine JOUX, Titouan COLADON
 */

#include <iostream>
#include <string>
#include <vector>
#include <libopenf4.h>

using namespace std;

int main (int argc, char **argv)
{
    cout << "#########################################################" << endl;
    cout << "#                 CHECK CYCLIC6 32 BITS                 #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Create polynomial array
    vector<string> polynomialArray;
    
    // Create variable name array
    vector<string> variableName;
    for(int i = 0; i < 6; i++)
    {
        variableName.push_back('x'+to_string(i));
    }
    
    // Fill the polynomial array
    polynomialArray.emplace_back("x0+x1+x2+x3+x4+x5");
    polynomialArray.emplace_back("x0*x1+x1*x2+x2*x3+x3*x4+x0*x5+x4*x5");
    polynomialArray.emplace_back("x0*x1*x2+x1*x2*x3+x2*x3*x4+x0*x1*x5+x0*x4*x5+x3*x4*x5");
    polynomialArray.emplace_back("x0*x1*x2*x3+x1*x2*x3*x4+x0*x1*x2*x5+x0*x1*x4*x5+x0*x3*x4*x5+x2*x3*x4*x5");
    polynomialArray.emplace_back("x0*x1*x2*x3*x4+x0*x1*x2*x3*x5+x0*x1*x2*x4*x5+x0*x1*x3*x4*x5+x0*x2*x3*x4*x5+x1*x2*x3*x4*x5");
    polynomialArray.emplace_back("x0*x1*x2*x3*x4*x5-1");
    
    // Compute a reduce groebner basis
    vector<string> basisInt = groebnerBasisF4(65521, 6, variableName, polynomialArray, 1, 0);
    
    // Fill reference vectors
    vector<string> groebnerBasisCyclic6Int;
    
    groebnerBasisCyclic6Int.push_back("(1*x0^1) + (1*x1^1) + (1*x2^1) + (1*x3^1) + (1*x4^1) + (1*x5^1)");
    groebnerBasisCyclic6Int.push_back("(1*x1^2) + (1*x1^1*x3^1) + (-1*x2^1*x3^1) + (1*x1^1*x4^1) + (-1*x3^1*x4^1) + (2*x1^1*x5^1) + (1*x2^1*x5^1) + (1*x3^1*x5^1) + (1*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x2^2) + (1*x2^1*x3^2) + (32756*x1^1*x2^1*x4^1) + (32759*x2^2*x4^1) + (3*x1^1*x3^1*x4^1) + (32759*x2^1*x3^1*x4^1) + (4*x3^2*x4^1) + (32759*x2^1*x4^2) + (1*x3^1*x4^2) + (3*x1^1*x2^1*x5^1) + (1*x2^2*x5^1) + (32759*x1^1*x3^1*x5^1) + (-4*x2^1*x3^1*x5^1) + (-3*x3^2*x5^1) + (-32757*x2^1*x4^1*x5^1) + (3*x3^1*x4^1*x5^1) + (1*x4^2*x5^1) + (-1*x1^1*x5^2) + (-4*x3^1*x5^2) + (2*x4^1*x5^2) + (-1*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x2^3) + (-28083*x2^1*x3^2) + (-28081*x3^3) + (28083*x1^1*x2^1*x4^1) + (18722*x2^2*x4^1) + (-9361*x1^1*x3^1*x4^1) + (28083*x2^1*x3^1*x4^1) + (28077*x3^2*x4^1) + (28077*x1^1*x4^2) + (18722*x2^1*x4^2) + (9361*x3^1*x4^2) + (-28081*x4^3) + (-18722*x1^1*x2^1*x5^1) + (-9358*x2^2*x5^1) + (18722*x1^1*x3^1*x5^1) + (9361*x3^2*x5^1) + (-9361*x1^1*x4^1*x5^1) + (-9361*x2^1*x4^1*x5^1) + (18716*x3^1*x4^1*x5^1) + (28083*x1^1*x5^2) + (18719*x2^1*x5^2) + (-28077*x3^1*x5^2) + (-18722*x4^1*x5^2) + (1*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x2^1*x3^1) + (1*x2^1*x3^2) + (32758*x1^1*x2^1*x4^1) + (32759*x2^2*x4^1) + (1*x1^1*x3^1*x4^1) + (32760*x2^1*x3^1*x4^1) + (2*x3^2*x4^1) + (32759*x2^1*x4^2) + (1*x3^1*x4^2) + (1*x1^1*x2^1*x5^1) + (1*x2^2*x5^1) + (-32760*x1^1*x3^1*x5^1) + (-1*x2^1*x3^1*x5^1) + (-1*x3^2*x5^1) + (-1*x1^1*x4^1*x5^1) + (32760*x2^1*x4^1*x5^1) + (1*x3^1*x4^1*x5^1) + (1*x2^1*x5^2) + (-1*x3^1*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x2^2*x3^1) + (2*x1^1*x2^1*x4^1) + (-2*x1^1*x3^1*x4^1) + (1*x2^1*x3^1*x4^1) + (-2*x3^2*x4^1) + (-2*x1^1*x2^1*x5^1) + (-1*x2^2*x5^1) + (2*x1^1*x3^1*x5^1) + (3*x2^1*x3^1*x5^1) + (2*x3^2*x5^1) + (-3*x2^1*x4^1*x5^1) + (-1*x3^1*x4^1*x5^1) + (-1*x2^1*x5^2) + (2*x3^1*x5^2) + (-1*x4^1*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x3^2) + (28081*x2^1*x3^2) + (28081*x3^3) + (4680*x1^1*x2^1*x4^1) + (14040*x2^2*x4^1) + (9360*x1^1*x3^1*x4^1) + (4680*x2^1*x3^1*x4^1) + (-28080*x3^2*x4^1) + (-28082*x1^1*x4^2) + (14040*x2^1*x4^2) + (-9360*x3^1*x4^2) + (28080*x4^3) + (18720*x1^1*x2^1*x5^1) + (9360*x2^2*x5^1) + (14041*x1^1*x3^1*x5^1) + (-9359*x3^2*x5^1) + (9360*x1^1*x4^1*x5^1) + (-23401*x2^1*x4^1*x5^1) + (-18721*x3^1*x4^1*x5^1) + (-28080*x1^1*x5^2) + (-18721*x2^1*x5^2) + (28081*x3^1*x5^2) + (18720*x4^1*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x2^1*x3^3) + (-30475*x4^4) + (7123*x2^1*x3^2*x5^1) + (7625*x3^3*x5^1) + (-27186*x1^1*x2^1*x4^1*x5^1) + (-25911*x2^2*x4^1*x5^1) + (20823*x1^1*x3^1*x4^1*x5^1) + (-19836*x2^1*x3^1*x4^1*x5^1) + (-19294*x3^2*x4^1*x5^1) + (-7077*x1^1*x4^2*x5^1) + (20569*x2^1*x4^2*x5^1) + (-7628*x3^1*x4^2*x5^1) + (-6084*x4^3*x5^1) + (-2027*x1^1*x2^1*x5^2) + (-3302*x2^2*x5^2) + (19033*x1^1*x3^1*x5^2) + (5341*x2^1*x3^1*x5^2) + (-1511*x3^2*x5^2) + (-22342*x1^1*x4^1*x5^2) + (-526*x2^1*x4^1*x5^2) + (22114*x3^1*x4^1*x5^2) + (-18283*x4^2*x5^2) + (18776*x1^1*x5^3) + (30731*x2^1*x5^3) + (10399*x3^1*x5^3) + (-2028*x4^1*x5^3) + (30975*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x3^4) + (18285*x4^4) + (6951*x2^1*x3^2*x5^1) + (-32656*x3^3*x5^1) + (18195*x1^1*x2^1*x4^1*x5^1) + (-24700*x2^2*x4^1*x5^1) + (4352*x1^1*x3^1*x4^1*x5^1) + (-32069*x2^1*x3^1*x4^1*x5^1) + (26543*x3^2*x4^1*x5^1) + (-7018*x1^1*x4^2*x5^1) + (19489*x2^1*x4^2*x5^1) + (13951*x3^1*x4^2*x5^1) + (-11332*x4^3*x5^1) + (-24065*x1^1*x2^1*x5^2) + (18830*x2^2*x5^2) + (-25437*x1^1*x3^1*x5^2) + (9890*x2^1*x3^1*x5^2) + (-2839*x3^2*x5^2) + (-2516*x1^1*x4^1*x5^2) + (30267*x2^1*x4^1*x5^2) + (31634*x3^1*x4^1*x5^2) + (17527*x4^2*x5^2) + (-29032*x1^1*x5^3) + (13383*x2^1*x5^3) + (-21200*x3^1*x5^3) + (-24065*x4^1*x5^3) + (7631*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x2^1*x3^2*x4^1) + (18723*x2^1*x3^2*x5^1) + (18722*x3^3*x5^1) + (-18723*x1^1*x2^1*x4^1*x5^1) + (9358*x2^2*x4^1*x5^1) + (28081*x1^1*x3^1*x4^1*x5^1) + (-18725*x2^1*x3^1*x4^1*x5^1) + (-18717*x3^2*x4^1*x5^1) + (-18714*x1^1*x4^2*x5^1) + (9358*x2^1*x4^2*x5^1) + (-28081*x3^1*x4^2*x5^1) + (18722*x4^3*x5^1) + (-9360*x1^1*x2^1*x5^2) + (28080*x2^2*x5^2) + (9359*x1^1*x3^1*x5^2) + (1*x2^1*x3^1*x5^2) + (-28076*x3^2*x5^2) + (28080*x1^1*x4^1*x5^2) + (28074*x2^1*x4^1*x5^2) + (9364*x3^1*x4^1*x5^2) + (-1*x4^2*x5^2) + (-18723*x1^1*x5^3) + (9361*x2^1*x5^3) + (18719*x3^1*x5^3) + (-9361*x4^1*x5^3) + (-1*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x3^3*x4^1) + (24380*x4^4) + (-9440*x2^1*x3^2*x5^1) + (3261*x3^3*x5^1) + (-719*x1^1*x2^1*x4^1*x5^1) + (12303*x2^2*x4^1*x5^1) + (-22273*x1^1*x3^1*x4^1*x5^1) + (-13152*x2^1*x3^1*x4^1*x5^1) + (32288*x3^2*x4^1*x5^1) + (9413*x1^1*x4^2*x5^1) + (-5225*x2^1*x4^2*x5^1) + (24823*x3^1*x4^2*x5^1) + (-11979*x4^3*x5^1) + (10046*x1^1*x2^1*x5^2) + (-2976*x2^2*x5^2) + (-10549*x1^1*x3^1*x5^2) + (-30481*x2^1*x3^1*x5^2) + (-19381*x3^2*x5^2) + (-20502*x1^1*x4^1*x5^2) + (-31409*x2^1*x4^1*x5^2) + (-32661*x3^1*x4^1*x5^2) + (-31238*x4^2*x5^2) + (-30937*x1^1*x5^3) + (-13354*x2^1*x5^3) + (-25171*x3^1*x5^3) + (10047*x4^1*x5^3) + (-11678*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x2^1*x4^2) + (-18285*x4^4) + (21115*x2^1*x3^2*x5^1) + (-4789*x3^3*x5^1) + (8343*x1^1*x2^1*x4^1*x5^1) + (5985*x2^2*x4^1*x5^1) + (5008*x1^1*x3^1*x4^1*x5^1) + (25868*x2^1*x3^1*x4^1*x5^1) + (-10955*x3^2*x4^1*x5^1) + (727*x1^1*x4^2*x5^1) + (-5443*x2^1*x4^2*x5^1) + (-1450*x3^1*x4^2*x5^1) + (-26121*x4^3*x5^1) + (-907*x1^1*x2^1*x5^2) + (1452*x2^2*x5^2) + (17667*x1^1*x3^1*x5^2) + (1015*x2^1*x3^1*x5^2) + (-28372*x3^2*x5^2) + (-31818*x1^1*x4^1*x5^2) + (-20897*x2^1*x4^1*x5^2) + (26050*x3^1*x4^1*x5^2) + (-17523*x4^2*x5^2) + (979*x1^1*x5^3) + (653*x2^1*x5^3) + (-5298*x3^1*x5^3) + (-907*x4^1*x5^3) + (-7619*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x2^2*x4^2) + (-12190*x4^4) + (10957*x2^1*x3^2*x5^1) + (-6313*x3^3*x5^1) + (-27718*x1^1*x2^1*x4^1*x5^1) + (-8491*x2^2*x4^1*x5^1) + (31420*x1^1*x3^1*x4^1*x5^1) + (-16033*x2^1*x3^1*x4^1*x5^1) + (-11462*x3^2*x4^1*x5^1) + (10881*x1^1*x4^2*x5^1) + (-16110*x2^1*x4^2*x5^1) + (76*x3^1*x4^2*x5^1) + (1305*x4^3*x5^1) + (-2685*x1^1*x2^1*x5^2) + (-21913*x2^2*x5^2) + (13860*x1^1*x3^1*x5^2) + (-6605*x2^1*x3^1*x5^2) + (11245*x3^2*x5^2) + (14148*x1^1*x4^1*x5^2) + (-29531*x2^1*x4^1*x5^2) + (-24235*x3^1*x4^1*x5^2) + (31999*x4^2*x5^2) + (25615*x1^1*x5^3) + (31634*x2^1*x5^3) + (7909*x3^1*x5^3) + (-2685*x4^1*x5^3) + (-5079*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x3^1*x4^2) + (-12190*x4^4) + (-14003*x2^1*x3^2*x5^1) + (12408*x3^3*x5^1) + (-13678*x1^1*x2^1*x4^1*x5^1) + (870*x2^2*x4^1*x5^1) + (-16942*x1^1*x3^1*x4^1*x5^1) + (8924*x2^1*x3^1*x4^1*x5^1) + (-30183*x3^2*x4^1*x5^1) + (-18755*x1^1*x4^2*x5^1) + (26012*x2^1*x4^2*x5^1) + (-17087*x3^1*x4^2*x5^1) + (-12734*x4^3*x5^1) + (-12045*x1^1*x2^1*x5^2) + (-26593*x2^2*x5^2) + (-20463*x1^1*x3^1*x5^2) + (-17523*x2^1*x3^1*x5^2) + (5007*x3^2*x5^2) + (31309*x1^1*x4^1*x5^2) + (20388*x2^1*x4^1*x5^2) + (6967*x3^1*x4^1*x5^2) + (-762*x4^2*x5^2) + (-14948*x1^1*x5^3) + (30075*x2^1*x5^3) + (-6132*x3^1*x5^3) + (-12045*x4^1*x5^3) + (-26920*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x2^1*x3^1*x4^2) + (-12190*x4^4) + (-4644*x2^1*x3^2*x5^1) + (21768*x3^3*x5^1) + (-12118*x1^1*x2^1*x4^1*x5^1) + (-27210*x2^2*x4^1*x5^1) + (-2902*x1^1*x3^1*x4^1*x5^1) + (-22276*x2^1*x3^1*x4^1*x5^1) + (-17704*x3^2*x4^1*x5^1) + (-17195*x1^1*x4^2*x5^1) + (30693*x2^1*x4^2*x5^1) + (12554*x3^1*x4^2*x5^1) + (29387*x4^3*x5^1) + (-5805*x1^1*x2^1*x5^2) + (9287*x2^2*x5^2) + (-26702*x1^1*x3^1*x5^2) + (-28444*x2^1*x3^1*x5^2) + (-19953*x3^2*x5^2) + (23509*x1^1*x4^1*x5^2) + (1667*x2^1*x4^1*x5^2) + (727*x3^1*x4^1*x5^2) + (31999*x4^2*x5^2) + (-24308*x1^1*x5^3) + (-30766*x2^1*x5^3) + (-7692*x3^1*x5^3) + (-5805*x4^1*x5^3) + (-26920*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x4^5) + (15227*x4^4*x5^1) + (-9943*x2^1*x3^2*x5^2) + (851*x3^3*x5^2) + (3715*x1^1*x2^1*x4^1*x5^2) + (-915*x2^2*x4^1*x5^2) + (1967*x1^1*x3^1*x4^1*x5^2) + (27335*x2^1*x3^1*x4^1*x5^2) + (-6181*x3^2*x4^1*x5^2) + (-1513*x1^1*x4^2*x5^2) + (-7775*x2^1*x4^2*x5^2) + (-30540*x3^1*x4^2*x5^2) + (7580*x4^3*x5^2) + (25805*x1^1*x2^1*x5^3) + (30435*x2^2*x5^3) + (-28700*x1^1*x3^1*x5^3) + (24634*x2^1*x3^1*x5^3) + (21240*x3^2*x5^3) + (-31054*x1^1*x4^1*x5^3) + (-20104*x2^1*x4^1*x5^3) + (-30903*x3^1*x4^1*x5^3) + (-23727*x4^2*x5^3) + (29780*x1^1*x5^4) + (8892*x2^1*x5^4) + (-27478*x3^1*x5^4) + (25804*x4^1*x5^4) + (-4433*x5^5)");
    groebnerBasisCyclic6Int.push_back("(1*x2^1*x3^2*x5^3) + (-14560*x3^3*x5^3) + (-14562*x1^1*x2^1*x4^1*x5^3) + (-32760*x2^2*x4^1*x5^3) + (-3644*x1^1*x3^1*x4^1*x5^3) + (3634*x2^1*x3^1*x4^1*x5^3) + (10918*x3^2*x4^1*x5^3) + (-21839*x1^1*x4^2*x5^3) + (-32760*x2^1*x4^2*x5^3) + (-29127*x3^1*x4^2*x5^3) + (-14560*x4^3*x5^3) + (-25477*x1^1*x2^1*x5^4) + (-7279*x2^2*x5^4) + (-10922*x1^1*x3^1*x5^4) + (-32757*x2^1*x3^1*x5^4) + (21840*x3^2*x5^4) + (10924*x1^1*x4^1*x5^4) + (-29117*x2^1*x4^1*x5^4) + (7281*x3^1*x4^1*x5^4) + (-1*x1^1*x5^5) + (-29118*x2^1*x5^5) + (21839*x3^1*x5^5) + (-25477*x4^1*x5^5) + (25480*x5^6) + (25480*1)");
    groebnerBasisCyclic6Int.push_back("(1*x4^4*x5^3) + (-20838*x1^1*x3^1*x4^1*x5^4) + (-5615*x2^1*x3^1*x4^1*x5^4) + (1072*x3^2*x4^1*x5^4) + (5247*x1^1*x4^2*x5^4) + (14199*x2^1*x4^2*x5^4) + (3074*x3^1*x4^2*x5^4) + (1298*x4^3*x5^4) + (-11534*x1^1*x2^1*x5^5) + (-11534*x2^2*x5^5) + (8010*x1^1*x3^1*x5^5) + (15223*x2^1*x3^1*x5^5) + (20070*x3^2*x5^5) + (-12767*x1^1*x4^1*x5^5) + (11120*x2^1*x4^1*x5^5) + (-14434*x3^1*x4^1*x5^5) + (-13286*x4^2*x5^5) + (10525*x1^1*x5^6) + (-20964*x2^1*x5^6) + (-32498*x3^1*x5^6) + (17925*x4^1*x5^6) + (26513*x5^7) + (21357*x1^1) + (16956*x2^1) + (604*x3^1) + (29047*x4^1) + (6750*x5^1)");
    groebnerBasisCyclic6Int.push_back("(1*x3^3*x5^4) + (17036*x1^1*x3^1*x4^1*x5^4) + (-23953*x2^1*x3^1*x4^1*x5^4) + (-21959*x3^2*x4^1*x5^4) + (4670*x1^1*x4^2*x5^4) + (-13052*x2^1*x4^2*x5^4) + (26689*x3^1*x4^2*x5^4) + (-32159*x4^3*x5^4) + (-23827*x1^1*x2^1*x5^5) + (-23827*x2^2*x5^5) + (-29723*x1^1*x3^1*x5^5) + (19655*x2^1*x3^1*x5^5) + (-24030*x3^2*x5^5) + (20274*x1^1*x4^1*x5^5) + (-19471*x2^1*x4^1*x5^5) + (-32370*x3^1*x4^1*x5^5) + (-32386*x4^2*x5^5) + (1515*x1^1*x5^6) + (9018*x2^1*x5^6) + (-13953*x3^1*x5^6) + (-25537*x4^1*x5^6) + (3846*x5^7) + (10055*x1^1) + (21843*x2^1) + (-16310*x3^1) + (23598*x4^1) + (-22206*x5^1)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x2^1*x4^1*x5^4) + (-10430*x1^1*x3^1*x4^1*x5^4) + (847*x2^1*x3^1*x4^1*x5^4) + (-19558*x3^2*x4^1*x5^4) + (29016*x1^1*x4^2*x5^4) + (24218*x2^1*x4^2*x5^4) + (29208*x3^1*x4^2*x5^4) + (-14864*x4^3*x5^4) + (-15997*x1^1*x2^1*x5^5) + (-15996*x2^2*x5^5) + (24559*x1^1*x3^1*x5^5) + (-29026*x2^1*x3^1*x5^5) + (7375*x3^2*x5^5) + (-14363*x1^1*x4^1*x5^5) + (-30277*x2^1*x4^1*x5^5) + (30363*x3^1*x4^1*x5^5) + (-22944*x4^2*x5^5) + (-15337*x1^1*x5^6) + (3169*x2^1*x5^6) + (-24150*x3^1*x5^6) + (-15842*x4^1*x5^6) + (-2108*x5^7) + (2551*x1^1) + (-16919*x2^1) + (-4582*x3^1) + (31052*x4^1) + (4513*x5^1)");
    groebnerBasisCyclic6Int.push_back("(1*x2^2*x4^1*x5^4) + (-19344*x1^1*x3^1*x4^1*x5^4) + (4877*x2^1*x3^1*x4^1*x5^4) + (-2809*x3^2*x4^1*x5^4) + (-19743*x1^1*x4^2*x5^4) + (25862*x2^1*x4^2*x5^4) + (-22625*x3^1*x4^2*x5^4) + (15208*x4^3*x5^4) + (-24379*x1^1*x2^1*x5^5) + (-24380*x2^2*x5^5) + (30991*x1^1*x3^1*x5^5) + (-21238*x2^1*x3^1*x5^5) + (10827*x3^2*x5^5) + (23600*x1^1*x4^1*x5^5) + (28454*x2^1*x4^1*x5^5) + (-6474*x3^1*x4^1*x5^5) + (-843*x4^2*x5^5) + (-25264*x1^1*x5^6) + (21564*x2^1*x5^6) + (-21713*x3^1*x5^6) + (-2430*x4^1*x5^6) + (9414*x5^7) + (-31382*x1^1) + (-1614*x2^1) + (22077*x3^1) + (-7932*x4^1) + (-26226*x5^1)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x3^1*x4^1*x5^5) + (-15787*x4^3*x5^5) + (-5468*x2^1*x3^1*x5^6) + (-29397*x3^2*x5^6) + (-4705*x1^1*x4^1*x5^6) + (28675*x2^1*x4^1*x5^6) + (-2742*x3^1*x4^1*x5^6) + (18762*x4^2*x5^6) + (12178*x1^1*x5^7) + (997*x2^1*x5^7) + (26429*x3^1*x5^7) + (-22504*x4^1*x5^7) + (-6979*x5^8) + (-7719*x1^1*x2^1) + (-7719*x2^2) + (16188*x1^1*x3^1) + (-15348*x2^1*x3^1) + (16300*x3^2) + (-4615*x1^1*x4^1) + (14526*x2^1*x4^1) + (-18439*x3^1*x4^1) + (-8685*x4^2) + (-11328*x1^1*x5^1) + (18982*x2^1*x5^1) + (12512*x3^1*x5^1) + (24832*x4^1*x5^1) + (-28947*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x2^1*x3^1*x4^1*x5^5) + (-6614*x4^3*x5^5) + (16502*x2^1*x3^1*x5^6) + (-4649*x3^2*x5^6) + (27006*x1^1*x4^1*x5^6) + (13614*x2^1*x4^1*x5^6) + (1180*x3^1*x4^1*x5^6) + (2522*x4^2*x5^6) + (30576*x1^1*x5^7) + (28129*x2^1*x5^7) + (-16312*x3^1*x5^7) + (20265*x4^1*x5^7) + (-4272*x5^8) + (12802*x1^1*x2^1) + (12802*x2^2) + (515*x1^1*x3^1) + (30340*x2^1*x3^1) + (-26831*x3^2) + (23379*x1^1*x4^1) + (23386*x2^1*x4^1) + (-9062*x3^1*x4^1) + (1235*x4^2) + (-28757*x1^1*x5^1) + (-8045*x2^1*x5^1) + (-18042*x3^1*x5^1) + (14469*x4^1*x5^1) + (-5097*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x3^2*x4^1*x5^5) + (9335*x4^3*x5^5) + (21257*x2^1*x3^1*x5^6) + (3093*x3^2*x5^6) + (6666*x1^1*x4^1*x5^6) + (-531*x2^1*x4^1*x5^6) + (23948*x3^1*x4^1*x5^6) + (25519*x4^2*x5^6) + (31246*x1^1*x5^7) + (-17392*x2^1*x5^7) + (-26155*x3^1*x5^7) + (-17814*x4^1*x5^7) + (17284*x5^8) + (-11394*x1^1*x2^1) + (-11394*x2^2) + (-14880*x1^1*x3^1) + (25636*x2^1*x3^1) + (-10328*x3^2) + (24695*x1^1*x4^1) + (-4180*x2^1*x4^1) + (-14209*x3^1*x4^1) + (19455*x4^2) + (29188*x1^1*x5^1) + (6549*x2^1*x5^1) + (20988*x3^1*x5^1) + (21301*x4^1*x5^1) + (-26842*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x4^2*x5^5) + (30096*x4^3*x5^5) + (10440*x2^1*x3^1*x5^6) + (-31920*x3^2*x5^6) + (12791*x1^1*x4^1*x5^6) + (-11737*x2^1*x4^1*x5^6) + (26065*x3^1*x4^1*x5^6) + (-25249*x4^2*x5^6) + (-21289*x1^1*x5^7) + (26278*x2^1*x5^7) + (11855*x3^1*x5^7) + (-30588*x4^1*x5^7) + (13152*x5^8) + (-16178*x1^1*x2^1) + (-16178*x2^2) + (-16551*x1^1*x3^1) + (13560*x2^1*x3^1) + (6577*x3^2) + (-16924*x1^1*x4^1) + (-18281*x2^1*x4^1) + (26445*x3^1*x4^1) + (-3174*x4^2) + (-7371*x1^1*x5^1) + (-3273*x2^1*x5^1) + (23377*x3^1*x5^1) + (-28731*x4^1*x5^1) + (-18714*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x2^1*x4^2*x5^5) + (32568*x4^3*x5^5) + (-7574*x2^1*x3^1*x5^6) + (8120*x3^2*x5^6) + (21973*x1^1*x4^1*x5^6) + (-1024*x2^1*x4^1*x5^6) + (-13680*x3^1*x4^1*x5^6) + (-5516*x4^2*x5^6) + (11103*x1^1*x5^7) + (26748*x2^1*x5^7) + (1694*x3^1*x5^7) + (5199*x4^1*x5^7) + (8446*x5^8) + (16778*x1^1*x2^1) + (16778*x2^2) + (20756*x1^1*x3^1) + (31610*x2^1*x3^1) + (-862*x3^2) + (-17*x1^1*x4^1) + (-18422*x2^1*x4^1) + (-3102*x3^1*x4^1) + (13658*x4^2) + (-5072*x1^1*x5^1) + (18995*x2^1*x5^1) + (-672*x3^1*x5^1) + (11772*x4^1*x5^1) + (6305*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x3^1*x4^2*x5^5) + (-27503*x4^3*x5^5) + (22150*x2^1*x3^1*x5^6) + (-29862*x3^2*x5^6) + (29079*x1^1*x4^1*x5^6) + (-26116*x2^1*x4^1*x5^6) + (5091*x3^1*x4^1*x5^6) + (21203*x4^2*x5^6) + (30823*x1^1*x5^7) + (3388*x2^1*x5^7) + (18041*x3^1*x5^7) + (-13747*x4^1*x5^7) + (-11702*x5^8) + (5064*x1^1*x2^1) + (5064*x2^2) + (-28726*x1^1*x3^1) + (17824*x2^1*x3^1) + (-749*x3^2) + (20898*x1^1*x4^1) + (24675*x2^1*x4^1) + (-30593*x3^1*x4^1) + (32440*x4^2) + (8383*x1^1*x5^1) + (-26729*x2^1*x5^1) + (4810*x3^1*x5^1) + (28139*x4^1*x5^1) + (-15825*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x2^1*x5^6) + (-10006*x2^1*x3^1*x5^6) + (-10007*x3^2*x5^6) + (3801*x1^1*x4^1*x5^6) + (-6206*x2^1*x4^1*x5^6) + (-16213*x3^1*x4^1*x5^6) + (-6206*x4^2*x5^6) + (18414*x1^1*x5^7) + (-651*x2^1*x5^7) + (8407*x3^1*x5^7) + (-6206*x4^1*x5^7) + (-26260*x5^8) + (-27461*x1^1*x2^1) + (-27460*x2^2) + (12389*x1^1*x3^1) + (28053*x2^1*x3^1) + (-10007*x3^2) + (1904*x1^1*x4^1) + (-8103*x2^1*x4^1) + (-18110*x3^1*x4^1) + (-8103*x4^2) + (-9048*x1^1*x5^1) + (-19946*x2^1*x5^1) + (-19055*x3^1*x5^1) + (-8103*x4^1*x5^1) + (-32381*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x2^2*x5^6) + (6149*x2^1*x3^1*x5^6) + (6149*x3^2*x5^6) + (-28734*x1^1*x4^1*x5^6) + (-22585*x2^1*x4^1*x5^6) + (-16436*x3^1*x4^1*x5^6) + (-22585*x4^2*x5^6) + (-23755*x1^1*x5^7) + (-31874*x2^1*x5^7) + (-17606*x3^1*x5^7) + (-22585*x4^1*x5^7) + (20478*x5^8) + (27575*x1^1*x2^1) + (27574*x2^2) + (-21685*x1^1*x3^1) + (-31797*x2^1*x3^1) + (6149*x3^2) + (-22742*x1^1*x4^1) + (-16593*x2^1*x4^1) + (-10444*x3^1*x4^1) + (-16593*x4^2) + (3820*x1^1*x5^1) + (-17138*x2^1*x5^1) + (9969*x3^1*x5^1) + (-16593*x4^1*x5^1) + (-30203*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x3^1*x5^6) + (-24870*x2^1*x3^1*x5^6) + (-24870*x3^2*x5^6) + (20319*x1^1*x4^1*x5^6) + (-4551*x2^1*x4^1*x5^6) + (-29421*x3^1*x4^1*x5^6) + (-4551*x4^2*x5^6) + (-23950*x1^1*x5^7) + (13036*x2^1*x5^7) + (16701*x3^1*x5^7) + (-4551*x4^1*x5^7) + (-2426*x5^8) + (-27302*x1^1*x2^1) + (-27302*x2^2) + (21840*x1^1*x3^1) + (13349*x2^1*x3^1) + (-24870*x3^2) + (-5177*x1^1*x4^1) + (-30047*x2^1*x4^1) + (10604*x3^1*x4^1) + (-30047*x4^2) + (14269*x1^1*x5^1) + (16698*x2^1*x5^1) + (-10601*x3^1*x5^1) + (-30047*x4^1*x5^1) + (-18797*x5^2)");
    groebnerBasisCyclic6Int.push_back("(1*x4^3*x5^6) + (9305*x4^1*x5^8) + (8470*x2^1*x3^2) + (-20357*x3^3) + (-4961*x1^1*x2^1*x4^1) + (-15292*x2^2*x4^1) + (21163*x1^1*x3^1*x4^1) + (10446*x2^1*x3^1*x4^1) + (13555*x3^2*x4^1) + (-16430*x1^1*x4^2) + (-26671*x2^1*x4^2) + (-5567*x3^1*x4^2) + (21262*x4^3) + (32682*x1^1*x2^1*x5^1) + (-22508*x2^2*x5^1) + (-20011*x1^1*x3^1*x5^1) + (15317*x2^1*x3^1*x5^1) + (-2593*x3^2*x5^1) + (15615*x1^1*x4^1*x5^1) + (-19458*x2^1*x4^1*x5^1) + (-27026*x3^1*x4^1*x5^1) + (22003*x4^2*x5^1) + (-28058*x1^1*x5^2) + (29538*x2^1*x5^2) + (10336*x3^1*x5^2) + (28086*x4^1*x5^2) + (-28847*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x2^1*x3^1*x5^7) + (15224*x4^1*x5^8) + (18759*x2^1*x3^2) + (-22072*x3^3) + (8764*x1^1*x2^1*x4^1) + (6701*x2^2*x4^1) + (-19015*x1^1*x3^1*x4^1) + (-5141*x2^1*x3^1*x4^1) + (2791*x3^2*x4^1) + (21511*x1^1*x4^2) + (-13255*x2^1*x4^2) + (-22152*x3^1*x4^2) + (-8489*x4^3) + (27991*x1^1*x2^1*x5^1) + (30054*x2^2*x5^1) + (-28612*x1^1*x3^1*x5^1) + (-2152*x2^1*x3^1*x5^1) + (-28774*x3^2*x5^1) + (-104*x1^1*x4^1*x5^1) + (1*x2^1*x4^1*x5^1) + (-6081*x3^1*x4^1*x5^1) + (21807*x4^2*x5^1) + (-10535*x1^1*x5^2) + (-18990*x2^1*x5^2) + (10698*x3^1*x5^2) + (-3620*x4^1*x5^2) + (24690*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x3^2*x5^7) + (-27227*x4^1*x5^8) + (-2767*x2^1*x3^2) + (1016*x3^3) + (-14775*x1^1*x2^1*x4^1) + (18905*x2^2*x4^1) + (-22600*x1^1*x3^1*x4^1) + (-28131*x2^1*x3^1*x4^1) + (17554*x3^2*x4^1) + (24619*x1^1*x4^2) + (-20794*x2^1*x4^2) + (549*x3^1*x4^2) + (-14281*x4^3) + (10271*x1^1*x2^1*x5^1) + (-23409*x2^2*x5^1) + (-6384*x1^1*x3^1*x5^1) + (12893*x2^1*x3^1*x5^1) + (-28008*x3^2*x5^1) + (306*x1^1*x4^1*x5^1) + (-2566*x2^1*x4^1*x5^1) + (26866*x3^1*x4^1*x5^1) + (-9492*x4^2*x5^1) + (8563*x1^1*x5^2) + (27853*x2^1*x5^2) + (30810*x3^1*x5^2) + (16438*x4^1*x5^2) + (3790*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x4^1*x5^7) + (17757*x4^1*x5^8) + (-11984*x2^1*x3^2) + (-24216*x3^3) + (-23379*x1^1*x2^1*x4^1) + (-4821*x2^2*x4^1) + (23373*x1^1*x3^1*x4^1) + (-26616*x2^1*x3^1*x4^1) + (3344*x3^2*x4^1) + (21132*x1^1*x4^2) + (-30179*x2^1*x4^2) + (-24343*x3^1*x4^2) + (22799*x4^3) + (4809*x1^1*x2^1*x5^1) + (-13749*x2^2*x5^1) + (-675*x1^1*x3^1*x5^1) + (-18536*x2^1*x3^1*x5^1) + (16484*x3^2*x5^1) + (-20889*x1^1*x4^1*x5^1) + (-23117*x2^1*x4^1*x5^1) + (15511*x3^1*x4^1*x5^1) + (30639*x4^2*x5^1) + (-4372*x1^1*x5^2) + (-10799*x2^1*x5^2) + (-26362*x3^1*x5^2) + (-10606*x4^1*x5^2) + (-12248*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x2^1*x4^1*x5^7) + (5754*x4^1*x5^8) + (-25189*x2^1*x3^2) + (-8005*x3^3) + (-10179*x1^1*x2^1*x4^1) + (8881*x2^2*x4^1) + (-28889*x1^1*x3^1*x4^1) + (23479*x2^1*x3^1*x4^1) + (27529*x3^2*x4^1) + (-14892*x1^1*x4^2) + (-10612*x2^1*x4^2) + (-16194*x3^1*x4^2) + (-28223*x4^3) + (12784*x1^1*x2^1*x5^1) + (-6276*x2^2*x5^1) + (27380*x1^1*x3^1*x5^1) + (3677*x2^1*x3^1*x5^1) + (13339*x3^2*x5^1) + (24517*x1^1*x4^1*x5^1) + (31378*x2^1*x4^1*x5^1) + (-2020*x3^1*x4^1*x5^1) + (9393*x4^2*x5^1) + (-10721*x1^1*x5^2) + (-14919*x2^1*x5^2) + (-1113*x3^1*x5^2) + (-28064*x4^1*x5^2) + (17184*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x3^1*x4^1*x5^7) + (-13627*x4^1*x5^8) + (31612*x2^1*x3^2) + (28812*x3^3) + (17601*x1^1*x2^1*x4^1) + (28588*x2^2*x4^1) + (-6403*x1^1*x3^1*x4^1) + (-24071*x2^1*x3^1*x4^1) + (920*x3^2*x4^1) + (11619*x1^1*x4^2) + (2369*x2^1*x4^2) + (-420*x3^1*x4^2) + (-16404*x4^3) + (-17967*x1^1*x2^1*x5^1) + (-28954*x2^2*x5^1) + (12541*x1^1*x3^1*x5^1) + (11919*x2^1*x3^1*x5^1) + (-25700*x3^2*x5^1) + (22221*x1^1*x4^1*x5^1) + (23121*x2^1*x4^1*x5^1) + (-28136*x3^1*x4^1*x5^1) + (12325*x4^2*x5^1) + (25909*x1^1*x5^2) + (19473*x2^1*x5^2) + (29734*x3^1*x5^2) + (16743*x4^1*x5^2) + (-2784*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x4^2*x5^7) + (-9883*x4^1*x5^8) + (-25365*x2^1*x3^2) + (-5917*x3^3) + (13846*x1^1*x2^1*x4^1) + (27598*x2^2*x4^1) + (9158*x1^1*x3^1*x4^1) + (-18155*x2^1*x3^1*x4^1) + (21891*x3^2*x4^1) + (-8294*x1^1*x4^2) + (-32374*x2^1*x4^2) + (869*x3^1*x4^2) + (12502*x4^3) + (-5739*x1^1*x2^1*x5^1) + (-19491*x2^2*x5^1) + (-22196*x1^1*x3^1*x5^1) + (24573*x2^1*x3^1*x5^1) + (-1175*x3^2*x5^1) + (-7116*x1^1*x4^1*x5^1) + (-1999*x2^1*x4^1*x5^1) + (-309*x3^1*x4^1*x5^1) + (13443*x4^2*x5^1) + (20341*x1^1*x5^2) + (12120*x2^1*x5^2) + (31846*x3^1*x5^2) + (15862*x4^1*x5^2) + (19484*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x5^8) + (22134*x4^1*x5^8) + (-15653*x2^1*x3^2) + (3371*x3^3) + (-17812*x1^1*x2^1*x4^1) + (19029*x2^2*x4^1) + (2506*x1^1*x3^1*x4^1) + (21233*x2^1*x3^1*x4^1) + (-5279*x3^2*x4^1) + (-3947*x1^1*x4^2) + (24034*x2^1*x4^2) + (28745*x3^1*x4^2) + (-25053*x4^3) + (5739*x1^1*x2^1*x5^1) + (-31102*x2^2*x5^1) + (9774*x1^1*x3^1*x5^1) + (-11532*x2^1*x3^1*x5^1) + (21513*x3^2*x5^1) + (1072*x1^1*x4^1*x5^1) + (26282*x2^1*x4^1*x5^1) + (-927*x3^1*x4^1*x5^1) + (24083*x4^2*x5^1) + (2667*x1^1*x5^2) + (-15062*x2^1*x5^2) + (-4122*x3^1*x5^2) + (30344*x4^1*x5^2) + (19004*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x2^1*x5^8) + (-14673*x2^1*x3^2) + (9027*x3^3) + (16548*x1^1*x2^1*x4^1) + (-23829*x2^2*x4^1) + (22400*x1^1*x3^1*x4^1) + (-29887*x2^1*x3^1*x4^1) + (29564*x3^2*x4^1) + (-11359*x1^1*x4^2) + (-23829*x2^1*x4^2) + (-6058*x3^1*x4^2) + (9027*x4^3) + (15915*x1^1*x2^1*x5^1) + (-9229*x2^2*x5^1) + (27759*x1^1*x3^1*x5^1) + (-9782*x2^1*x3^1*x5^1) + (-23767*x3^2*x5^1) + (13965*x1^1*x4^1*x5^1) + (-6004*x2^1*x4^1*x5^1) + (21692*x3^1*x4^1*x5^1) + (31121*x4^2*x5^1) + (-19707*x1^1*x5^2) + (-15000*x2^1*x5^2) + (22047*x3^1*x5^2) + (15895*x4^1*x5^2) + (23684*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x3^1*x5^8) + (-22134*x4^1*x5^8) + (20667*x2^1*x3^2) + (5327*x3^3) + (-16258*x1^1*x2^1*x4^1) + (-5939*x2^2*x4^1) + (-26290*x1^1*x3^1*x4^1) + (4781*x2^1*x3^1*x4^1) + (-13496*x3^2*x4^1) + (-9386*x1^1*x4^2) + (-10944*x2^1*x4^2) + (-15821*x3^1*x4^2) + (-31770*x4^3) + (-12895*x1^1*x2^1*x5^1) + (-23214*x2^2*x5^1) + (30132*x1^1*x3^1*x5^1) + (-18990*x2^1*x3^1*x5^1) + (-433*x3^2*x5^1) + (12820*x1^1*x4^1*x5^1) + (24071*x2^1*x4^1*x5^1) + (-3951*x3^1*x4^1*x5^1) + (5447*x4^2*x5^1) + (21877*x1^1*x5^2) + (23998*x2^1*x5^2) + (-15775*x3^1*x5^2) + (27994*x4^1*x5^2) + (-15340*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x5^9) + (-19818*x2^1*x3^2) + (-26831*x3^3) + (-20037*x1^1*x2^1*x4^1) + (-30106*x2^2*x4^1) + (314*x1^1*x3^1*x4^1) + (-9597*x2^1*x3^1*x4^1) + (-19483*x3^2*x4^1) + (-29944*x1^1*x4^2) + (-30106*x2^1*x4^2) + (20509*x3^1*x4^2) + (-26831*x4^3) + (18711*x1^1*x2^1*x5^1) + (28780*x2^2*x5^1) + (14088*x1^1*x3^1*x5^1) + (28424*x2^1*x3^1*x5^1) + (4846*x3^2*x5^1) + (23686*x1^1*x4^1*x5^1) + (-16154*x2^1*x4^1*x5^1) + (4182*x3^1*x4^1*x5^1) + (-26528*x4^2*x5^1) + (-6818*x1^1*x5^2) + (-22248*x2^1*x5^2) + (-1972*x3^1*x5^2) + (18819*x4^1*x5^2) + (-6929*x5^3)");
    groebnerBasisCyclic6Int.push_back("(1*x3^2*x4^2) + (-12190*x4^4) + (4715*x2^1*x3^2*x5^1) + (31127*x3^3*x5^1) + (364*x1^1*x2^1*x4^1*x5^1) + (10231*x2^2*x4^1*x5^1) + (-21624*x1^1*x3^1*x4^1*x5^1) + (-9793*x2^1*x3^1*x4^1*x5^1) + (16614*x3^2*x4^1*x5^1) + (-4718*x1^1*x4^2*x5^1) + (2614*x2^1*x4^2*x5^1) + (-12404*x3^1*x4^2*x5^1) + (-26774*x4^3*x5^1) + (-21406*x1^1*x2^1*x5^2) + (-31273*x2^2*x5^2) + (-11099*x1^1*x3^1*x5^2) + (15236*x2^1*x3^1*x5^2) + (-23076*x3^2*x5^2) + (26628*x1^1*x4^1*x5^2) + (-17052*x2^1*x4^1*x5^2) + (16323*x3^1*x4^1*x5^2) + (32001*x4^2*x5^2) + (31855*x1^1*x5^3) + (6675*x2^1*x5^3) + (-20171*x3^1*x5^3) + (-21406*x4^1*x5^3) + (-26919*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x1^1*x4^3) + (7619*x4^4) + (9142*x2^1*x3^2*x5^1) + (-18285*x3^3*x5^1) + (14984*x1^1*x2^1*x4^1*x5^1) + (22856*x2^2*x4^1*x5^1) + (22095*x1^1*x3^1*x4^1*x5^1) + (-22346*x2^1*x3^1*x4^1*x5^1) + (-27935*x3^2*x4^1*x5^1) + (-20064*x1^1*x4^2*x5^1) + (-29713*x2^1*x4^2*x5^1) + (7367*x3^1*x4^2*x5^1) + (-31237*x4^3*x5^1) + (507*x1^1*x2^1*x5^2) + (-7365*x2^2*x5^2) + (-2031*x1^1*x3^1*x5^2) + (-14985*x2^1*x3^1*x5^2) + (-5079*x3^2*x5^2) + (16506*x1^1*x4^1*x5^2) + (27427*x2^1*x4^1*x5^2) + (29966*x3^1*x4^1*x5^2) + (-28189*x4^2*x5^2) + (-31998*x1^1*x5^3) + (11427*x2^1*x5^3) + (-27173*x3^1*x5^3) + (507*x4^1*x5^3) + (30475*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x2^1*x4^3) + (1524*x4^4) + (3699*x2^1*x3^2*x5^1) + (11319*x3^3*x5^1) + (-31633*x1^1*x2^1*x4^1*x5^1) + (18613*x2^2*x4^1*x5^1) + (-5879*x1^1*x3^1*x4^1*x5^1) + (13317*x2^1*x3^1*x4^1*x5^1) + (31851*x3^2*x4^1*x5^1) + (-25543*x1^1*x4^2*x5^1) + (-5004*x2^1*x4^2*x5^1) + (-14438*x3^1*x4^2*x5^1) + (2177*x4^3*x5^1) + (-837*x1^1*x2^1*x5^2) + (14438*x2^2*x5^2) + (26742*x1^1*x3^1*x5^2) + (3555*x2^1*x3^1*x5^2) + (-10374*x3^2*x5^2) + (12660*x1^1*x4^1*x5^2) + (1738*x2^1*x4^1*x5^2) + (-25835*x3^1*x4^1*x5^2) + (20571*x4^2*x5^2) + (24490*x1^1*x5^3) + (16326*x2^1*x5^3) + (-23216*x3^1*x5^3) + (-837*x4^1*x5^3) + (6096*x5^4)");
    groebnerBasisCyclic6Int.push_back("(1*x3^1*x4^3) + (-10666*x4^4) + (-16544*x2^1*x3^2*x5^1) + (-4354*x3^3*x5^1) + (-28152*x1^1*x2^1*x4^1*x5^1) + (5442*x2^2*x4^1*x5^1) + (22423*x1^1*x3^1*x4^1*x5^1) + (17563*x2^1*x3^1*x4^1*x5^1) + (7910*x3^2*x4^1*x5^1) + (-5301*x1^1*x4^2*x5^1) + (6964*x2^1*x4^2*x5^1) + (10596*x3^1*x4^2*x5^1) + (-5878*x4^3*x5^1) + (23000*x1^1*x2^1*x5^2) + (-10594*x2^2*x5^2) + (-7762*x1^1*x3^1*x5^2) + (18791*x2^1*x3^1*x5^2) + (3990*x3^2*x5^2) + (12769*x1^1*x4^1*x5^2) + (-30910*x2^1*x4^1*x5^2) + (-147*x3^1*x4^1*x5^2) + (19808*x4^2*x5^2) + (-16977*x1^1*x5^3) + (-11320*x2^1*x5^3) + (19012*x3^1*x5^3) + (23000*x4^1*x5^3) + (22857*x5^4)");

    /* 32 bits */
    bool testCyclic6 = true;
    size_t i = 0;
    while(i < basisInt.size() && testCyclic6 == true )
    {
        testCyclic6 = testCyclic6 && (groebnerBasisCyclic6Int[i].compare(basisInt[i])==0);
        i++;
    }
    if(testCyclic6==true)
    {
        cout << "Test cyclic6 32 bits pass" << endl;
        return 0;
    }
    else
    {
        cout << "Test cyclic6 32 bits failed" << endl;
        return -1;
    }
}


