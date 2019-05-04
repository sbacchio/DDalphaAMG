/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */
 
#ifndef CLIFFORD_HEADER
  #define CLIFFORD_HEADER
  // assertion: gamma5 = (+/-) diag( 1, 1, -1, -1 )
  
  // choose basis:
  // BASIS0: Basis used within the OPENQCD/DD-HMC code
  // BASIS1: Basis used for BMW-c lattices
  // BASIS2: Basis used for QCDSF lattices
  // BASIS3: Basis used in the QOPQDP Code
  // BASIS4: Basis used in the tmLQCD Code
  #define BASIS4 // change here
  
  enum { T, Z, Y, X };
  
  #ifndef I
    #define I _Complex_I
  #endif
  
  #ifdef BASIS0
    // Basis used within the OPENQCD/DD-HMC code
    #define CLIFFORD_BASIS "BASIS0:OPENQCD/DD-HMC BASIS"
    /* gamma_T = 
     *  0  0 -1  0
     *  0  0  0 -1
     * -1  0  0  0
     *  0 -1  0  0
     */
    #define GAMMA_T_SPIN0_CO   2
    #define GAMMA_T_SPIN0_VAL -1
    #define GAMMA_T_SPIN1_CO   3
    #define GAMMA_T_SPIN1_VAL -1
    #define GAMMA_T_SPIN2_CO   0
    #define GAMMA_T_SPIN2_VAL -1
    #define GAMMA_T_SPIN3_CO   1
    #define GAMMA_T_SPIN3_VAL -1
    
    /* gamma_Z =
     *  0  0  0 -I
     *  0  0 -I  0
     *  0  I  0  0
     *  I  0  0  0
     */
    #define GAMMA_Z_SPIN0_CO   3
    #define GAMMA_Z_SPIN0_VAL -I
    #define GAMMA_Z_SPIN1_CO   2
    #define GAMMA_Z_SPIN1_VAL -I
    #define GAMMA_Z_SPIN2_CO   1
    #define GAMMA_Z_SPIN2_VAL  I
    #define GAMMA_Z_SPIN3_CO   0
    #define GAMMA_Z_SPIN3_VAL  I
    
    /* gamma_Y =
     *  0  0  0 -1
     *  0  0  1  0
     *  0  1  0  0
     * -1  0  0  0
     */
    #define GAMMA_Y_SPIN0_CO   3
    #define GAMMA_Y_SPIN0_VAL -1
    #define GAMMA_Y_SPIN1_CO   2
    #define GAMMA_Y_SPIN1_VAL  1
    #define GAMMA_Y_SPIN2_CO   1
    #define GAMMA_Y_SPIN2_VAL  1
    #define GAMMA_Y_SPIN3_CO   0
    #define GAMMA_Y_SPIN3_VAL -1   
    
    /* gamma_X = 
     *  0  0 -I  0
     *  0  0  0  I
     *  I  0  0  0
     *  0 -I  0  0
     */
    #define GAMMA_X_SPIN0_CO   2
    #define GAMMA_X_SPIN0_VAL -I
    #define GAMMA_X_SPIN1_CO   3
    #define GAMMA_X_SPIN1_VAL  I
    #define GAMMA_X_SPIN2_CO   0
    #define GAMMA_X_SPIN2_VAL  I
    #define GAMMA_X_SPIN3_CO   1
    #define GAMMA_X_SPIN3_VAL -I
  
  /* ------------------------------------------------- */
  #else
    #ifdef BASIS1
      // Basis used for BMW-c lattices
      #define CLIFFORD_BASIS "BASIS1:BMW-c BASIS"
      // Stored negatively due to different interpretations
      // of the Wilson Dirac operator
      /* gamma_T =
      *  0  0  1  0 
      *  0  0  0  1
      *  1  0  0  0
      *  0  1  0  0  
      */
      #define GAMMA_T_SPIN0_CO   2
      #define GAMMA_T_SPIN0_VAL -1
      #define GAMMA_T_SPIN1_CO   3
      #define GAMMA_T_SPIN1_VAL -1
      #define GAMMA_T_SPIN2_CO   0
      #define GAMMA_T_SPIN2_VAL -1
      #define GAMMA_T_SPIN3_CO   1
      #define GAMMA_T_SPIN3_VAL -1

      /* gamma_Z =
      *  0  0  I  0
      *  0  0  0 -I
      * -I  0  0  0
      *  0  I  0  0  
      */
      #define GAMMA_Z_SPIN0_CO   2
      #define GAMMA_Z_SPIN0_VAL -I
      #define GAMMA_Z_SPIN1_CO   3
      #define GAMMA_Z_SPIN1_VAL  I
      #define GAMMA_Z_SPIN2_CO   0
      #define GAMMA_Z_SPIN2_VAL  I
      #define GAMMA_Z_SPIN3_CO   1
      #define GAMMA_Z_SPIN3_VAL -I

      /* gamma_Y =
      *  0  0  0 -1
      *  0  0  1  0
      *  0  1  0  0
      * -1  0  0  0  
      */
      #define GAMMA_Y_SPIN0_CO   3
      #define GAMMA_Y_SPIN0_VAL  1
      #define GAMMA_Y_SPIN1_CO   2
      #define GAMMA_Y_SPIN1_VAL -1
      #define GAMMA_Y_SPIN2_CO   1
      #define GAMMA_Y_SPIN2_VAL -1
      #define GAMMA_Y_SPIN3_CO   0
      #define GAMMA_Y_SPIN3_VAL  1

      /* gamma_X =
      *  0  0  0  I
      *  0  0  I  0
      *  0 -I  0  0
      * -I  0  0  0  
      */
      #define GAMMA_X_SPIN0_CO   3
      #define GAMMA_X_SPIN0_VAL -I
      #define GAMMA_X_SPIN1_CO   2
      #define GAMMA_X_SPIN1_VAL -I
      #define GAMMA_X_SPIN2_CO   1
      #define GAMMA_X_SPIN2_VAL  I
      #define GAMMA_X_SPIN3_CO   0
      #define GAMMA_X_SPIN3_VAL  I
  /* ------------------------------------------------- */
    #else
      #ifdef BASIS2
        // Bais used for QCDSF lattices
        #define CLIFFORD_BASIS "BASIS2:QCDSF BASIS"
        /* gamma_T =
        *  0  0  1  0
        *  0  0  0  1
        *  1  0  0  0
        *  0  1  0  0  
        */
        #define GAMMA_T_SPIN0_CO   2
        #define GAMMA_T_SPIN0_VAL  1
        #define GAMMA_T_SPIN1_CO   3
        #define GAMMA_T_SPIN1_VAL  1
        #define GAMMA_T_SPIN2_CO   0
        #define GAMMA_T_SPIN2_VAL  1
        #define GAMMA_T_SPIN3_CO   1
        #define GAMMA_T_SPIN3_VAL  1

        /* gamma_Z =
        *  0  0  I  0
        *  0  0  0 -I
        * -I  0  0  0
        *  0  I  0  0  
        */
        #define GAMMA_Z_SPIN0_CO   2
        #define GAMMA_Z_SPIN0_VAL  I
        #define GAMMA_Z_SPIN1_CO   3
        #define GAMMA_Z_SPIN1_VAL -I
        #define GAMMA_Z_SPIN2_CO   0
        #define GAMMA_Z_SPIN2_VAL -I
        #define GAMMA_Z_SPIN3_CO   1
        #define GAMMA_Z_SPIN3_VAL  I

        /* gamma_Y =
        *  0  0  0 -1
        *  0  0  1  0
        *  0  1  0  0
        * -1  0  0  0  
        */
        #define GAMMA_Y_SPIN0_CO   3
        #define GAMMA_Y_SPIN0_VAL -1
        #define GAMMA_Y_SPIN1_CO   2
        #define GAMMA_Y_SPIN1_VAL  1
        #define GAMMA_Y_SPIN2_CO   1
        #define GAMMA_Y_SPIN2_VAL  1
        #define GAMMA_Y_SPIN3_CO   0
        #define GAMMA_Y_SPIN3_VAL -1

        /* gamma_X =
        *  0  0  0  I
        *  0  0  I  0
        *  0 -I  0  0
        * -I  0  0  0  
        */
        #define GAMMA_X_SPIN0_CO   3
        #define GAMMA_X_SPIN0_VAL  I
        #define GAMMA_X_SPIN1_CO   2
        #define GAMMA_X_SPIN1_VAL  I
        #define GAMMA_X_SPIN2_CO   1
        #define GAMMA_X_SPIN2_VAL -I
        #define GAMMA_X_SPIN3_CO   0
        #define GAMMA_X_SPIN3_VAL -I
      #else
        #ifdef BASIS3
          // Basis used in the QOPQDP Code (by James Osborn/USQCD)
          #define CLIFFORD_BASIS "BASIS3:QOPQDP BASIS"
          /* gamma_T =
          *  0  0  1  0 
          *  0  0  0  1
          *  1  0  0  0
          *  0  1  0  0  
          */
          #define GAMMA_T_SPIN0_CO   2
          #define GAMMA_T_SPIN0_VAL  1
          #define GAMMA_T_SPIN1_CO   3
          #define GAMMA_T_SPIN1_VAL  1
          #define GAMMA_T_SPIN2_CO   0
          #define GAMMA_T_SPIN2_VAL  1
          #define GAMMA_T_SPIN3_CO   1
          #define GAMMA_T_SPIN3_VAL  1

          /* gamma_Z =
          *  0  0  0  I
          *  0  0  I  0
          *  0 -I  0  0
          * -I  0  0  0  
          */
          #define GAMMA_Z_SPIN0_CO   3
          #define GAMMA_Z_SPIN0_VAL  I
          #define GAMMA_Z_SPIN1_CO   2
          #define GAMMA_Z_SPIN1_VAL  I
          #define GAMMA_Z_SPIN2_CO   1
          #define GAMMA_Z_SPIN2_VAL -I
          #define GAMMA_Z_SPIN3_CO   0
          #define GAMMA_Z_SPIN3_VAL -I
          
          /* gamma_Y =
          *  0  0  0 -1
          *  0  0  1  0
          *  0  1  0  0
          * -1  0  0  0  
          */
          #define GAMMA_Y_SPIN0_CO   3
          #define GAMMA_Y_SPIN0_VAL -1
          #define GAMMA_Y_SPIN1_CO   2
          #define GAMMA_Y_SPIN1_VAL  1
          #define GAMMA_Y_SPIN2_CO   1
          #define GAMMA_Y_SPIN2_VAL  1
          #define GAMMA_Y_SPIN3_CO   0
          #define GAMMA_Y_SPIN3_VAL -1
          
          /* gamma_X =
          *  0  0  I  0
          *  0  0  0 -I
          * -I  0  0  0
          *  0  I  0  0  
          */
          #define GAMMA_X_SPIN0_CO   2
          #define GAMMA_X_SPIN0_VAL  I
          #define GAMMA_X_SPIN1_CO   3
          #define GAMMA_X_SPIN1_VAL -I
          #define GAMMA_X_SPIN2_CO   0
          #define GAMMA_X_SPIN2_VAL -I
          #define GAMMA_X_SPIN3_CO   1
          #define GAMMA_X_SPIN3_VAL  I
#else
    #ifdef BASIS4
      // tmLQCD BASIS with an addition change of sign in gamma5
      //applied changing the order of spin component in the rhs and solution
      #define CLIFFORD_BASIS "BASIS4: tmLQCD basis"
      /* gamma_T =
      *  0  0 -1  0 
      *  0  0  0 -1
      * -1  0  0  0
      *  0 -1  0  0  
      */
      #define GAMMA_T_SPIN0_CO   2
      #define GAMMA_T_SPIN0_VAL -1
      #define GAMMA_T_SPIN1_CO   3
      #define GAMMA_T_SPIN1_VAL -1
      #define GAMMA_T_SPIN2_CO   0
      #define GAMMA_T_SPIN2_VAL -1
      #define GAMMA_T_SPIN3_CO   1
      #define GAMMA_T_SPIN3_VAL -1

      /* gamma_Z =
      *  0  0 -I  0
      *  0  0  0  I
      *  I  0  0  0
      *  0 -I  0  0  
      */
      #define GAMMA_Z_SPIN0_CO   2
      #define GAMMA_Z_SPIN0_VAL -I
      #define GAMMA_Z_SPIN1_CO   3
      #define GAMMA_Z_SPIN1_VAL  I
      #define GAMMA_Z_SPIN2_CO   0
      #define GAMMA_Z_SPIN2_VAL  I
      #define GAMMA_Z_SPIN3_CO   1
      #define GAMMA_Z_SPIN3_VAL -I

      /* gamma_Y =
      *  0  0  0 -1
      *  0  0  1  0
      *  0  1  0  0
      * -1  0  0  0  
      */
      #define GAMMA_Y_SPIN0_CO   3
      #define GAMMA_Y_SPIN0_VAL -1
      #define GAMMA_Y_SPIN1_CO   2
      #define GAMMA_Y_SPIN1_VAL  1
      #define GAMMA_Y_SPIN2_CO   1
      #define GAMMA_Y_SPIN2_VAL  1
      #define GAMMA_Y_SPIN3_CO   0
      #define GAMMA_Y_SPIN3_VAL -1

      /* gamma_X =
      *  0  0  0  I
      *  0  0  I  0
      *  0 -I  0  0
      * -I  0  0  0  
      */
      #define GAMMA_X_SPIN0_CO   3
      #define GAMMA_X_SPIN0_VAL  I
      #define GAMMA_X_SPIN1_CO   2
      #define GAMMA_X_SPIN1_VAL  I
      #define GAMMA_X_SPIN2_CO   1
      #define GAMMA_X_SPIN2_VAL -I
      #define GAMMA_X_SPIN3_CO   0
      #define GAMMA_X_SPIN3_VAL -I
      /* ------------------------------------------------- */
          #endif
        #endif
      #endif
    #endif
  #endif
#endif
