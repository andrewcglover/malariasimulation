/*
 * adult_mosquito_eqs.h
 *
 *  Created on: 17 Mar 2021
 *      Author: gc1610
 */

#ifndef SRC_ADULT_MOSQUITO_EQS_H_
#define SRC_ADULT_MOSQUITO_EQS_H_

#include <vector>
#include "aquatic_mosquito_eqs.h"

/*
 * Index helpers for the 2-D ATN-exposure state vector (0-based).
 *
 * Layout:
 *   [0..2]                                     aquatic E/L/P
 *   [3 .. 3+deltaqp1-1]                        Sv[0..deltaq]
 *   [3+deltaqp1 .. 3+deltaqp1*(1+spor_len)-1]  Ev[q,j] row-major (q outer)
 *   [3+deltaqp1*(1+spor_len) ..]               Iv[0..deltaq]
 *
 * q=0 is the unexposed (baseline) compartment; q=1..deltaq are ATN-exposed.
 * j=0..spor_len-1 are Erlang EIP stages.
 */
inline size_t sv_idx(size_t q) {
    return 3u + q;
}
inline size_t ev_idx(size_t q, size_t j, size_t deltaqp1, size_t spor_len) {
    return 3u + deltaqp1 + q * spor_len + j;
}
inline size_t iv_idx(size_t q, size_t deltaqp1, size_t spor_len) {
    return 3u + deltaqp1 * (1u + spor_len) + q;
}
inline size_t adult_state_size(size_t deltaqp1, size_t spor_len) {
    return 3u + deltaqp1 * (2u + spor_len);
}

/*
 * AdultMosquitoModel
 *
 * Structural constants (deltaq, spor_len) are fixed at construction.
 *
 * Per-step scalars/vectors are updated each timestep by
 * adult_mosquito_model_update(). Conventions:
 *   av_da      = a * delta_atn  (ATN exposure rate, day^-1; precomputed in R)
 *   delta_atn  = exposure probability per bite (fraction; needed for (1-da) terms)
 *   dn_atn     = extra mortality fraction on ATN exposure
 *   Lambda0_t  = coverage-weighted reduced baseline FOI (day^-1; from R kernels)
 *   Lambda_i   = per-compartment raw FOI (R applies Hill decay; C++ applies (1-da))
 *                Lambda_i[0] == foim always (baseline carries no Hill decay)
 *   rho_i      = per-compartment EIP stage rate (rho_i[0] unused; C++ uses scalar rho)
 *   B_post     = per-stage post-infection blocking probability (length spor_len)
 */
struct AdultMosquitoModel {
    AquaticMosquitoModel growth_model;
    // structural constants (fixed at construction)
    size_t deltaq;
    size_t deltaqp1;  // = deltaq + 1
    size_t spor_len;
    double kappa;     // waning rate out of each exposure compartment (day^-1); default 1.0
    double rho;       // = spor_len / dem  (baseline EIP stage rate)
    // per-step scalars
    double mu;
    double foim;      // baseline FOI (Lambda in v3)
    double av_da;     // a * delta_atn  (ATN exposure rate, day^-1)
    double delta_atn; // ATN exposure probability per bite (fraction, 0–1)
    double dn_atn;    // extra mortality fraction from ATN exposure
    double Lambda0_t; // coverage-weighted reduced baseline FOI
    // per-step vectors
    std::vector<double> Lambda_i;  // length deltaqp1
    std::vector<double> rho_i;     // length deltaqp1
    std::vector<double> B_post;    // length spor_len
    AdultMosquitoModel(
        AquaticMosquitoModel growth_model,
        size_t deltaq,
        size_t spor_len,
        double mu,
        double dem,
        double kappa,
        double foim
    );
};

// create a system of equations for the solver
integration_function_t create_eqs(AdultMosquitoModel& model);

#endif /* SRC_ADULT_MOSQUITO_EQS_H_ */
