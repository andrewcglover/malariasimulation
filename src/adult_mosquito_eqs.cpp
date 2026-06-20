/*
 * adult_mosquito_eqs.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "adult_mosquito_eqs.h"

AdultMosquitoModel::AdultMosquitoModel(
    AquaticMosquitoModel growth_model,
    size_t deltaq,
    size_t spor_len,
    double mu,
    double dem,
    double kappa,
    double foim
    ) : growth_model(growth_model),
        deltaq(deltaq),
        deltaqp1(deltaq + 1u),
        spor_len(spor_len),
        kappa(kappa),
        rho(static_cast<double>(spor_len) / dem),
        mu(mu),
        foim(foim),
        av_da(0.0),
        delta_atn(0.0),
        dn_atn(0.0),
        Lambda0_t(foim),
        Lambda_i(deltaq + 1u, foim),
        rho_i(deltaq + 1u, static_cast<double>(spor_len) / dem),
        B_post(spor_len, 0.0)
{}

integration_function_t create_eqs(AdultMosquitoModel& model) {
    auto growth_eqs = create_eqs(model.growth_model);
    return [&model, growth_eqs](const state_t& x, state_t& dxdt, double t) {
        const size_t dp1 = model.deltaqp1;
        const size_t dq  = model.deltaq;
        const size_t sl  = model.spor_len;

        // Non-negativity floor for adult compartment reads.
        // DOPRI5 with high a_tol can overshoot near-zero states to slightly negative
        // values; a negative compartment feeds back into rate terms with the wrong sign
        // and can NaN-cascade. Clamping reads to zero breaks that cascade.
        // This is a robustness guard, not a speed fix (the kink at zero may cost
        // a few extra steps near a crossing, but prevents NaN propagation entirely).
        // Aquatic states (x[0..2]) are handled by the aquatic ODE; we clamp only
        // adult indices (3+). The nn lambda is used for all adult x[] reads below.
        auto nn = [&x](size_t i) -> double { return x[i] < 0.0 ? 0.0 : x[i]; };

        // --- total_M for aquatic sub-model ---
        double total_M_d = 0.0;
        for (size_t q = 0u; q < dp1; ++q) {
            total_M_d += nn(sv_idx(q));
            for (size_t j = 0u; j < sl; ++j)
                total_M_d += nn(ev_idx(q, j, dp1, sl));
            total_M_d += nn(iv_idx(q, dp1, sl));
        }
        model.growth_model.total_M = static_cast<size_t>(total_M_d);

        // run aquatic ODE (uses total_M set above)
        growth_eqs(x, dxdt, t);

        // --- pre-compute aggregates ---
        // Ecol[j] = sum_q Ev[q, j]  (column sum across ATN-exposure rows)
        std::vector<double> Ecol(sl, 0.0);
        for (size_t j = 0u; j < sl; ++j)
            for (size_t q = 0u; q < dp1; ++q)
                Ecol[j] += nn(ev_idx(q, j, dp1, sl));

        // Svtot = sum_q Sv[q]
        double Svtot = 0.0;
        for (size_t q = 0u; q < dp1; ++q)
            Svtot += nn(sv_idx(q));

        // Ivtot = sum_q Iv[q]
        double Ivtot = 0.0;
        for (size_t q = 0u; q < dp1; ++q)
            Ivtot += nn(iv_idx(q, dp1, sl));

        // B_post_Ecol_sum = sum_j B_post[j] * Ecol[j]  (used in dSv[1])
        // Ecol[j] already non-negative from nn reads above.
        double B_post_Ecol_sum = 0.0;
        for (size_t j = 0u; j < sl; ++j)
            B_post_Ecol_sum += model.B_post[j] * Ecol[j];

        // --- scalar shorthands ---
        const double mu      = model.mu;
        const double foim    = model.foim;       // Lambda (baseline FOI)
        const double av_da   = model.av_da;      // a * delta_atn  (rate, day^-1)
        const double da      = model.delta_atn;  // fraction, for (1-da) terms
        const double dn      = model.dn_atn;
        const double L0      = model.Lambda0_t;  // coverage-weighted reduced FOI
        const double kappa   = model.kappa;
        const double rho     = model.rho;

        const double av_da_s = av_da * (1.0 - dn);   // a*delta_atn*(1-dn)
        const double da_s    = da    * (1.0 - dn);    // delta_atn*(1-dn)

        // Clamp pupa count via nn for betaa (aquatic state P is index 2).
        const double betaa = 0.5 * nn(static_cast<size_t>(AquaticState::P))
                                 / model.growth_model.dp;

        // ==================== dSv ====================

        // q = 0 (baseline / unexposed)
        // dSv[0] = betaa + kappa*Sv[dq] - (av_da + (1-da)*Lambda_i[0] + mu)*Sv[0]
        dxdt[sv_idx(0)] =
            betaa
            + kappa * nn(sv_idx(dq))
            - (av_da + (1.0 - da) * model.Lambda_i[0] + mu) * nn(sv_idx(0));

        // q = 1 (first ATN-exposed compartment)
        // Inflow: newly exposed susceptibles minus those infected (Lambda0_t term)
        //         plus post-infection-blocked Ev mosquitoes returning here
        // v3: delta_atn*(av*(1-dn) - Lambda0_t)*Svtot + av*delta_atn*(1-dn)*sum(B_post*Ecol)
        if (dp1 > 1u) {
            dxdt[sv_idx(1)] =
                (av_da_s - da_s * L0) * Svtot
                + av_da_s * B_post_Ecol_sum
                - (av_da + (1.0 - da) * model.Lambda_i[1] + kappa + mu) * nn(sv_idx(1));
        }

        // q = 2..deltaq (conveyor: waning ATN effect)
        for (size_t q = 2u; q < dp1; ++q) {
            dxdt[sv_idx(q)] =
                kappa * nn(sv_idx(q - 1u))
                - (av_da + (1.0 - da) * model.Lambda_i[q] + kappa + mu) * nn(sv_idx(q));
        }

        // ==================== dEv ====================

        // (q=0, j=0)  baseline, first EIP stage
        dxdt[ev_idx(0u, 0u, dp1, sl)] =
            kappa * nn(ev_idx(dq, 0u, dp1, sl))
            + (1.0 - da) * model.Lambda_i[0] * nn(sv_idx(0u))
            - (av_da + rho + mu) * nn(ev_idx(0u, 0u, dp1, sl));

        // (q=1, j=0)  first exposed, first EIP stage
        // Inflow: surviving re-exposed non-blocked Ev + Lambda0_t term from Sv
        //         + non-ATN infections from Sv[1]
        if (dp1 > 1u) {
            dxdt[ev_idx(1u, 0u, dp1, sl)] =
                (1.0 - model.B_post[0]) * av_da_s * Ecol[0]
                + da_s * L0 * Svtot
                + (1.0 - da) * model.Lambda_i[1] * nn(sv_idx(1u))
                - (av_da + kappa + model.rho_i[1] + mu) * nn(ev_idx(1u, 0u, dp1, sl));
        }

        // (q=2..deltaq, j=0)
        for (size_t q = 2u; q < dp1; ++q) {
            dxdt[ev_idx(q, 0u, dp1, sl)] =
                kappa * nn(ev_idx(q - 1u, 0u, dp1, sl))
                + (1.0 - da) * model.Lambda_i[q] * nn(sv_idx(q))
                - (av_da + kappa + model.rho_i[q] + mu) * nn(ev_idx(q, 0u, dp1, sl));
        }

        // (q=0, j=1..spor_len-1)  baseline, later EIP stages
        for (size_t j = 1u; j < sl; ++j) {
            dxdt[ev_idx(0u, j, dp1, sl)] =
                kappa * nn(ev_idx(dq, j, dp1, sl))
                + rho * nn(ev_idx(0u, j - 1u, dp1, sl))
                - (av_da + rho + mu) * nn(ev_idx(0u, j, dp1, sl));
        }

        // (q=1, j=1..spor_len-1)  first exposed, later EIP stages
        if (dp1 > 1u) {
            for (size_t j = 1u; j < sl; ++j) {
                dxdt[ev_idx(1u, j, dp1, sl)] =
                    (1.0 - model.B_post[j]) * av_da_s * Ecol[j]
                    + model.rho_i[1] * nn(ev_idx(1u, j - 1u, dp1, sl))
                    - (av_da + kappa + model.rho_i[1] + mu) * nn(ev_idx(1u, j, dp1, sl));
            }
        }

        // (q=2..deltaq, j=1..spor_len-1)
        for (size_t q = 2u; q < dp1; ++q) {
            for (size_t j = 1u; j < sl; ++j) {
                dxdt[ev_idx(q, j, dp1, sl)] =
                    kappa * nn(ev_idx(q - 1u, j, dp1, sl))
                    + model.rho_i[q] * nn(ev_idx(q, j - 1u, dp1, sl))
                    - (av_da + kappa + model.rho_i[q] + mu) * nn(ev_idx(q, j, dp1, sl));
            }
        }

        // ==================== dIv ====================

        // q = 0  (baseline)
        dxdt[iv_idx(0u, dp1, sl)] =
            kappa * nn(iv_idx(dq, dp1, sl))
            + rho * nn(ev_idx(0u, sl - 1u, dp1, sl))
            - (av_da + mu) * nn(iv_idx(0u, dp1, sl));

        // q = 1  (first exposed)
        if (dp1 > 1u) {
            dxdt[iv_idx(1u, dp1, sl)] =
                av_da_s * Ivtot
                + model.rho_i[1] * nn(ev_idx(1u, sl - 1u, dp1, sl))
                - (av_da + kappa + mu) * nn(iv_idx(1u, dp1, sl));
        }

        // q = 2..deltaq
        for (size_t q = 2u; q < dp1; ++q) {
            dxdt[iv_idx(q, dp1, sl)] =
                kappa * nn(iv_idx(q - 1u, dp1, sl))
                + model.rho_i[q] * nn(ev_idx(q, sl - 1u, dp1, sl))
                - (av_da + kappa + mu) * nn(iv_idx(q, dp1, sl));
        }
    };
}

//[[Rcpp::export]]
Rcpp::XPtr<AdultMosquitoModel> create_adult_mosquito_model(
    Rcpp::XPtr<AquaticMosquitoModel> growth_model,
    double mu,
    int deltaq,
    int spor_len,
    double dem,
    double kappa,
    double foim
    ) {
    auto model = new AdultMosquitoModel(
        *growth_model,
        static_cast<size_t>(deltaq),
        static_cast<size_t>(spor_len),
        mu,
        dem,
        kappa,
        foim
    );
    return Rcpp::XPtr<AdultMosquitoModel>(model, true);
}

//[[Rcpp::export]]
void adult_mosquito_model_update(
    Rcpp::XPtr<AdultMosquitoModel> model,
    double mu,
    double foim,
    double av_da,
    double delta_atn,
    double dn_atn,
    double Lambda0_t,
    std::vector<double> Lambda_i,
    std::vector<double> rho_i,
    std::vector<double> B_post,
    double f
    ) {
    model->mu        = mu;
    model->foim      = foim;
    model->av_da     = av_da;
    model->delta_atn = delta_atn;
    model->dn_atn    = dn_atn;
    model->Lambda0_t = Lambda0_t;
    model->Lambda_i  = Lambda_i;
    model->rho_i     = rho_i;
    model->B_post    = B_post;
    model->growth_model.f   = f;
    model->growth_model.mum = mu;
}

//[[Rcpp::export]]
std::vector<double> adult_mosquito_model_save_state(
    Rcpp::XPtr<AdultMosquitoModel> model
    ) {
    // No deque to checkpoint; ODE solver owns the state vector.
    return {};
}

//[[Rcpp::export]]
void adult_mosquito_model_restore_state(
    Rcpp::XPtr<AdultMosquitoModel> model,
    std::vector<double> state
    ) {
    // Nothing to restore; see adult_mosquito_model_save_state.
}

//[[Rcpp::export]]
Rcpp::XPtr<Solver> create_adult_solver(
    Rcpp::XPtr<AdultMosquitoModel> model,
    std::vector<double> init,
    double r_tol,
    double a_tol,
    size_t max_steps
    ) {
    return Rcpp::XPtr<Solver>(
        new Solver(
            init,
            create_eqs(*model),
            r_tol,
            a_tol,
            max_steps,
            "adult mosquito"
        ),
        true
    );
}
