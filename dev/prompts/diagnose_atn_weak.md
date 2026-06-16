# Task: diagnose why the pure-ATN arm under-performs

I've run `dev/atn_local_test.R` and looked at the plots. The **pure-ATN** arm is
averting fewer clinical cases and reducing PfPR by less than I'd expect from my
experience with the malariasimple version. The **Pyr-ATN** arm looks about right.

Please **diagnose only — don't change anything until we've confirmed the cause.**

## Why this pattern is informative

Pyr-ATN's impact is dominated by pyrethroid *killing* (`dn0 = only$dn0`), which
would mask any weakness in the ATN drug terms. The pure-ATN arm has no killing
(`dn0 = 0`, `dn0_atn` unset) — its entire effect is the drug's transmission
blocking. So a weak pure-ATN arm + a healthy Pyr-ATN arm points at the **ATN
drug-effect path**, not the net/killing path.

Note `p_atn = 1.0` is set in both ATN arms, so drug *exposure* is on — that's not
the cause.

## Prime suspect: `Lambda00sf` is never set

The script sets every ATN parameter **except `Lambda00sf`**. In the reference
(`dev/reference/malariasimple_aitn_deterministic_v3.R`) the magnitude of
**pre-infection blocking** is governed by `Lambda00sf`:

```
Lambda00   = Lambda * Lambda00sf
Lambda0_t  starts at Lambda00 just after distribution, relaxes back to Lambda
b_lab_pre  = (1 - Lambda0_t / Lambda) * Hill(...)
Lambda_i   = Lambda * (1 - b_field_pre)
```

At peak effect the pre-infection blocking is `(1 - Lambda00sf)`. If the
`get_parameters()` default is `Lambda00sf = 1`, then `(1 - 1) = 0` and
pre-infection blocking is **completely off** — exposed mosquitoes keep the
baseline force of infection. Pre-infection blocking is the ATN's dominant
mechanism, so losing it would make the pure-ATN arm under-avert while leaving
Pyr-ATN (killing-dominated) looking fine.

This would parallel the other no-op ATN defaults: neither the script nor the
malariasimple run scripts set `Lambda00sf` explicitly, so malariasimple must rely
on a working default `< 1`; if the fork's default is `1`, that difference is the
bug.

## Checks (report findings, don't fix yet)

1. Report the fork's `get_parameters()` default for `Lambda00sf`, and find the
   malariasimple default (in `dev/reference/...` or the original
   `get_parameters.R` if available). Flag if fork = 1 and malariasimple < 1.
2. Print `Lambda_i` across the exposure compartments during the pure-ATN run and
   compare it to the baseline `foim`. If `Lambda_i == foim` everywhere,
   pre-infection blocking is confirmed off.
3. The script also omits `m_bompard` / `k_bompard` — confirm the fork's defaults
   match malariasimple's (the lab-TRA → field-TBA transform depends on them).

## If checks 1–2 don't explain it

Fall back to a head-to-head against malariasimple under matched ATN parameters
and bisect: compare the computed kernels (`Lambda_i`, `rho_i`, `B_post`) for
identical inputs (watch for an off-by-one from 0-based C++ vs 1-based reference in
the `s = i - 0.5` term), then the ATN-exposure compartment trajectories
(`Sv`/`Ev`/`Iv` by exposure index). Report where they first diverge.
