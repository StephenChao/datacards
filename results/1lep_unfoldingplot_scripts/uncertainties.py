import math

total_unc_up = {
            'ggf1':{    'pdf':    1.0190, 
                        'scale':  1.0390, 
                        'alphaS': 1.0260, 
                        'isr':    1.0000, 
                        'fsr':    1.0000},
             'ggf2':{   'pdf':    1.0190, 
                        'scale':  1.0390, 
                        'alphaS': 1.0260, 
                        'isr':    1.0000, 
                        'fsr':    1.0000},
             'ggf3':{   'pdf':    1.0190, 
                        'scale':  1.0390, 
                        'alphaS': 1.0260, 
                        'isr':    1.0000, 
                        'fsr':    1.0000},
             'vbf':{    'pdf':    1.0210, 
                        'scale':  1.0040, 
                        'alphaS': 1.0000, 
                        'isr':    1.0000, 
                        'fsr':    1.0000}
             }

total_unc_down = {
            'ggf1':{    'pdf':    1.0190, 
                        'scale':  1.0390, 
                        'alphaS': 1.0260, 
                        'isr':    1.0000, 
                        'fsr':    1.0000},
             'ggf2':{   'pdf':    1.0190, 
                        'scale':  1.0390, 
                        'alphaS': 1.0260, 
                        'isr':    1.0000, 
                        'fsr':    1.0000},
             'ggf3':{   'pdf':    1.0190, 
                        'scale':  1.0390, 
                        'alphaS': 1.0260, 
                        'isr':    1.0000, 
                        'fsr':    1.0000},
             'vbf':{    'pdf':    1.0210, 
                        'scale':  0.9970, #0.9970/1.0040 
                        'alphaS': 1.0000, 
                        'isr':    1.0000, 
                        'fsr':    1.0000}
             }


acc_uncs = {
    "nom": {
	"ggf1": 0.059673,
	"ggf2": 0.018935,
	"ggf3": 0.003297,
	"vbf": 0.017611
    },
    "down": {
	"ggf1": 0.996240,
	"ggf2": 0.997987,
	"ggf3": 0.969822,
	"vbf": 0.979847
    },
    "up": {
	"ggf1": 1.004013,
	"ggf2": 1.011878,
	"ggf3": 1.020116,
	"vbf": 1.002024 
    }
}


COMPONENTS = ("pdf", "scale", "alphaS", "isr", "fsr")

def quad_sum_by_proc(kappa_dict):
    """
    kappa_dict[proc][source] = multiplicative factor (e.g. 1.019)
    Returns per-proc totals with deltas and combined kappas.
    """
    out = {}
    for proc, sources in kappa_dict.items():
        deltas = [abs(sources[k] - 1.0) for k in COMPONENTS if k in sources]
        delta = math.sqrt(sum(d*d for d in deltas))
        out[proc] = {
            "delta": delta,                # absolute fractional uncertainty
            "kappa_up": 1.0 + delta,       # convenient combined factors
            "kappa_down": 1.0 - delta,
        }
    return out

up_tot  = quad_sum_by_proc(total_unc_up)
down_tot = quad_sum_by_proc(total_unc_down)

# Pretty print
for proc in sorted(up_tot.keys()):
    du = up_tot[proc]["delta"]
    dd = down_tot[proc]["delta"]
    print(f"{proc:4s}  Δ_up={du:.6f}  κ_up={1+du:.6f}   Δ_down={dd:.6f}  κ_down={1-dd:.6f}")

def residual_nonacceptance(up_tot, down_tot, acc_uncs):
    """
    Compute the quadrature difference (total ⊖ acceptance) per process.

    up_tot[proc]["delta"]   = |kappa_up_total - 1|
    down_tot[proc]["delta"] = |1 - kappa_down_total|
    acc_uncs["up"][proc]    = acceptance kappa_up (>=1)
    acc_uncs["down"][proc]  = acceptance kappa_down (<=1)
    """
    out = {}
    for proc in acc_uncs["nom"].keys():
        # totals as deltas
        d_tot_up   = float(up_tot[proc]["delta"])
        d_tot_down = float(down_tot[proc]["delta"])

        # acceptance as deltas
        d_acc_up   = abs(acc_uncs["up"][proc]   - 1.0)
        d_acc_down = abs(1.0 - acc_uncs["down"][proc])

        # residual (protect against tiny negatives)
        d_res_up   = math.sqrt(max(0.0, d_tot_up**2   - d_acc_up**2))
        d_res_down = math.sqrt(max(0.0, d_tot_down**2 - d_acc_down**2))

        out[proc] = {
            "delta_up_total":   d_tot_up,
            "delta_down_total": d_tot_down,
            "delta_up_acc":     d_acc_up,
            "delta_down_acc":   d_acc_down,
            "delta_up_res":     d_res_up,
            "delta_down_res":   d_res_down,
            "kappa_up_res":     1.0 + d_res_up,
            "kappa_down_res":   1.0 - d_res_down,
        }
    return out

residuals = residual_nonacceptance(up_tot, down_tot, acc_uncs)

# Example printout
for proc, r in residuals.items():
    print(
        f"{proc:4s}  κ_res_up={r['kappa_up_res']:.6f}  κ_res_down={r['kappa_down_res']:.6f}  "
        f"(Δres_up={r['delta_up_res']:.6f}, Δres_down={r['delta_down_res']:.6f})"
    )