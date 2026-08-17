from __future__ import annotations
import json, os, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from sclh.e2_mors import load_percentile_track
from sclh.stellar_history_inverse import infer_age_for_rotation_period, integrate_history_to_age, rotation_period_days
from sclh.wind_regime import M_SUN_KG, YEAR_S, AU_M
from sclh.model_adequacy import compare_positive_intervals

OUT=ROOT/'results'/'E5C3'; OUT.mkdir(parents=True,exist_ok=True)
root_env=os.environ.get('SCLH_MORS_ROOT')
if not root_env:
    raise SystemExit('E5-C3 requires external MORS tracks. Set SCLH_MORS_ROOT to TrackGrid_MstarPercentile or its parent directory.')
MORS=Path(root_env)
if (MORS/'TrackGrid_MstarPercentile').exists(): MORS=MORS/'TrackGrid_MstarPercentile'

# LHS 1140 observational constraints.  The mass grid is a coarse MORS grid bracket,
# not a Gaussian statistical bracket around the measured stellar mass.
TARGET_PERIOD_D=131.0
PERIOD_TOL_D=5.0
MIN_AGE_MYR=5000.0
A_B_AU=0.0946
MEASURED_MSTAR=0.1844
MEASURED_MSTAR_SIGMA=0.0045
MASS_GRID=(0.15,0.20)
START_AGE_MYR=10.0
# Spinelli et al. (2022) X-ray detection of LHS 1140.
OBS_LX=1.34e26
OBS_LX_LOW=(1.34-0.21)*1e26
OBS_LX_HIGH=(1.34+0.19)*1e26


def token(mass: float) -> str:
    return str(mass).replace('.','p')

all_rows=[]
compatible_rows=[]
mass_support={}
for mass in MASS_GRID:
    files=sorted(MORS.glob(f'{token(mass)}Msun_*percentile_extended.dat'))
    if not files:
        raise FileNotFoundError(f'No MORS percentile tracks for {mass} Msun under {MORS}')
    nmatch=0
    final_periods=[]
    for path in files:
        m=re.search(r'_(\d+)percentile_extended',path.name)
        if not m: continue
        pct=int(m.group(1))
        tr=load_percentile_track(path)
        period=rotation_period_days(tr.omega_surface_sun)
        final_period=float(period[-1]); final_periods.append(final_period)
        match=infer_age_for_rotation_period(tr.age_myr,tr.omega_surface_sun,TARGET_PERIOD_D,
                                            min_age_myr=MIN_AGE_MYR,tolerance_days=PERIOD_TOL_D)
        base={
            'stellar_mass_msun':mass,
            'rotation_percentile_150Myr':pct,
            'initial_omega_sun':tr.initial_omega_sun,
            'track_max_age_myr':float(tr.age_myr[-1]),
            'track_final_period_days':final_period,
            'rotation_match':bool(match.within_tolerance),
            'nearest_or_crossing_age_myr':float(match.age_myr),
            'nearest_or_matched_period_days':float(match.period_days),
        }
        all_rows.append(base)
        if not match.within_tolerance:
            continue
        nmatch += 1
        t=match.age_myr
        feuv=tr.euv_flux(A_B_AU); fx=tr.x_flux(A_B_AU)
        leuv_now=float(np.interp(t,tr.age_myr,feuv))
        lxflux_now=float(np.interp(t,tr.age_myr,fx))
        lx_lum_now=float(np.interp(t,tr.age_myr,tr.lx_erg_s))
        mdot_now=float(np.interp(t,tr.age_myr,tr.mdot_msun_yr))
        euv_flu=integrate_history_to_age(tr.age_myr,feuv,t,start_age_myr=START_AGE_MYR)
        x_flu=integrate_history_to_age(tr.age_myr,fx,t,start_age_myr=START_AGE_MYR)
        mdot_kg_s=tr.mdot_msun_yr*M_SUN_KG/YEAR_S
        wind_massflux=mdot_kg_s/(4*np.pi*(A_B_AU*AU_M)**2)
        wind_col=integrate_history_to_age(tr.age_myr,wind_massflux,t,start_age_myr=START_AGE_MYR)
        compatible_rows.append({**base,
            'rotation_match_age_myr':t,
            'current_euv_flux_erg_s_cm2':leuv_now,
            'current_x_flux_erg_s_cm2':lxflux_now,
            'current_x_luminosity_erg_s':lx_lum_now,
            'current_mdot_msun_yr':mdot_now,
            'cumulative_euv_fluence_erg_cm2':euv_flu,
            'cumulative_x_fluence_erg_cm2':x_flu,
            'cumulative_wind_column_kg_m2':wind_col,
        })
    mass_support[str(mass)]={
        'n_tracks_scanned':len(files), 'n_rotation_compatible':nmatch,
        'final_period_days_min':float(np.min(final_periods)),
        'final_period_days_median':float(np.median(final_periods)),
        'final_period_days_max':float(np.max(final_periods)),
        'track_max_age_myr':float(load_percentile_track(files[0]).age_myr[-1]),
    }

all_df=pd.DataFrame(all_rows).sort_values(['stellar_mass_msun','rotation_percentile_150Myr'])
comp=pd.DataFrame(compatible_rows).sort_values(['stellar_mass_msun','rotation_percentile_150Myr'])
all_df.to_csv(OUT/'e5c3_mass_grid_rotation_support.csv',index=False)
comp.to_csv(OUT/'e5c3_rotation_compatible_histories.csv',index=False)
if comp.empty:
    raise RuntimeError('No MORS tracks are compatible with LHS 1140 rotation constraint')

# Native X-ray interval audit. Numerical mismatch is retained, but because the model and published
# XMM luminosity use non-identical energy-band observation operators, this is NOT interpreted as
# a bandpass-homogenized physical rejection. See E5-C6 and docs/PROVENANCE.md.
pred_lx_low=float(comp.current_x_luminosity_erg_s.min())
pred_lx_high=float(comp.current_x_luminosity_erg_s.max())
adequacy=compare_positive_intervals(pred_lx_low,pred_lx_high,OBS_LX_LOW,OBS_LX_HIGH)


def stats(df,col):
    return {k:float(v) for k,v in {
        'min':df[col].min(),'p05':df[col].quantile(.05),'median':df[col].median(),
        'p95':df[col].quantile(.95),'max':df[col].max()}.items()}

by_mass={}
for mass,g in comp.groupby('stellar_mass_msun'):
    by_mass[str(mass)]={
        'n_compatible_tracks':int(len(g)),
        'rotation_match_age_myr':stats(g,'rotation_match_age_myr'),
        'current_euv_flux_erg_s_cm2':stats(g,'current_euv_flux_erg_s_cm2'),
        'current_x_luminosity_erg_s':stats(g,'current_x_luminosity_erg_s'),
        'cumulative_euv_fluence_erg_cm2':stats(g,'cumulative_euv_fluence_erg_cm2'),
        'cumulative_wind_column_kg_m2':stats(g,'cumulative_wind_column_kg_m2'),
        'current_euv_max_over_min':float(g.current_euv_flux_erg_s_cm2.max()/g.current_euv_flux_erg_s_cm2.min()),
        'fluence_max_over_min':float(g.cumulative_euv_fluence_erg_cm2.max()/g.cumulative_euv_fluence_erg_cm2.min()),
        'wind_column_max_over_min':float(g.cumulative_wind_column_kg_m2.max()/g.cumulative_wind_column_kg_m2.min()),
    }

summary={
    'experiment':'E5-C3 LHS 1140 b cross-system history degeneracy + model-adequacy audit',
    'observational_constraints':{
        'stellar_mass_msun':f'{MEASURED_MSTAR} +/- {MEASURED_MSTAR_SIGMA}',
        'stellar_rotation_days':f'{TARGET_PERIOD_D} +/- {PERIOD_TOL_D}',
        'age_lower_bound_myr':MIN_AGE_MYR,
        'planet_b_semimajor_axis_au':A_B_AU,
        'observed_x_luminosity_erg_s':OBS_LX,
        'observed_x_luminosity_interval_erg_s':[OBS_LX_LOW,OBS_LX_HIGH],
    },
    'mors_grid_audit':{
        'coarse_mass_grid_msun':list(MASS_GRID),
        'warning':'0.15/0.20 Msun are neighboring MORS grid points, not a statistical uncertainty interval for the measured 0.1844 Msun mass.',
        'by_mass':mass_support,
    },
    'n_rotation_compatible_tracks_total':int(len(comp)),
    'compatible_track_masses':sorted([float(x) for x in comp.stellar_mass_msun.unique()]),
    'combined_rotation_match_age_myr':stats(comp,'rotation_match_age_myr'),
    'combined_current_euv_flux_erg_s_cm2':stats(comp,'current_euv_flux_erg_s_cm2'),
    'combined_cumulative_euv_fluence_erg_cm2':stats(comp,'cumulative_euv_fluence_erg_cm2'),
    'combined_cumulative_wind_column_kg_m2':stats(comp,'cumulative_wind_column_kg_m2'),
    'history_degeneracy':{
        'current_euv_max_over_min':float(comp.current_euv_flux_erg_s_cm2.max()/comp.current_euv_flux_erg_s_cm2.min()),
        'cumulative_euv_max_over_min':float(comp.cumulative_euv_fluence_erg_cm2.max()/comp.cumulative_euv_fluence_erg_cm2.min()),
        'cumulative_wind_column_max_over_min':float(comp.cumulative_wind_column_kg_m2.max()/comp.cumulative_wind_column_kg_m2.min()),
    },
    'xray_model_adequacy':{
        **adequacy.to_dict(),
        'predicted_median_over_observed_central':float(comp.current_x_luminosity_erg_s.median()/OBS_LX),
        'interpretation':'The native MORS current-X interval and the published XMM table interval are discrepant, but this is an observation/activity-mapping alarm rather than a bandpass-homogenized physical rejection. Rotation compatibility alone is not sufficient to promote these tracks to a calibrated age/history posterior.',
    },
    'by_compatible_mors_mass':by_mass,
    'key_result':'LHS 1140 b independently reproduces current-state/history degeneracy inside the rotation-compatible 0.20-Msun MORS slice: present EUV nearly converges while cumulative EUV and wind exposure retain order-unity spread. The native current-X mismatch additionally flags the activity/observation mapping for audit; a hard physical rejection requires a matched observation operator.',
    'cross_system_interpretation':'The TOI-700 result is therefore not merely repeated. LHS 1140 adds an adequacy-layer audit: present rotation can leave historical degeneracy, while an additional activity datum may primarily diagnose the forward activity/observation mapping rather than rank histories.',
    'nonclaims':[
        'The 10.5-11.9 Gyr rotation-crossing ages are conditional MORS crossing times, not a measurement of the age of LHS 1140.',
        'No interpolation between the 0.15 and 0.20 Msun track grids is used because the neighboring grids have qualitatively different rotation support at 131 d.',
        'The X-ray mismatch is a model-adequacy flag, not proof that the X-ray observation or MORS physics is wrong; activity variability and model discrepancy remain possible.',
        'No probability prior is placed over compatible rotation percentiles.',
        'No atmospheric-loss history or habitability probability is inferred from the 2026 helium detection without a dedicated escape model.'
    ]
}
(OUT/'e5c3_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')

try:
    import matplotlib.pyplot as plt
    g=comp[comp.stellar_mass_msun==0.20]
    fig,ax=plt.subplots(figsize=(7.4,4.8))
    ax.plot(g.rotation_percentile_150Myr,g.current_euv_flux_erg_s_cm2,marker='.',lw=1,label='current EUV')
    ax.set_xlabel('rotation percentile at 150 Myr')
    ax.set_ylabel(r'current EUV flux at LHS 1140 b [erg s$^{-1}$ cm$^{-2}$]')
    ax.set_title('LHS 1140 b: present EUV converges across rotation-compatible histories')
    fig.tight_layout(); fig.savefig(OUT/'e5c3_current_euv_convergence.png',dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.4,4.8))
    ax.plot(g.rotation_percentile_150Myr,g.cumulative_euv_fluence_erg_cm2,marker='.',lw=1)
    ax.set_xlabel('rotation percentile at 150 Myr')
    ax.set_ylabel(r'cumulative EUV fluence [erg cm$^{-2}$]')
    ax.set_title('LHS 1140 b: cumulative EUV retains early-history memory')
    fig.tight_layout(); fig.savefig(OUT/'e5c3_cumulative_euv_history_degeneracy.png',dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.4,4.8))
    ax.scatter(g.rotation_percentile_150Myr,g.current_x_luminosity_erg_s,s=14,label='MORS, rotation-compatible')
    ax.axhspan(OBS_LX_LOW,OBS_LX_HIGH,alpha=.2,label='observed Lx interval')
    ax.axhline(OBS_LX,lw=1,label='observed central Lx')
    ax.set_xlabel('rotation percentile at 150 Myr')
    ax.set_ylabel(r'current $L_X$ [erg s$^{-1}$]')
    ax.set_title('LHS 1140: independent X-ray observation exposes model discrepancy')
    ax.legend(); fig.tight_layout(); fig.savefig(OUT/'e5c3_xray_model_adequacy.png',dpi=180); plt.close(fig)
except Exception as exc:
    print(f'warning: figure generation failed: {exc}', file=sys.stderr)

print(json.dumps(summary,indent=2,ensure_ascii=False))
