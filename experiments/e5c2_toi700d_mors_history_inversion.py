from __future__ import annotations
import json,os,re,sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from sclh.e2_mors import load_percentile_track
from sclh.stellar_history_inverse import infer_age_for_rotation_period, integrate_history_to_age, pressure_support_window
from sclh.wind_regime import ram_pressure_from_mdot_npa, M_SUN_KG, YEAR_S, AU_M

OUT=ROOT/'results'/'E5C2'; OUT.mkdir(parents=True,exist_ok=True)
root_env=os.environ.get('SCLH_MORS_ROOT')
if not root_env:
    raise SystemExit('E5-C2 requires external MORS tracks. Set SCLH_MORS_ROOT to TrackGrid_MstarPercentile or its parent directory.')
MORS=Path(root_env)
# Accept either the direct grid or the extracted package root.
if (MORS/'TrackGrid_MstarPercentile').exists(): MORS=MORS/'TrackGrid_MstarPercentile'

# TOI-700 constraints used here are deliberately minimal and source-bounded.
TARGET_PERIOD_D=54.0; PERIOD_TOL_D=0.8; MIN_AGE_MYR=1500.0
A_D_AU=0.1633
# 0.415 +0.021/-0.020 Msun is bracketed by available MORS masses 0.40 and 0.45 Msun.
MASS_GRID=(0.40,0.45)
WIND_SPEEDS=(300.0,400.0,650.0,800.0)
# Exact Dong et al. 2020 Table-1 ram pressures reconstructed in E3-A.
P_DONG_MIN=68.86374908578253
P_DONG_MAX=169.42235348294153
START_AGE_MYR=10.0

rows=[]; support_rows=[]
for mass in MASS_GRID:
    token=str(mass).replace('.','p')
    files=sorted(MORS.glob(f'{token}Msun_*percentile_extended.dat'))
    for path in files:
        m=re.search(r'_(\d+)percentile_extended',path.name)
        if not m: continue
        pct=int(m.group(1))
        tr=load_percentile_track(path)
        match=infer_age_for_rotation_period(tr.age_myr,tr.omega_surface_sun,TARGET_PERIOD_D,min_age_myr=MIN_AGE_MYR,tolerance_days=PERIOD_TOL_D)
        if not match.within_tolerance: continue
        t=match.age_myr
        feuv=tr.euv_flux(A_D_AU)
        fx=tr.x_flux(A_D_AU)
        leuv_now=float(np.interp(t,tr.age_myr,feuv))
        lx_now=float(np.interp(t,tr.age_myr,fx))
        mdot_now=float(np.interp(t,tr.age_myr,tr.mdot_msun_yr))
        euv_flu=integrate_history_to_age(tr.age_myr,feuv,t,start_age_myr=START_AGE_MYR)
        x_flu=integrate_history_to_age(tr.age_myr,fx,t,start_age_myr=START_AGE_MYR)
        # Integrated spherical wind mass column at the planet, kg/m^2.
        mdot_kg_s=tr.mdot_msun_yr*M_SUN_KG/YEAR_S
        wind_massflux=mdot_kg_s/(4*np.pi*(A_D_AU*AU_M)**2)
        wind_col=integrate_history_to_age(tr.age_myr,wind_massflux,t,start_age_myr=START_AGE_MYR)
        rows.append({
            'stellar_mass_msun':mass,'rotation_percentile_150Myr':pct,'initial_omega_sun':tr.initial_omega_sun,
            'rotation_match_age_myr':t,'matched_period_days':match.period_days,
            'current_euv_flux_erg_s_cm2':leuv_now,'current_x_flux_erg_s_cm2':lx_now,
            'current_mdot_msun_yr':mdot_now,'cumulative_euv_fluence_erg_cm2':euv_flu,
            'cumulative_x_fluence_erg_cm2':x_flu,'cumulative_wind_column_kg_m2':wind_col,
        })
        for v in WIND_SPEEDS:
            p=ram_pressure_from_mdot_npa(tr.mdot_msun_yr,v,A_D_AU)
            pnow=float(np.interp(t,tr.age_myr,p))
            enter,exit_,frac=pressure_support_window(tr.age_myr,p,t,P_DONG_MIN,P_DONG_MAX,start_age_myr=START_AGE_MYR)
            support_rows.append({
                'stellar_mass_msun':mass,'rotation_percentile_150Myr':pct,'rotation_match_age_myr':t,
                'wind_speed_kms':v,'current_ram_pressure_nPa':pnow,
                'dong_support_entry_age_myr':enter,'dong_support_exit_age_myr':exit_,
                'fraction_10Myr_to_current_inside_Dong_pressure_support':frac,
                'dong_pressure_support_nPa_low':P_DONG_MIN,'dong_pressure_support_nPa_high':P_DONG_MAX,
            })

matches=pd.DataFrame(rows).sort_values(['stellar_mass_msun','rotation_percentile_150Myr'])
support=pd.DataFrame(support_rows).sort_values(['stellar_mass_msun','rotation_percentile_150Myr','wind_speed_kms'])
if matches.empty: raise RuntimeError('No MORS tracks matched TOI-700 rotation constraint')
matches.to_csv(OUT/'e5c2_mors_rotation_conditioned_histories.csv',index=False)
support.to_csv(OUT/'e5c2_dong_pressure_support_overlap.csv',index=False)

# Publication-oriented diagnostics: present-state convergence versus cumulative-history spread.
try:
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7.4,4.8))
    for mass,g in matches.groupby('stellar_mass_msun'):
        ax.plot(g.rotation_percentile_150Myr,g.cumulative_euv_fluence_erg_cm2,marker='.',lw=1,label=fr'{mass:.2f} $M_\odot$ MORS')
    ax.set_xlabel('rotation percentile at 150 Myr')
    ax.set_ylabel(r'cumulative EUV fluence at TOI-700 d [erg cm$^{-2}$]')
    ax.set_title('TOI-700 d: compatible present-day rotation histories retain different EUV doses')
    ax.legend(); fig.tight_layout(); fig.savefig(OUT/'e5c2_cumulative_euv_history_degeneracy.png',dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.4,4.8))
    for mass,g in matches.groupby('stellar_mass_msun'):
        ax.plot(g.rotation_percentile_150Myr,g.current_euv_flux_erg_s_cm2,marker='.',lw=1,label=fr'{mass:.2f} $M_\odot$ MORS')
    ax.set_xlabel('rotation percentile at 150 Myr')
    ax.set_ylabel(r'current EUV flux at TOI-700 d [erg s$^{-1}$ cm$^{-2}$]')
    ax.set_title('TOI-700 d: present EUV converges despite different early histories')
    ax.legend(); fig.tight_layout(); fig.savefig(OUT/'e5c2_current_euv_convergence.png',dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.4,4.8))
    for v,g in support.groupby('wind_speed_kms'):
        # median over the two bracketing stellar masses at each percentile when both exist
        q=g.groupby('rotation_percentile_150Myr',as_index=False)['fraction_10Myr_to_current_inside_Dong_pressure_support'].median()
        ax.plot(q.rotation_percentile_150Myr,q.fraction_10Myr_to_current_inside_Dong_pressure_support,label=f'{v:g} km/s')
    ax.set_xlabel('rotation percentile at 150 Myr')
    ax.set_ylabel('fraction of 10 Myr→current history inside Dong pressure support')
    ax.set_ylim(0,1)
    ax.set_title('Published TOI-700 d MHD anchors cover only part of admissible wind histories')
    ax.legend(); fig.tight_layout(); fig.savefig(OUT/'e5c2_dong_pressure_support_coverage.png',dpi=180); plt.close(fig)
except Exception as exc:
    print(f'warning: figure generation failed: {exc}', file=sys.stderr)

# Compact summaries by mass and combined model-mass bracket.
def stats(df,col):
    return {k:float(v) for k,v in {'min':df[col].min(),'p05':df[col].quantile(.05),'median':df[col].median(),'p95':df[col].quantile(.95),'max':df[col].max()}.items()}
summary_by_mass={}
for mass,g in matches.groupby('stellar_mass_msun'):
    summary_by_mass[str(mass)]={
        'n_compatible_tracks':int(len(g)),
        'rotation_match_age_myr':stats(g,'rotation_match_age_myr'),
        'current_euv_flux_erg_s_cm2':stats(g,'current_euv_flux_erg_s_cm2'),
        'cumulative_euv_fluence_erg_cm2':stats(g,'cumulative_euv_fluence_erg_cm2'),
        'cumulative_wind_column_kg_m2':stats(g,'cumulative_wind_column_kg_m2'),
        'fluence_max_over_min':float(g.cumulative_euv_fluence_erg_cm2.max()/g.cumulative_euv_fluence_erg_cm2.min()),
        'wind_column_max_over_min':float(g.cumulative_wind_column_kg_m2.max()/g.cumulative_wind_column_kg_m2.min()),
        'current_euv_max_over_min':float(g.current_euv_flux_erg_s_cm2.max()/g.current_euv_flux_erg_s_cm2.min()),
    }

pressure_summary={}
for v,g in support.groupby('wind_speed_kms'):
    pressure_summary[str(int(v))]={
        'current_ram_pressure_nPa':stats(g,'current_ram_pressure_nPa'),
        'support_fraction':stats(g,'fraction_10Myr_to_current_inside_Dong_pressure_support'),
        'fraction_of_compatible_tracks_currently_inside_Dong_pressure_support':float(((g.current_ram_pressure_nPa>=P_DONG_MIN)&(g.current_ram_pressure_nPa<=P_DONG_MAX)).mean()),
        'fraction_currently_below_Dong_pressure_support':float((g.current_ram_pressure_nPa<P_DONG_MIN).mean()),
    }

summary={
    'experiment':'E5-C2 TOI-700 d MORS history inversion + Dong MHD support-overlap audit',
    'observational_constraints':{
        'stellar_mass_msun':'0.415 +0.021/-0.020 (Nelson & Becker 2026 adopting Gilbert et al. 2023)',
        'stellar_rotation_days':'54 +/- 0.8 (Gilbert et al. 2020; adopted by Nelson & Becker 2026)',
        'age_lower_bound_myr':1500,
        'planet_d_semimajor_axis_au':A_D_AU,
    },
    'model_conditioning':{
        'mors_mass_bracket_msun':list(MASS_GRID),'rotation_percentiles_scanned':'all available integer percentile tracks in assets',
        'condition':'track must reach P_rot=54 d after 1.5 Gyr; no probability prior over compatible tracks is imposed',
    },
    'n_compatible_tracks_total':int(len(matches)),
    'combined_rotation_match_age_myr':stats(matches,'rotation_match_age_myr'),
    'combined_current_euv_flux_erg_s_cm2':stats(matches,'current_euv_flux_erg_s_cm2'),
    'combined_cumulative_euv_fluence_erg_cm2':stats(matches,'cumulative_euv_fluence_erg_cm2'),
    'combined_cumulative_wind_column_kg_m2':stats(matches,'cumulative_wind_column_kg_m2'),
    'combined_history_degeneracy':{
        'cumulative_euv_max_over_min':float(matches.cumulative_euv_fluence_erg_cm2.max()/matches.cumulative_euv_fluence_erg_cm2.min()),
        'cumulative_wind_column_max_over_min':float(matches.cumulative_wind_column_kg_m2.max()/matches.cumulative_wind_column_kg_m2.min()),
        'current_euv_max_over_min':float(matches.current_euv_flux_erg_s_cm2.max()/matches.current_euv_flux_erg_s_cm2.min()),
    },
    'by_mors_mass':summary_by_mass,
    'dong_mhd_pressure_support_nPa':[P_DONG_MIN,P_DONG_MAX],
    'wind_speed_sensitivity_kms':list(WIND_SPEEDS),
    'pressure_support_overlap':pressure_summary,
    'key_result':'Present-day slow rotation strongly narrows the current activity state but does not identify the early rotation percentile: compatible histories converge to ~54 d at different ages while retaining roughly factor-of-two differences in integrated EUV and wind exposure.',
    'mhd_bridge_result':'Dong et al. 2020 provides absolute escape anchors only at two ram-pressure states. MORS-conditioned histories spend only a subset of their lifetime inside that pressure interval, so the published MHD grid cannot be legitimately integrated over the full history without additional simulations.',
    'nonclaims':[
        'The inferred crossing age is conditional on MORS and is not an independent stellar-age measurement.',
        '0.40/0.45 Msun tracks bracket, rather than interpolate, the measured 0.415 Msun stellar mass.',
        'Wind speed is not supplied by these percentile track files; 300/400/650/800 km/s are explicit sensitivity values.',
        'No interpolation of Dong escape rates outside or between the two published wind states is used to emit a full atmospheric-loss history.',
        'No full TOI-700 d habitability probability is emitted.'
    ]
}
(OUT/'e5c2_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(summary,indent=2,ensure_ascii=False))
