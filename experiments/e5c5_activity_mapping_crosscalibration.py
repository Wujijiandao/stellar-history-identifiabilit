from __future__ import annotations
import json, math, sys
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from sclh.activity_mapping import wright2018_fully_convective_lx, stellar_lbol_erg_s

IN=ROOT/'results'/'E5C3'/'e5c3_rotation_compatible_histories.csv'
OUT=ROOT/'results'/'E5C5'; OUT.mkdir(parents=True,exist_ok=True)
if not IN.exists():
    raise SystemExit('E5-C5 requires frozen E5-C3 output.')

# LHS 1140 present-day constraints used only for the external activity-mapping audit.
MSTAR=0.1844
RSTAR=0.216
TEFF=3096.0
PROT=131.0
OBS_LX=1.34e26
OBS_LOW=1.13e26
OBS_HIGH=1.53e26

df=pd.read_csv(IN)
if df.empty:
    raise RuntimeError('E5-C3 compatible-history table is empty')

mors_low=float(df.current_x_luminosity_erg_s.min())
mors_med=float(df.current_x_luminosity_erg_s.median())
mors_high=float(df.current_x_luminosity_erg_s.max())
w18=wright2018_fully_convective_lx(MSTAR,RSTAR,TEFF,PROT)
lbol=stellar_lbol_erg_s(RSTAR,TEFF)

def dex_ratio(a,b):
    return float(math.log10(a/b))

def min_interval_gap_to_point(x,lo,hi):
    if lo <= x <= hi:
        return {'factor':1.0,'dex':0.0,'direction':'inside'}
    if x < lo:
        f=lo/x; return {'factor':float(f),'dex':float(math.log10(f)),'direction':'below'}
    f=x/hi; return {'factor':float(f),'dex':float(math.log10(f)),'direction':'above'}

rows=[
    {'mapping':'Spinelli2023 observed central','lx_erg_s':OBS_LX,'log10_lx_over_lbol':math.log10(OBS_LX/lbol),'role':'observation'},
    {'mapping':'MORS native rotation-compatible median','lx_erg_s':mors_med,'log10_lx_over_lbol':math.log10(mors_med/lbol),'role':'history-model native activity mapping'},
    {'mapping':'Wright2018 fully-convective central','lx_erg_s':w18.lx_erg_s,'log10_lx_over_lbol':w18.log10_lx_over_lbol,'role':'independent present-day empirical activity mapping'},
]
pd.DataFrame(rows).to_csv(OUT/'e5c5_activity_model_comparison.csv',index=False)

summary={
    'experiment':'E5-C5 independent fully-convective rotation-activity cross-calibration',
    'target':'LHS 1140',
    'input':'frozen E5-C3 rotation-compatible histories; no MORS raw tracks required',
    'stellar_inputs':{'mass_msun':MSTAR,'radius_rsun':RSTAR,'teff_k':TEFF,'rotation_period_days':PROT},
    'observed_xray':{'central_erg_s':OBS_LX,'interval_erg_s':[OBS_LOW,OBS_HIGH],'log10_lx_over_lbol_central':math.log10(OBS_LX/lbol)},
    'mors_native':{
        'range_erg_s':[mors_low,mors_high], 'median_erg_s':mors_med,
        'median_over_observed_central':mors_med/OBS_LX,
        'median_residual_dex':dex_ratio(mors_med,OBS_LX),
        'minimum_support_gap_factor':mors_low/OBS_HIGH,
        'minimum_support_gap_dex':dex_ratio(mors_low,OBS_HIGH),
    },
    'wright2018_central':{
        **w18.to_dict(),
        'prediction_over_observed_central':w18.lx_erg_s/OBS_LX,
        'residual_dex':dex_ratio(w18.lx_erg_s,OBS_LX),
        'distance_to_observed_interval':min_interval_gap_to_point(w18.lx_erg_s,OBS_LOW,OBS_HIGH),
    },
    'mapping_contrast':{
        'mors_median_over_wright2018_central':mors_med/w18.lx_erg_s,
        'mors_vs_wright_dex':dex_ratio(mors_med,w18.lx_erg_s),
    },
    'decision':'LOCALIZE_DISCREPANCY_TO_ACTIVITY_MAPPING_FIRST',
    'interpretation':(
        'The deterministic E5-C4 empty set is real for the native MORS present-day X-ray mapping, '
        'but the observed 131-d rotation and measured Lx are not intrinsically contradictory: an '
        'independent fully-convective empirical rotation-activity relation predicts an Lx close to the '
        'observed interval. The minimal model repair should therefore expand/audit the present-activity '
        'observation mapping before rejecting the rotation-history trajectories themselves.'
    ),
    'required_next_model':[
        'Represent stellar activity mapping as an explicit model-class index m_activity rather than a fixed deterministic MORS mapping.',
        'Carry intrinsic/epoch activity scatter explicitly in the present-day likelihood or support model; do not silently widen the measurement error.',
        'Keep the MORS current-state/history degeneracy result distinct from the current-activity mapping adequacy result.',
        'Do not multiply an entire historical XUV track by the present-day Lx calibration factor except as a labelled sensitivity experiment; a present normalization does not identify early activity history.'
    ],
    'nonclaims':[
        'The Wright2018 central relation is not a new age estimate and is not used as a zero-scatter likelihood.',
        'This experiment does not prove which activity prescription is physically correct.',
        'This experiment does not alter the frozen 1.567x cumulative-EUV and 1.543x wind-history spreads from E5-C3.',
        'No Gyr XUV history is reconstructed from the Wright2018 field-star relation.'
    ]
}
(OUT/'e5c5_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')

try:
    import matplotlib.pyplot as plt
    labels=['Observed','MORS native','Wright+2018']
    vals=[OBS_LX,mors_med,w18.lx_erg_s]
    fig,ax=plt.subplots(figsize=(7.3,4.7))
    ax.scatter(labels,vals,s=55)
    ax.axhspan(OBS_LOW,OBS_HIGH,alpha=.2,label='Spinelli et al. observational interval')
    ax.set_yscale('log')
    ax.set_ylabel(r'present $L_X$ [erg s$^{-1}$]')
    ax.set_title('LHS 1140: the E5-C4 mismatch is activity-mapping dependent')
    ax.legend()
    fig.tight_layout(); fig.savefig(OUT/'e5c5_lx_mapping_crosscalibration.png',dpi=180); plt.close(fig)
except Exception as exc:
    print(f'warning: figure generation failed: {exc}',file=sys.stderr)

print(json.dumps(summary,indent=2,ensure_ascii=False))
