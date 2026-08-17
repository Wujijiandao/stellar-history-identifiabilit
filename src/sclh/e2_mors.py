"""Adapters for Johnstone-Bartel-Guedel rotation/XUV track assets (Zenodo 4266670)."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import re
import numpy as np

AU_CM = 1.495978707e13

@dataclass(frozen=True)
class MORSTrack:
    stellar_mass_msun: float
    percentile: float
    initial_omega_sun: float
    age_myr: np.ndarray
    omega_surface_sun: np.ndarray
    lx_erg_s: np.ndarray
    leuv_erg_s: np.ndarray
    mdot_msun_yr: np.ndarray
    dipole_g: np.ndarray

    def euv_flux(self, a_au: float) -> np.ndarray:
        r = float(a_au) * AU_CM
        return self.leuv_erg_s / (4.0 * np.pi * r * r)

    def x_flux(self, a_au: float) -> np.ndarray:
        r = float(a_au) * AU_CM
        return self.lx_erg_s / (4.0 * np.pi * r * r)


def _parse_header(path: Path) -> tuple[float, float, float]:
    mass = pct = omega0 = None
    with path.open(errors="replace") as f:
        for line in f:
            stripped=line.lstrip()
            if not stripped.startswith("#"):
                break
            if "Stellar mass =" in line:
                mass = float(re.search(r"Stellar mass\s*=\s*([0-9.]+)", line).group(1))
            elif "Initial rotation rate =" in line:
                omega0 = float(re.search(r"Initial rotation rate\s*=\s*([0-9.]+)", line).group(1))
            elif "Percentile =" in line:
                m = re.search(r"Percentile\s*=\s*([0-9.]+)", line)
                if m:
                    pct = float(m.group(1))
    if mass is None or pct is None or omega0 is None:
        raise ValueError(f"Could not parse MORS track header: {path}")
    return mass, pct, omega0


def load_percentile_track(path: str | Path) -> MORSTrack:
    path = Path(path)
    mass, pct, omega0 = _parse_header(path)
    a = np.loadtxt(path, comments="#")
    if a.ndim != 2 or a.shape[1] < 25:
        raise ValueError(f"Unexpected MORS track shape {a.shape}: {path}")
    # Published asset column map: 1 age, 2 surface Omega, 7 dipole, 8 wind mass loss,
    # 16 Lx, 22 total EUV luminosity (10--92 nm).
    return MORSTrack(
        stellar_mass_msun=mass,
        percentile=pct,
        initial_omega_sun=omega0,
        age_myr=a[:, 0],
        omega_surface_sun=a[:, 1],
        lx_erg_s=a[:, 15],
        leuv_erg_s=a[:, 21],
        mdot_msun_yr=a[:, 7],
        dipole_g=a[:, 6],
    )


def find_percentile_track(root: str | Path, stellar_mass_msun: float, percentile: int) -> Path:
    root = Path(root)
    mass_token = str(stellar_mass_msun).replace(".", "p")
    if mass_token.endswith("p0"):
        mass_token = mass_token[:-2] + "p0"
    pattern = f"{mass_token}Msun_{int(percentile)}percentile_extended.dat"
    hits = list(root.rglob(pattern))
    if len(hits) != 1:
        raise FileNotFoundError(f"Expected one {pattern} under {root}, found {len(hits)}")
    return hits[0]
