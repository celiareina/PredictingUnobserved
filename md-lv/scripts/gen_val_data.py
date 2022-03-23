"""Generate validation dataset for Hertzian and LJ systems"""

import os
import pathlib
import subprocess


data_dir = pathlib.Path(os.environ["MDLV_DATA_DIR"])

# configure parameters

lj_input = data_dir / "input-data/equil_n-10_na-4_l-3.015113445777636_dt-1e-3_visc-5_seed-203_phi-1.1000_pot-lj_rs-1_vs-0.01.json"
hertz_input = data_dir / "input-data/equil_n-10_na-5_l-3.1622776601683795_dt-1e-3_visc-5_seed-1_phi-1.0000_pot-hertz_rs-1_vs-1.json"