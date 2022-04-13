"""Generate validation dataset for Hertzian and LJ systems"""

import os
import pathlib
import shutil
import subprocess
import sys

import numpy as np

# configure parameters

idx = int(sys.argv[1]) + 4000

# Overlap function parameter list
As = np.linspace(0.1, 0.3, 15)
as_str = ",".join([str(a) for a in As])

# vscales
vs = -np.logspace(0, 1, 11) + 1.0
vs_str = ",".join([str(a) for a in vs])

mdlv_dir = pathlib.Path(os.environ["MDLV_DIR"])
data_dir = pathlib.Path(os.environ["MDLV_DATA_DIR"])
output_dir = data_dir / "pred"
lj_input = data_dir / "input-data/equil_n-10_na-4_l-3.015113445777636_dt-1e-3_visc-5_seed-203_phi-1.1000_pot-lj_rs-1_vs-0.01.json"
hertz_input = data_dir / "input-data/equil_n-10_na-5_l-3.1622776601683795_dt-1e-3_visc-5_seed-1_phi-1.0000_pot-hertz_rs-1_vs-1.json"

mdlv_bin = mdlv_dir / "target/release/md-lv"

# check that the mdlv binary is valid, if not look in path
# if that fails, raise an exception
if not mdlv_bin.is_file():
    if shutil.which("md-lv") is not None:
        mdlv_bin = "md-lv"
    else:
        raise shutil.ExecError

# lennard jones
command = f'/home/igraham/Documents/md-lv/target/release/md-lv --unwrap --potential lj \
    --init-config {lj_input} --temp 0.1 --seed {idx} --dt 1e-3 --out-time 1e-2 \
    --dir {output_dir}/lj --time 1.0 \
    gen-variant --realizations 100000 --del-var={vs_str} --calc-msd --calc-pos --calc-q={as_str}'

input = command.split()
subprocess.run(input)

command = f'/home/igraham/Documents/md-lv/target/release/md-lv --unwrap --potential hertz \
    --init-config {lj_input} --temp 0.1 --seed {idx} --dt 1e-3 --out-time 1e-2 \
    --dir {output_dir}/hertz --time 1.0 \
    gen-variant --realizations 100000 --del-var={vs_str} --calc-msd --calc-pos --calc-q={as_str}'

input = command.split()
subprocess.run(input)