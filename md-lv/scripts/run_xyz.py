"""Example script that runs an MD simulation and dumps an XYZ file."""

import glob
import os
import pathlib
import shutil
import subprocess
import numpy as np

# list of potential scales that we are going to apply
vs = np.linspace(1, 11, 11)

mdlv_dir = pathlib.Path(os.environ["MDLV_DIR"])
data_dir = pathlib.Path(os.environ["MDLV_DATA_DIR"])

# fetch the example Hertzian starting configuration
input_data = (mdlv_dir / "prepared-data").glob("equil*hertz*.json")
output_dir = mdlv_dir / "xyz"
mdlv_bin = mdlv_dir / "target/release/md-lv"

# check that the mdlv binary is valid, if not look in path
# if that fails, raise an exception
if not mdlv_bin.is_file():
    if shutil.which("md-lv") is not None:
        mdlv_bin = "md-lv"
    else:
        raise shutil.ExecError

# run 
for i, f in enumerate(input_data):
    print(i)
    for v in vs:
        command = f'{mdlv_bin} --potential hertz --time 10.0 \
            --init-config {f} --temp 0.01 --vscale {v} --seed {2000+i} --dt 1e-3 --out-time 1e-1 \
            --dir {output_dir}'

        input = command.split()
        subprocess.run(input)