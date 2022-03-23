"""Generates equilibrated input states to for the Hertzian and LJ systems"""

import os
import pathlib
import shutil
import subprocess

mdlv_dir = pathlib.Path(os.environ["MDLV_DIR"])
data_dir = pathlib.Path(os.environ["MDLV_DATA_DIR"])

output_dir = data_dir / ""

# configure parameters

hertz_seed = 1

lj_prep_seed = 202
lj_finish_seed = 203

# loop over process calls

mdlv_bin = mdlv_dir / "target/release/md-lv"

# check that the mdlv binary is valid, if not look in path
# if that fails, raise an exception
if not mdlv_bin.is_file():
    if shutil.which("md-lv") is not None:
        mdlv_bin = "md-lv"
    else:
        raise shutil.ExecError

# hertzian case
command = f'{mdlv_bin} --potential hertz --time 1000000 --phi 1.0 --num 10 --num-a 5\
    --seed {hertz_seed} --dir {output_dir} equil-gd '

input = command.split()
subprocess.run(input)

# hertzian preprocess for lj case
command = f'{mdlv_bin} --potential hertz --time 1000000 --phi 1.1 --num 10 --num-a 4\
    --seed {lj_prep_seed} --dir {output_dir} equil-gd '

input = command.split()
subprocess.run(input)

# get previous output from glob
init_file = output_dir.glob(f"equil*seed-{lj_prep_seed}*pot-hertz*.json").__next__()

command = f'{mdlv_bin} --potential lj --time 1000000 \
    --init-config {init_file} --seed {lj_finish_seed} --dir {output_dir} equil-gd '

input = command.split()
subprocess.run(input)