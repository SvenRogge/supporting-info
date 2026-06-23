#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script was used to run the MLIP MD jobs using the "mace-mp-0b3-medium" model, available at https://github.com/ACEsuit/mace-foundations
# This script will generate a single trajectory file, a checkpoint file and a log file containing information of the temperature and potential energy.
# You can run this script by calling "python stability_run_mace.py STRUCTURE_PATH TEMPERATURE OUTDIR", with STRUCTURE_PATH the path to the starting structure, TEMPERATURE the desired NVT temperature and OUTDIR the desired folder for the trajectory

import sys
import os
from ase.io import read, Trajectory
from ase.md.langevin import Langevin
from ase import units
from mace.calculators import MACECalculator

#setup Arguments
if len(sys.argv) < 4:
    print("Usage: python script.py <path_to_structure> <temperature> <outdir>")
    sys.exit(1)

STRUCTUREPATH = sys.argv[1]
temperature = float(sys.argv[2])
filename_base = os.path.basename(STRUCTUREPATH).rsplit('.', 1)[0]
output_dir =  sys.argv[3]
os.makedirs(output_dir, exist_ok=True) #creates folder if it does not yet exists.

traj_file = os.path.join(output_dir, f"{filename_base}.traj")
checkpoint_file = os.path.join(output_dir, f"checkpoint_{filename_base}.traj")
log_file = f"md_evolution_{filename_base}.log" # Clean filename
total_nsteps = 250000

#restart logic: start either from a fresh run, or a checkpoint file if it exists.
if os.path.exists(checkpoint_file):
    try:
        atoms = read(checkpoint_file, format='traj')
        if os.path.exists(traj_file):
            steps_completed = len(read(traj_file, index=':', format='traj')) * 10 #currently hardcoded for trajectories that are written each 10 steps
        else:
            steps_completed = 0
        print(f"Successfully loaded checkpoint at step {steps_completed}.")
        mode = 'a' 
    except Exception as e:
        print(f"Checkpoint corrupted: {e}. Starting fresh.")
        atoms = read(STRUCTUREPATH)
        steps_completed = 0
        mode = 'w'
else:
    print(f"--- Fresh start for {filename_base} ---")
    atoms = read(STRUCTUREPATH)
    steps_completed = 0
    mode = 'w'

#MLIP ASE Calculator
calc = MACECalculator(model_paths="mace-mp-0b3-medium.model", device="cuda", dispersion=True) #Disperion setting enabled.
atoms.calc = calc

#Langevin engine setup.
dyn = Langevin(atoms, timestep=0.5 * units.fs, temperature_K=temperature, friction=0.01 / units.fs)

#definition of the log file logger. This can be modified to different properties.
def write_info(dyn, atoms, steps_completed, checkpoint_file, log_filename):
    curr_step = dyn.get_number_of_steps() + steps_completed
    epot = atoms.get_potential_energy() / len(atoms)
    temp = atoms.get_temperature()
    
    log_line = f"Step: {curr_step:8d} | T: {temp:6.1f} K | Epot: {epot:12.4f} eV/at\n"
    
    with open(log_filename, "a") as f:
        if os.path.exists(log_filename) and os.stat(log_filename).st_size == 0:
            f.write("Simulation Log\n" + "="*45 + "\n")
        f.write(log_line)

    print(log_line.strip())
    tmp_path = checkpoint_file + ".tmp"
    atoms.write(tmp_path, format='traj')
    os.replace(tmp_path, checkpoint_file)

#log the initial state of the structure (reports 0K).
if steps_completed == 0:
    write_info(dyn, atoms, steps_completed, checkpoint_file, log_file)

#write to log file each 100 steps.
dyn.attach(lambda: write_info(dyn, atoms, steps_completed, checkpoint_file, log_file), interval=100)

#write trajectory file each 10 steps.
traj = Trajectory(traj_file, mode, atoms)
dyn.attach(traj.write, interval=10)

#defines steps that needs to be run, either from the checkpoint file, or a fresh run.
remaining_steps = total_nsteps - steps_completed
#runs the MD job.
if remaining_steps > 0:
    print(f"Running {remaining_steps} steps to reach total of {total_nsteps}...")
    dyn.run(remaining_steps)
    print("Simulation finished successfully.")
else:
    print("Target steps already reached or exceeded.")