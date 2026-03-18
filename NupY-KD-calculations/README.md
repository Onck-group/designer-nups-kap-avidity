# NupY-KD-calculations

This repository contains codes for simulating the binding affinity between NupY-variants and yeast NTR Kap95.

Dependencies:
- GROMACS 2018.4, installed on a Linux operating system
- Python > 2.8 with numpy, MDanalysis
- MATLAB R2018b

To reproduce results / run the flow, perform the following steps:

- extract the files from `reference_structures` such that individual NupY protein files are prepared.
- Use the `prepare_NupY_variants.sh` script, which should run _as is_ to prepare folders with the simulation files for all the NupY-variants.
- Enter the folders for individual protein simulations just created and run the `prep_replicas.sh` shell script to produce 20 replicas.
- Use the example shell scripts such as `js_PG1/2` to run the GROMACS simulations. Specific commands will differ per cluster type, operating system or load management software (slurm, PBS, etc.)

Post simulation, the following commands can be used to extract the KD-values:
`combine_trajectories.sh` to generate one aggregate simulation file (this approach works - on the basis of a statistical analysis)
- use the `do_kd_calculation_alt.sh` script to perform an automated calculation (based on Python) to generate the KD-values.

  
