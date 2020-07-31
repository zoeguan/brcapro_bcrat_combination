
### folders

sim_code
- code to simulate families and run models

families
- model inputs for simulated families

model_results
- model outputs for simulated families

tables
- code to generate performance tables

### instructions

1) run sim_code/simFam.sh \
$ cd sim_code \
$ sh simFam.sh
- generates 100000 families (10 sets of 10000), saves model inputs to families/simFam_1_*.RData

2) run sim_code/run_penmod_sim_array.sh, sim_code/run_brcapro_sim_array.sh, sim_code/run_gail_sim.sh (use cluster for BRCAPRO and BRCAPRO+BCRAT (M)) \
$ cd sim_code \
$ sbatch --array=1-10 run_penmod_sim_array.sh \
$ sbatch --array=1-10 run_brcapro_sim_array.sh \
$ sh run_gail_sim.sh
- runs models on simulated families
- model outputs saved in model_results: run_penmod_sim_1_*.RData, run_brcapro_sim_1_*.RData, run_gail_sim_1_*.RData

3) run tables/combine_results.R
- calculates performance measures




