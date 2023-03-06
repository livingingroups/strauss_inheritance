### Demographic turnover can be a leading driver of hierarchy dynamics, and social inheritance modifies its effects — Data and code 
#### By Eli Strauss, March 6th 2023

<!-- Reminder to update link once the paper is published -->
This repository contains scripts and data for conducting the analysis for [this paper](https://royalsocietypublishing.org/journal/rstb). 
For questions or comments, please contact straussed at gmail.com or submit an issue. 

----

### File description and instructions

Files are numbered and listed in the order that they are used for analysis.  

**0_create_ranks** - directory used for creating ranks used in subsequent analyses  
|     **1.hyena_intx_data.RData** - Interaction data used for inferring ranks  
|     **2.get_ranks_dynarank.R** - Script for inferring ranks. Note that this takes a long time to run  
|     **female_ranks.RData** - Ranks output by prior script  
**00_pull_data.R** - Script for assembling relevant data from the Mara Hyena Project database and combinng it with rank data. This script interacts with a private database and is thus not usable.  
**0_hyena_data.RData** - Complete dataset used in the remainder of analyses.**<-----I recommend starting analysis here**  
**1_define_functions.R** - Script defining functions used in simulations.   
**2_ars.R** - Analysis of annual reproductive success in hyenas, used to parameterize simulations
**2_repro_function.RData** - Function relating rank to reproductive success, created in previous script.   
**3_markov_model_rank_hyenas.R** - Script analyzing hierarchy dynamics in empirical hyena data.  
**4_simulate_societies.R** - Script that simulates societies that differ in rank inheritance and rank-related reproductive success. 
**4_rank_data_simulated.Rdata** - Dynamic hierarchies created by prior script. 
**5_markov_model_rank_simulated.R** - Script analyzing hierarchy dynamics in simulated data.  
**5_sim_markov_chains.Rdata** - Markov models, transition matrices, and plotting objects created by prior script.
**6_sim_plot_queens_inequality.R** - Models including queenship, plots of simulated data, analysis of inequality in reproductive success. 
