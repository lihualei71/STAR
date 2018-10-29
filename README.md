# Paper Repository

This repository contains the code to exactly reproduce all results and figures in our paper: [STAR: A general interactive framework for FDR control under structural constraints](https://arxiv.org/pdf/1710.02776.pdf). 

## R files
The folder `R/` contains all R files:

- `generic_STAR.R` implements a generic version of STAR with conditional one-group model as the working model. The function `generic.STAR` allows any accumulation functions (h), masking functions (g), model specification, and constraints. It also allows to produce movies of the updating process;

- `STAR_convex.R`, `STAR_bump_hunting.R`, `STAR_tree.R`, `STAR_wavelet.R` and `STAR_DAG.R` implement STAR in specific settings. Each of them include at least three components: (1) the function to find all possible candidate hypotheses (or subsets of hypotheses) to reveal; (2) the function to update the masked set; (3) the function to compute the scores;

- `STAR_convex_expr.R`, `STAR_bump_hunting_expr.R`, `STAR_tree_expr.R`, `STAR_wavelet_expr.R`, `STAR_DAG_expr.R` and `STAR_factorial_expr.R` implement all experiments in the paper;

- `preview_fun_plot.R`, `STAR_convex_plot.R`, `STAR_bump_hunting_expr.R`, `STAR_tree_plot.R`, `STAR_wavelet_plot.R`, `STAR_DAG_plot.R` and `STAR_factorial_expr.R` produce all figures in the paper;

- `AdaPT.R`, `AdaPT_gam.R`, `Lynch.R` and `Ramdas.R` implement other methods that we compare STAR with;

- `FDR_power_plot.R`, `expr_fun.R` and `summarize_methods.R` give a nuch of helper functions. 

## Replicating the experiments
- `STAR_bump_hunting_expr.R`, `STAR_wavelet_expr.R` and `STAR_factorial_expr.R` can be executed in a laptop;

- `STAR_convex_expr.R`, `STAR_tree_expr.R` and `STAR_DAG_expr.R` must be executed on clusters. To do so, submit `gen_job_XXX.sh` to a cluster, respectively (with "sbatch" system; otherwise the bash files need to be modified correspondingly). 

- The final results to reproduce the figures for above cases are already contained in `data/` folder. If one wants to reproduce these results, run the commented code chunks in `STAR_convex_plot.R`, `STAR_tree_plot.R` and `STAR_DAG_plot.R`, started with "## Read Patch Data from Clusters". Then the .RData objects required for figures will be generated into `data/` folder.
=======
# STAR
Selectively Traversed Accumulation Rules

This is the repository for "STAR: A general interactive framework for FDR control under structural constraints". Several movies are included for illustration. (1) naiveconvex.gif and gamconvex.gif are movies for Section 4 (convex region detection). The former uses the canonical score and the latter uses the GAM-assisted score; (2) naivetree.gif and isotree.gif are movies for Section 5 (hierarchical testing). The former uses the canonical score and the latter uses the isotonic-regression-assisted score; (4) DAG_strong.gif and DAG_weak.gif are movies for Section 6 (testing on DAGs). The former considers the strong heredity principle and the latter considers the weak heredity principle. Both movies use the canonical score. 
