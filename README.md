# Processing scripts for Miller et al. 2020
Repository of scripts and data for reproducing analysis in Miller et al. 2020

## Quickstart

Copy and paste the following code to recreate the paper's Conda environment and run all analyses:

```
git clone https://github.com/millerh1/Ewing-sarcoma-cell-of-origin-paper-2020.git
cd Ewing-sarcoma-cell-of-origin-paper-2020/
conda env create -f environment.yml
conda activate ewsPaperEnv
(Rscript generateFigures.R) |& tee generateFigures_logFile.txt
```

## Additional details

### Conda installation

Conda can be used to quickly and easily recreate the analysis environment from this paper. It must be installed and available to the `$PATH` variable. To do this, download [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) and follow the installation instructions. 

### Manipulation of analysis scripts

Open the R project file `ewsBioinformaticsPaperRepo.Rproj` in R-Studio to easily interact with the analysis scripts. 

### Citation and Issues

You are welcome to re-use the analysis code and helper functions from this paper, but please cite it appropriately. Also, any issues or bugs in the analysis script should be reported here and/or you can email the code maintainer at millerh1@livemail.uthscsa.edu.

