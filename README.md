# Terminal stuff

- Install R: I guess you have admin rights right? 

so on Ubuntu/Debian:

sudo apt update && sudo apt install -y r-base build-essential

- on fedora or so:

sudo dnf install -y R gcc-c++ make

- or we do:

cd ~

curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

./bin/micromamba create -y -p ~/ddm-env -c conda-forge r-base r-dplyr r-posterior cxx-compiler make
eval "$(./bin/micromamba shell hook -s bash)"

micromamba activate ~/ddm-env

which Rscript g++ make       
(then they should all have a path)

- get the files:

git clone https://github.com/annanerum/ddm7.git

cd ddm7

- start the run:

nohup Rscript 02_fit_ddm7_hier_linux.R . > log_ddm7.txt 2>&1 &

- logging:

tail -f log_ddm7.txt

# To check 

data_sim_svhier --> Daten aus Simulation --> put inside the BASE_DIR_RSTUDIO in 02_fit_ddm7_hier_linux.R 

*see* 

DATA_DIR  <- file.path(BASE_DIR, "data_sim_svhier")

**important**
BASE_DIR_RSTUDIO  <- "~/svhier" --> change to where everything is saved 

---

if CmdStan cannot be found / installed --> install manually in the terminal of RStudio by: 


install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

library(cmdstanr)

---

check_cmdstan_toolchain(fix = TRUE)   # check for C++ compiler and make


install_cmdstan(cores = 4)            # downloads and builds CmdStan


cmdstan_version()                     # should print a version number, if not sth is wrong 



