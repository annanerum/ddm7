# To check 

data_sim_svhier --> Daten aus Simulation --> put inside the BASE_DIR_RSTUDIO in 02_fit_ddm7_hier_linux.R 

*see* 

DATA_DIR  <- file.path(BASE_DIR, "data_sim_svhier")

**important**
BASE_DIR_RSTUDIO  <- "~/svhier" --> change to where everything is saved 

--

if CmdStan cannot be found / installed --> install manually in the terminal of RStudio by: 


install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

library(cmdstanr)

--

check_cmdstan_toolchain(fix = TRUE)   # check for C++ compiler and make


install_cmdstan(cores = 4)            # downloads and builds CmdStan


cmdstan_version()                     # should print a version number, if not sth is wrong 



