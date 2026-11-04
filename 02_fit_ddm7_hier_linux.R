options(error = NULL)


BASE_DIR_RSTUDIO  <- "~/svhier"   # gucken ob das stimmt dann
REPS              <- 1:5
N_PERSON          <- 100         
ITER_WARMUP       <- 500
ITER_SAMPLING     <- 500
PRECISION         <- 1e-3         # if its not working we can take Stan default 1e-4
ADAPT_DELTA       <- 0.90
THREADS_PER_CHAIN <- NULL         # NULL = all cores / 4
TEST_MODE         <- FALSE


args     <- commandArgs(trailingOnly = TRUE)
BASE_DIR <- normalizePath(if (length(args) >= 1) args[1] else BASE_DIR_RSTUDIO, mustWork = TRUE)
if (length(args) >= 2) REPS <- as.integer(args[-1])



# packages + CmdStan 



user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
cat("R library:", user_lib, "\n")

repos <- c(stan = "https://stan-dev.r-universe.dev",   
           CRAN = "https://cloud.r-project.org")
for (pkg in c("cmdstanr", "posterior", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("installing", pkg, "...\n")
    install.packages(pkg, lib = user_lib, repos = repos,
                     Ncpus = max(1, parallel::detectCores() - 1))
  }
}
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
})

# CmdStan: needs >= 2.35 for the 7-parameter Wiener density
# otherwise we need to download the latest release from the official Stan GitHub page
cmdstan_ok <- function() {
  v <- tryCatch(cmdstan_version(error_on_NA = FALSE), error = function(e) NULL)
  !is.null(v) && !is.na(v) && utils::compareVersion(v, "2.35.0") >= 0
}
CMDSTAN_DIR <- file.path(BASE_DIR, "cmdstan")
if (!cmdstan_ok() && dir.exists(CMDSTAN_DIR)) {             # installed here before
  inst <- sort(list.dirs(CMDSTAN_DIR, recursive = FALSE), decreasing = TRUE)
  if (length(inst) > 0) try(set_cmdstan_path(inst[1]), silent = TRUE)
}
if (!cmdstan_ok()) {
  cat("\nCmdStan not found (or older than 2.35) --> installing into", CMDSTAN_DIR, "\n")
  check_cmdstan_toolchain(fix = TRUE)
  dir.create(CMDSTAN_DIR, showWarnings = FALSE)
  install_cmdstan(dir = CMDSTAN_DIR, cores = max(1, parallel::detectCores() - 1))
}
cat("CmdStan", cmdstan_version(), "at", cmdstan_path(), "\n\n")




DATA_DIR  <- file.path(BASE_DIR, "data_sim_svhier")
DDM5_DIR  <- file.path(BASE_DIR, "fits_ddm5_2conditions")
FIT_DIR   <- file.path(BASE_DIR, if (TEST_MODE) "fits_ddm7_hier_TEST" else "fits_ddm7_hier")
STAN_FILE <- file.path(BASE_DIR, "ddm7_hier.stan")
dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(STAN_FILE)) stop("Stan file not found: ", STAN_FILE)
missing_reps <- REPS[!file.exists(file.path(DATA_DIR, sprintf("rep%02d", REPS), "sim_data_full_ddm.csv"))]
if (length(missing_reps) > 0) stop("No data for reps: ", paste(missing_reps, collapse = ", "))

n_chains <- 4
n_cores  <- parallel::detectCores()
threads_per_chain <- if (is.null(THREADS_PER_CHAIN)) max(1, floor(n_cores / n_chains)) else THREADS_PER_CHAIN
cat(sprintf("BASE_DIR: %s\nreps: %s | %d persons | %d + %d iterations | precision %g\ncores: %d --> %d chains x %d threads\n",
            BASE_DIR, paste(REPS, collapse = " "), N_PERSON, ITER_WARMUP, ITER_SAMPLING,
            PRECISION, n_cores, n_chains, threads_per_chain))
cat("DDM5 summaries for inits:", if (dir.exists(DDM5_DIR)) "found" else "not found (default inits)", "\n")


mod <- tryCatch(cmdstan_model(STAN_FILE, cpp_options = list(stan_threads = TRUE)),
                error = function(e) { message("compile error: ", conditionMessage(e)); NULL })
HAS_PRECISION <- !is.null(mod)
if (HAS_PRECISION) {
  cat("\nCOMPILE CHECK: precision argument accepted\n\n")
} else {
  cat("\nCOMPILE CHECK: precision argument not accepted --> Stan default 1e-4 is used\n\n")
  code <- gsub(", st0[p], prec)", ", st0[p])", readLines(STAN_FILE), fixed = TRUE)
  f <- file.path(BASE_DIR, "ddm7_hier_noprec.stan"); writeLines(code, f)
  mod <- cmdstan_model(f, cpp_options = list(stan_threads = TRUE))
  PRECISION <- 1e-4
}

priors <- list(
  prior_beta_nu_intercept     = c(1.5, 1.0),
  prior_beta_nu_c             = c(0.0, 1.0),
  prior_beta_alpha            = c(log(1.5), 0.5),     # log scale
  prior_beta_tau              = c(0.25, 0.10),        # seconds (bounded 0.05 - 0.45)
  prior_beta_w                = c(0.0, 0.5),          # logit scale
  prior_beta_sv               = c(log(0.5), 0.5),     # log scale
  prior_sigma_nu_intercept    = c(0, 0.5),
  prior_sigma_nu_c            = c(0, 0.5),
  prior_sigma_alpha_intercept = c(0, 0.2),
  prior_sigma_tau_intercept   = c(0, 0.1),
  prior_sigma_w               = c(0, 0.5),
  prior_sigma_sv              = c(0, 0.3),
  # sw and st0 
  # beta_sw = logit of the fraction of each person's max sw (sw_max = 0.99*2*min(w,1-w))
  prior_beta_sw               = c(qlogis(0.15), 0.5),
  prior_beta_st0              = c(log(0.08), 0.5),    # log seconds
  prior_sigma_sw              = c(0, 0.07),
  prior_sigma_st0             = c(0, 0.07),
  guard_a_max = 5, guard_sv_max = 3, guard_v_max = 8, guard_w_min = 0.05, guard_st0_max = 0.5
)

make_data <- function(r, n_person) {
  dat <- read.csv(file.path(DATA_DIR, sprintf("rep%02d", r), "sim_data_full_ddm.csv"))
  all_subj  <- sort(unique(dat$subj))
  keep_subj <- seq_len(min(n_person, length(all_subj)))
  dat <- dat[dat$subj %in% all_subj[keep_subj], ]
  dat$pid <- match(dat$subj, sort(unique(dat$subj)))
  min_rt <- as.numeric(tapply(dat$rt, dat$pid, min))
  list(stan_data = c(list(N = nrow(dat), I = max(dat$pid), pid = dat$pid, rt = dat$rt,
                          acc = as.integer(dat$acc), condition = dat$condition,
                          min_rt = min_rt, grainsize = 1),
                     priors, list(precision_derivatives = PRECISION)),
       subj_ids = all_subj[keep_subj], keep_subj = keep_subj, min_rt = min_rt)
}

# inits from the DDM5 posterior means of the same rep or else defaults
make_init_fun <- function(r, d) {
  I <- d$stan_data$I; min_rt <- d$min_rt; keep <- d$keep_subj
  f5 <- file.path(DDM5_DIR, sprintf("rep%02d", r), "summary_ddm5.csv")
  if (file.exists(f5)) {
    s5 <- read.csv(f5)
    g  <- function(v) s5$mean[s5$variable == v]
    pv <- function(base) s5$mean[grepl(paste0("^", base, "\\["), s5$variable)][keep]
    zf <- function(x, mu, sig) pmin(pmax((x - mu) / max(sig, 0.05), -3), 3)
    b  <- list(nu1 = g("beta_nu_intercept"), nuc = g("beta_nu_c"), al = g("beta_alpha"),
               tau = g("beta_tau"), w = g("beta_w"), sv = g("beta_sv"),
               s_nu1 = g("sigma_nu_intercept"), s_nuc = g("sigma_nu_c"),
               s_al = g("sigma_alpha_intercept"), s_tau = g("sigma_tau_intercept"),
               s_w = g("sigma_w"), s_sv = g("sigma_sv"))
    z0 <- list(nu1 = zf(pv("nu_intercept"), b$nu1, b$s_nu1), nuc = zf(pv("nu_c"), b$nuc, b$s_nuc),
               al = zf(pv("alpha_intercept"), b$al, b$s_al), w = zf(qlogis(pv("w")), b$w, b$s_w),
               sv = zf(log(pv("sv")), b$sv, b$s_sv))
    t0_0 <- pmax(pmin(pv("t0") - 0.02, 0.9 * min_rt), 0.05)
    cat("inits from DDM5 rep", r, "\n")
  } else {
    b  <- list(nu1 = 1.5, nuc = 0.5, al = log(1.5), tau = 0.2, w = 0, sv = log(0.5),
               s_nu1 = 0.3, s_nuc = 0.2, s_al = 0.1, s_tau = 0.05, s_w = 0.2, s_sv = 0.1)
    z0 <- list(nu1 = rep(0, I), nuc = rep(0, I), al = rep(0, I), w = rep(0, I), sv = rep(0, I))
    t0_0 <- pmin(0.15, 0.8 * min_rt)
    cat("default inits (no DDM5 summary for rep", r, ")\n")
  }
  jit <- function(x, sd = 0.05) x + rnorm(length(x), 0, sd)
  function() list(
    beta_nu_intercept = jit(b$nu1), beta_nu_c = jit(b$nuc), beta_alpha = jit(b$al, 0.02),
    beta_tau = min(max(jit(b$tau, 0.01), 0.06), 0.44), beta_w = jit(b$w, 0.02), beta_sv = jit(b$sv),
    beta_sw = qlogis(0.15), beta_st0 = log(0.06),
    sigma_nu_intercept = max(b$s_nu1, 0.05), sigma_nu_c = max(b$s_nuc, 0.05),
    sigma_alpha_intercept = max(b$s_al, 0.03), sigma_tau_intercept = max(b$s_tau, 0.02),
    sigma_w = max(b$s_w, 0.05), sigma_sv = max(b$s_sv, 0.05),
    sigma_sw = 0.03, sigma_st0 = 0.03,
    z_nu_intercept = jit(z0$nu1), z_nu_c = jit(z0$nuc), z_alpha_intercept = jit(z0$al),
    z_w = jit(z0$w), z_sv = jit(z0$sv),
    z_sw = rnorm(I, 0, 0.1), z_st0 = rnorm(I, 0, 0.1),
    t0 = pmax(pmin(t0_0 + rnorm(I, 0, 0.005), 0.95 * min_rt), 0.02)
  )
}



if (TEST_MODE) {
  d <- make_data(REPS[1], N_PERSON)
  precs <- if (HAS_PRECISION) unique(c(1e-4, PRECISION)) else 1e-4
  res <- do.call(rbind, lapply(precs, function(pr) {
    d$stan_data$precision_derivatives <- pr
    set.seed(1)
    fit <- mod$sample(data = d$stan_data, init = list(make_init_fun(REPS[1], d)()), seed = 1,
                      chains = 1, threads_per_chain = threads_per_chain,
                      iter_warmup = 0, iter_sampling = 3, adapt_engaged = FALSE,
                      step_size = 0.0005, max_treedepth = 3,
                      refresh = 1, show_messages = FALSE, show_exceptions = FALSE,
                      output_dir = FIT_DIR)
    n_eval <- sum(fit$sampler_diagnostics(format = "df")$n_leapfrog__)
    sec    <- fit$time()$chains$sampling / n_eval
    it     <- ITER_WARMUP + ITER_SAMPLING
    data.frame(precision = pr, persons = d$stan_data$I, threads_per_chain = threads_per_chain,
               sec_per_eval = round(sec, 2),
               hours_per_rep_31_steps  = round(sec * 31  * it / 3600, 1),
               hours_per_rep_63_steps  = round(sec * 63  * it / 3600, 1),
               hours_per_rep_127_steps = round(sec * 127 * it / 3600, 1))
  }))
  cat("\n TIMING TEST \n")
  print(res, row.names = FALSE)
}



pop_pars <- c("beta_nu_intercept", "beta_nu_c", "beta_alpha", "alpha_pop_median",
              "beta_tau", "beta_w", "w_pop_median", "beta_sv", "sv_pop_median",
              "beta_sw", "sw_pop_mean", "beta_st0", "st0_pop_median",
              "sigma_nu_intercept", "sigma_nu_c", "sigma_alpha_intercept",
              "sigma_tau_intercept", "sigma_w", "sigma_sv", "sigma_sw", "sigma_st0")
subj_pars <- c("nu_intercept", "nu_c", "a", "t0", "w", "sv", "sw", "st0", "mean_ndt")

fit_one_rep <- function(r) {
  rep_tag <- sprintf("rep%02d", r)
  out_dir <- file.path(FIT_DIR, rep_tag)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(out_dir, "pop_summary.rds"))) {
    cat(rep_tag, "already fitted --> skipped\n"); return(invisible(NULL))
  }
  d <- make_data(r, N_PERSON)
  set.seed(9000 + r)
  init_fun <- make_init_fun(r, d)

  t_start <- Sys.time()
  fit <- mod$sample(
    data = d$stan_data, init = init_fun, seed = 9000 + r,
    chains = n_chains, parallel_chains = n_chains, threads_per_chain = threads_per_chain,
    iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
    adapt_delta = ADAPT_DELTA, max_treedepth = 8, step_size = 0.1,
    refresh = 20, show_exceptions = FALSE,
    output_dir = out_dir, output_basename = sprintf("ddm7hier_%s", rep_tag)
  )
  runtime_h <- as.numeric(difftime(Sys.time(), t_start, units = "hours"))

  pop_summ <- fit$summary(pop_pars, "mean", "median", "sd", "quantile2",
                          "rhat", "ess_bulk", "ess_tail") %>% mutate(rep = r)
  subj_summ <- fit$summary(subj_pars, "mean", "median", "sd", "quantile2",
                           "rhat", "ess_bulk", "ess_tail") %>%
    mutate(param = sub("\\[.*", "", variable),
           pid   = as.integer(gsub(".*\\[|\\]", "", variable)),
           subj  = d$subj_ids[pid], rep = r)

  diag <- fit$diagnostic_summary(quiet = TRUE)
  diag_df <- data.frame(rep = r, model = "ddm7_hier",
                        divergences   = sum(diag$num_divergent),
                        max_treedepth = sum(diag$num_max_treedepth),
                        min_ebfmi     = min(diag$ebfmi),
                        max_rhat      = max(c(pop_summ$rhat, subj_summ$rhat), na.rm = TRUE),
                        min_ess_bulk  = min(c(pop_summ$ess_bulk, subj_summ$ess_bulk), na.rm = TRUE),
                        runtime_h = runtime_h, n_person = d$stan_data$I, precision = PRECISION)

  saveRDS(pop_summ,   file.path(out_dir, "pop_summary.rds"))
  saveRDS(subj_summ,  file.path(out_dir, "subj_summary.rds"))
  saveRDS(diag_df,    file.path(out_dir, "diagnostics.rds"))
  saveRDS(d$subj_ids, file.path(out_dir, "subj_ids.rds"))

  cat(sprintf("%s done | %.1f h | div=%d | treedepth=%d | max Rhat=%.3f | min ESS=%.0f\n",
              rep_tag, runtime_h, diag_df$divergences, diag_df$max_treedepth,
              diag_df$max_rhat, diag_df$min_ess_bulk))
  print(pop_summ %>% select(variable, mean, q5, q95, rhat), n = Inf)
  invisible(diag_df)
}

if (!TEST_MODE) {
  for (r in REPS) {
    tryCatch(fit_one_rep(r),
             error = function(e) cat(sprintf("rep%02d FAILED: %s\n", r, conditionMessage(e))))
  }
  diag_all <- bind_rows(lapply(REPS, function(r) {
    f <- file.path(FIT_DIR, sprintf("rep%02d", r), "diagnostics.rds")
    if (file.exists(f)) readRDS(f)
  }))
  print(diag_all)
}
