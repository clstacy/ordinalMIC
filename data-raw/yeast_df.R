## code to prepare `yeast_df` dataset goes here
set.seed(2025)

strains  <- c("WT", "Mut")
treatments <- c("None", "Salt")
conc_seq <- seq(0, 4, length.out = 11)
reps     <- 3

yeast_df <- expand.grid(
  strain = factor(strains, levels = strains),
  treatment = factor(treatments, levels = treatments),
  conc   = conc_seq,
  rep    = seq_len(reps)
)

## ------------------------------------------------------------------
## 1.  Define α (intercept) and β (negative slope) for each combo
## ------------------------------------------------------------------
# Baselines
alpha_base <-  4.0      # WT / None
beta_base  <- -6.0      # WT / None  (negative ⇒ score ↓ as conc ↑)

# Additive adjustments
alpha_adj_strain     <- c(WT = 0.2,  Mut =  -0.3)   # Mut grows worse
alpha_adj_treatment  <- c(None = 0,  Salt = -0.2)   # Pretreatment worse start
beta_adj_strain      <- c(WT = 0.0,  Mut =  -0.8)   # WT slope a bit flatter
beta_adj_treatment   <- c(None = 0,  Salt = 2)   # Pretreated slope flatter
beta_intxn           <- -0.7                    # extra interaction effect

with(yeast_df, {
  alpha <- alpha_base +
    alpha_adj_strain[strain] +
    alpha_adj_treatment[treatment]

  beta  <- beta_base +
    beta_adj_strain[strain] +
    beta_adj_treatment[treatment] +
    beta_intxn * (strain == "Mut" & treatment == "Salt")

  linpred <- alpha + beta * log1p(conc)         # note: beta is already negative
  linpred
}) -> linpred

# negative slope (↓ as conc ↑)
# beta  <- ifelse(yeast_df$strain == "Mut", 1.1, 0.9)
# alpha <- ifelse(yeast_df$strain == "Mut", 3.0, 2.0)
# linpred <- alpha - beta * yeast_df$conc

cutpoints <- c(-3, -1, 1, 3)
latent    <- linpred + 0.5*stats::rlogis(nrow(yeast_df))

yeast_df$score <- ordered(
  findInterval(latent, c(-Inf, cutpoints, Inf)) - 1,
  levels = 0:4
)

# yeast_df$rep <- NULL


usethis::use_data(yeast_df, overwrite = TRUE)
