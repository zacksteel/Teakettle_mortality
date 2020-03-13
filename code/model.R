## Purpose: Teakettle Mortality Final Model

## Bring in data
d <- read.csv("data/modeldata_stand.csv")

## Set up submodels
## density model
dens_mod <- bf(mvbind(n_size1_s, n_lg_s,
                      nba_eng_hs, nba_eng_ns,
                      nba_rtb_hs, nba_rtb_ns,
                      nba_mpb_hs, nba_mpb_ns,
                      nba_jpb_hs, nba_jpb_ns) ~
                 0 + treat + (1|p|plot),
               family = gaussian("identity"))

## Vigor model
vigor_mod <- bf(growth_z ~ 
                  (1 + n_size1_s + n_size1_s:dbh_s +
                     n_lg_s + n_lg_s:dbh_s | species)  +
                  (1|p|plot),
                family = gaussian("identity"))

## Infestation model

## Note the use of subset to ignore irrelevant species

eng_mod <- bf(f_eng | subset(abma_abco) ~ 
                (1 + burned + growth_z + dbh_s + nba_eng_hs + nba_eng_ns):abma  + 
                (1 + burned + growth_z + dbh_s + nba_eng_hs + nba_eng_ns):abco  + 
                (1|p|plot),
              family = bernoulli("logit"))

mpb_mod <- bf(p_mpb | subset(pila) ~ 
                (1 + burned + growth_z + dbh_s + nba_mpb_hs + nba_mpb_ns):pila + 
                (1|p|plot),
              family = bernoulli("logit"))

jpb_mod <- bf(p_jpb | subset(pije) ~ 
                (1 + burned + growth_z + dbh_s + nba_jpb_hs + nba_jpb_ns):pije + 
                (1|p|plot),
              family = bernoulli("logit"))

rtb_mod <- bf(p_rtb | subset(pije_pila) ~ 
                (1 + burned + growth_z + dbh_s + nba_rtb_hs + nba_rtb_ns):pije + 
                (1 + burned + growth_z + dbh_s + nba_rtb_hs + nba_rtb_ns):pila + 
                (1|p|plot),
              family = bernoulli("logit"))

## Mortality model
mort_mod <- bf(dmort ~ 0 +
                 ## only interested in species-specific effects since beetles are specialists
                 f_eng:abma + f_eng:abco + 
                 p_rtb:pije + p_rtb:pila + 
                 p_mpb:pila + 
                 p_jpb:pije +
                 ## Interactions for non-additive effects of multiple beetles
                 p_rtb:p_mpb:pila +
                 p_rtb:p_jpb:pije +
                 (1 + growth_z + 
                    n_size1_s + n_size1_s:dbh_s +
                    n_lg_s + n_lg_s:dbh_s +
                    burned +
                    srad_s + twi_s | species) +
                 (1|p|plot),
               family = bernoulli("logit"))

## Run full model
priorb <- prior(normal(0,2), class = b)

mfull <- brm(dens_mod + #Change in density
               vigor_mod + #Growth/vigor
               eng_mod + #Beetle models
               mpb_mod +
               jpb_mod +
               rtb_mod +
               mort_mod, #Mortality model
             set_rescor(FALSE),
             prior = priorb,
             data = d, 
             control = list(adapt_delta = 0.95, max_treedepth = 15),
             chains = 3, 
             iter = 3000)

## Save data and model
save(mfull, d, file = "fit_models/full_model.RData")