rm(list = ls())
library(tidyverse)
library(dplyr)
library(stargazer)
library(gridExtra)
library(panelView)
library(lfe)
library(interflex)
library(gsynth)
library(sf)
library(rnaturalearthdata)
library(rnaturalearth)

####################################################################################
# Replication code for Right, Springman, Wibbels - Pushing Back or Backing Down
# 
# Please direct questions to Lucy Right (lucille.right@yale.edu)
####################################################################################

# Set WD
wd = NA
if(Sys.info()["sysname"] == "Linux"){ wd = Sys.getenv("PWD") } else{ wd = Sys.getenv("USERPROFILE") }
setwd(paste0(wd, "/Users/nourkebbi/Desktop/Trinity ASDS/Semester 2/Applied Statistical Analysis II POP77004-202526/dataverse_files"))


## PREP DATA #######################################################################
# ----------------------------------------------------------------------------------
## Load data  ----------------------------------------------------------------------

df = read.csv("/Users/nourkebbi/Desktop/Trinity ASDS/Semester 2/Applied Statistical Analysis II POP77004-202526/dataverse_files/rsw_dataset.csv") # dyadic dataset

# ----------------------------------------------------------------------------------
## Transform Variables -------------------------------------------------------------

# -- Lag covariates ------------------------------------------------------------------

cov = names(df %>% select(r_pop.ln, r_gdppc.ln, d_gdppc.ln, tradeflow.ln, dis_death.ihs, 
                          human_rights, civil_war, cpj.ihs, oppseat, v2x_libdem, 
                          v2xel_frefair, v2x_freexp_altinf, v2meharjrn, v2csreprss))

df = df %>% 
  group_by(dyad_id) %>%
  mutate(across(all_of(cov), ~dplyr::lag(.x), .names = "{col}.l")) %>%
  ungroup()


# -- Identify small countries -----
small_countries = (df %>% 
                     select(RecipientName, r_iso, year, r_pop) %>% 
                     unique() %>%
                     group_by(RecipientName, r_iso) %>%
                     # Indicator for if pop < 2mil every year in the dataset
                     summarize(small_all_years = ifelse(sum(r_pop < 2000000 | is.na(r_pop)) == n_distinct(year), 1, 0)) %>%
                     # Retain only small countries (countries < 2mil every year)
                     filter(small_all_years == 1))$RecipientName

df = df %>% mutate(r_small = ifelse(RecipientName %in% small_countries, 1, 0))

# Print small countries for appendix
# name = "Small Countries: " ; for(i in 1:length(small_countries)){name = paste0(name, small_countries[i], ", ")}; print(name)

rm(small_countries)
# ----------------------------------------------------------------------------------
## Create measures of democracy orientation ----------------------------------------

# -- DOR: Varies within donor, constant within dyad and over time --------

# DOR measures relationship orientation within dyad for *treated dyads* only (treatment provides reference year)
# 4 countries pass first law before 2007 - can't get relationship orientation for Afghanistan, Eritrea, Peru, Sudan
# Measured in the 3 years prior to law *onset* (if passage == 2008, treatment "turns on" in 2009, DOR is based off of 2008 + 2007 + 2006)

rel_orient = df %>% 
  filter(treated_dyad == 1 & first.dir.sig > 2006) %>%  # Keep only dyad years containing recipients countries that pass a law later than 2006 
  select(dyad_id, d_iso, r_iso, year, first.dir.sig, commit_tot, commit_dem, commit_dem_v2) %>%
  filter(year <= first.dir.sig & year >= first.dir.sig - 2) %>% 
  mutate(ref_year = paste0("t", first.dir.sig - year)) %>% # reference year index: year of, year before, and two years before law passage
  select(-year) %>%
  pivot_wider(names_from = ref_year, values_from = c(commit_tot, commit_dem, commit_dem_v2)) %>% # spread to wide using the ref_year
  mutate(TOTAL = commit_tot_t2 + commit_tot_t1 + commit_tot_t0,
         DEM_TOT = commit_dem_t2 + commit_dem_t1 + commit_dem_t0,
         DEM_TOT_v2 = commit_dem_v2_t2 + commit_dem_v2_t1 + commit_dem_v2_t0, 
         DEM_PCT = ifelse(TOTAL == 0, NA, DEM_TOT/TOTAL), # our coding
         DEM_PCT_v2 = ifelse(TOTAL == 0, NA, DEM_TOT_v2/TOTAL)) # robustness check coding

# define advocacy oriented relationship as >= 75th percentile among all donors to single country
rel_orient = rel_orient %>%
  group_by(r_iso) %>%
  mutate(quartile.75 = quantile(DEM_PCT, na.rm = T)[4],
         quartile.75_v2 = quantile(DEM_PCT_v2, na.rm = T)[4]) %>% 
  ungroup() %>%
  mutate(DOR = ifelse(DEM_PCT >= quartile.75, 1, 0),
         DOR_v2 = ifelse(DEM_PCT_v2 >= quartile.75_v2, 1, 0)) %>%   
  select(dyad_id, DEM_PCT:DOR_v2)


# -- DOD: Constant within donor, varies over time and within dyad ---------

# DOD measures donor's overall orientation by year, taking into account aid to *all* recipients in the three years prior
# Donor's Democracy Orientation in year t is based off spending in T-1 + T-2 + T-3 
# if passage == 2008, indicator turns on in 2009, we want AOD in 2009 measured as 2008 + 2007 + 2006 (equivalent temporal range as DOR)

dnr_orient <- setNames(data.frame(matrix(nrow = 0, ncol = 11)), 
                       c("d_iso","DEM_DNR", "DEM_DNR_v2","TOT","year","DEM_PCT_DNR", "DEM_PCT_DNR_v2", 
                         "quartile.75.DNR", "quartile.75.DNR_v2", "DOD", "DOD_v2"))

for(i in 2008:2019){
  temp <- df %>%
    filter(year < i & year >= (i - 3) 
           & zero_amt == 0) %>%       # exclude 0 commitment years - they will be NA
    group_by(d_iso) %>%
    summarize(DEM_DNR = sum(commit_dem, na.rm = T),
              DEM_DNR_v2 = sum(commit_dem_v2, na.rm = T),
              TOT = sum(commit_tot, na.rm = T)) %>%
    mutate(year = i,
           DEM_PCT_DNR = DEM_DNR/TOT,
           DEM_PCT_DNR_v2 = DEM_DNR_v2/TOT) 
  
  quart <- temp %>%
    summarize(quartile.75.DNR = quantile(DEM_PCT_DNR, na.rm = T)[4],
              quartile.75.DNR_v2 = quantile(DEM_PCT_DNR_v2, na.rm = T)[4],
              year = unique(year))
  
  temp <- left_join(temp, quart) %>%
    mutate(# our coding
           DOD = case_when(DEM_PCT_DNR >= quartile.75.DNR & quartile.75.DNR > 0 ~ 1,
                           quartile.75.DNR == 0 & DEM_PCT_DNR > 0 ~ 1, 
                           DEM_PCT_DNR < quartile.75.DNR & quartile.75.DNR > 0 ~ 0, 
                           TRUE ~ NA),
           # robustness check coding
           DOD_v2 = case_when(DEM_PCT_DNR_v2 >= quartile.75.DNR_v2 & quartile.75.DNR_v2 > 0 ~ 1,
                           quartile.75.DNR_v2 == 0 & DEM_PCT_DNR_v2 > 0 ~ 1, 
                           DEM_PCT_DNR_v2 < quartile.75.DNR_v2 & quartile.75.DNR_v2 > 0 ~ 0, 
                           TRUE ~ NA))
  
  dnr_orient <- rbind(dnr_orient, temp)
}


rm(temp, quart)

# -- Join with main datasets -------

df = left_join(df, rel_orient)
df = left_join(df, dnr_orient %>% select(d_iso, year, DEM_PCT_DNR:DOD_v2))

# Create "event history" dataset without control units
event_hist = df %>% filter(treated_dyad == 1)

# ----------------------------------------------------------------------------------
## Create descriptive figures and tables -------------------------------------------

# -- Figure 1 ----------
map = df %>% 
  group_by(r_iso) %>% 
  summarize(treated_dyad = max(treated_dyad),
            max = max(year)) %>% 
  filter(max == 2019) %>% # remove countries that transition out of recipient status
  select(r_iso, treated_dyad)

world = ne_countries(scale = "medium", 
                      continent = c("africa", "north america", "south america", 
                                    "europe", "asia", "Oceania"), 
                      returnclass = "sf")

# Merge datasets by ISO country code
map_data <- world %>%
  left_join(map, by = c("iso_a3" = "r_iso")) 

jpeg("./fig_1_map_of_laws.jpg", width = 8, height = 5, units = "in", res = 300)
ggplot() +
  geom_sf(data = map_data,
          aes(fill = as.factor(treated_dyad)),
          color = "darkgrey")  +
  scale_fill_manual(labels = c("No NGO Law", "NGO Law", "Non-Recipient"),
                      name = "",
                      values = c("white", "#f56f42"))+
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "bottom",            
    legend.direction = "horizontal",
    legend.title.align = 0,   
    legend.key.width = unit(1.5, "cm") 
  )
dev.off()
  

# -- Figures 2, 3 ------ 

# Figure 2 ---
do = df %>% 
  group_by(d_iso, DonorName) %>% 
  summarize(share_dem = sum(commit_dem, na.rm = T)/sum(commit_tot, na.rm =T))

jpeg("./fig_2_democracy aid.jpeg", width = 8, height = 5, units = "in", res = 300)
ggplot(do, aes(y = reorder(DonorName, share_dem), x = share_dem)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  theme_bw()+
  labs(title = "Democracy Aid as Share of Total Aid, 2005 - 2019",
       subtitle = "Among Large OECD Donors",
       x = "Democracy Aid as % of Total Aid",
       y = "Donor")
dev.off()

# Figure 3 ---

df_r = df %>% 
  filter(r_small ==0) %>%
  group_by(RecipientName, treated_dyad, r_iso, year, first.dir.sig) %>%
  summarize(commit = sum(commit_tot, na.rm = T)) %>%
  mutate(treatment = ifelse(treated_dyad == 0, -1,
                            ifelse(first.dir.sig == year, 1,
                                   ifelse(year < first.dir.sig, 0,
                                          ifelse(year > first.dir.sig, 2, NA))))) 

# treated
jpeg("./fig_3_treated.jpg", width = 6, height = 11, units = "in", res = 300)
panelview(commit ~ treatment, 
          data = df_r %>% filter(treated_dyad == 1) %>% arrange(first.dir.sig), 
          index = c("r_iso", "year"), 
          ylab = "Recipient Country", xlab = "Year", by.timing = T,
          main = "Aid Recipients with NGO Law",
          cex.main = 14, cex.legend = 10, background = "white",
          legend.labs = c("Pre-Law", "Law Passage", "Post-Law"),
          axis.lab.gap = c(0,0), color = c("lightblue", "darkorange", "darkblue"))
dev.off()

# control
jpeg("./fig_3_control.jpg", width = 6, height = 11, units = "in", res = 300)
panelview(commit ~ treatment, 
          data = df_r %>% filter(treated_dyad == 0) %>% arrange(first.dir.sig), 
          index = c("r_iso", "year"), 
          ylab = "Recipient Country", xlab = "Year",  
          main = "Aid Recipients with No NGO Law",
          cex.main = 14, cex.legend = 10, background = "white",
          legend.labs = c("No NGO Law", "Not Aid-Receiving"),
          axis.lab.gap = c(0,0), color = c("Grey", "White"))
dev.off()

rm(df_r)

# -- Tables A3, A4, A5 -----

# Table A3: Relationship Orientation ---
DOR = df %>%
  filter(treated_dyad == 1 & first.dir.sig == year & r_small == 0) %>%
  group_by(RecipientName, year) %>%
  summarize(quartile.75 = round(max(quartile.75, na.rm = T), 3)) %>%
  mutate(quartile.75 = ifelse(quartile.75 == "-Inf", NA, quartile.75))

DOR = left_join(DOR, df %>%
  filter(DOR == 1 & year == first.dir.sig & r_small == 0) %>%
  group_by(RecipientName, first.dir.sig) %>%
  summarize(DOR = paste(d_iso, collapse = ", "))) %>%
  select(RecipientName, year, quartile.75, DOR) %>%
  dplyr::rename(`Recipient Country` = RecipientName,
         `Passsage Year` = year,
         `75th percentile` = quartile.75,
         `Donors in DO relationship` = DOR)

stargazer(DOR, summary = F, rownames = F, type = "latex", 
          title = "Democracy-Oriented Relationships in Countries that Pass NGO Laws",
          digits = 3, column.sep.width = "3pt", font.size = "scriptsize", label = "app:tab_relorient", 
          out = "./tab_A3_relationship_orientation.tex")

rm(DOR)

# Table A4: Donor Orientation ---
DOD <- dnr_orient %>%
  filter(DOD == 1) %>%
  group_by(year) %>%
  summarize(`75th percentile` = round(unique(quartile.75.DNR), 3),
            `Democracy-Oriented Donors` = paste(d_iso, collapse = ", ")) %>%
  dplyr::rename(Year = year)

stargazer(DOD, summary = F, rownames = F, type = "latex",
          title = "Democracy-Oriented Donors, 2008-2019", label = "app:tab_dnrorient",
          out = "./tab_A4_donor_orientation.tex")

# Table A5: Mean Annual Commitments ---

descriptives <- function(x) {
  cbind.data.frame(
    'n' = sum(is.na(x)==F),
    'Mean' = round(mean(x, na.rm=T), 0))
}


sum <- df %>% filter(r_small == 0) %>% select(commit_tot, commit_dem, commit_econ)
sum_do <- df %>% filter(DOD == 1 & r_small == 0) %>% select(commit_tot, commit_dem, commit_econ)
sum_ndo <- df %>% filter(DOD == 0 & r_small == 0) %>% select(commit_tot, commit_dem, commit_econ)

# Extracting descriptive statistics 
sum_table <- cbind(c(do.call(rbind, lapply(sum, descriptives)) %>% rownames_to_column(var = "Variable"),
                     do.call(rbind, lapply(sum_do, descriptives))),
                   do.call(rbind, lapply(sum_ndo, descriptives)))

stargazer(sum_table, type = "text", summary = F, rownames = F, label = "app:tab_summary",
          out = "./tab_A5_mean_annual_aid.tex", title = "Average Annual Aid Flows, by Donor-Orientation")


## TWFE MODELS ######################################################################
# ----------------------------------------------------------------------------------
## Main Results --------------------------------------------------------------------

# -- Overall Response --------------------------------------------------------------
mod.0 <- felm(tot.log ~ dir.sig.event.l  | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))
mod.1 <- felm(dem.log ~ dir.sig.event.l  | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))
mod.2 <- felm(econ.log ~ dir.sig.event.l  | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))

models <- list(mod.0, mod.1, mod.2)

number_fe <- NA

for(i in 1:length(models)){
  number_fe[i] <- length(levels((models)[[i]]$fe$dyad_id[[1]]))
}

stargazer(models, 
          title = "Homogenous Response to NGO Law Enactment",
          label = "tab:event_response",
          dep.var.labels = c("Total Aid", "Democracy Aid", "Economic Aid"),
          covariate.labels = c("NGO Law Passage"),
          add.lines = list(c("Year FEs", rep("\\checkmark", length(models))), 
                           c("Dyad FEs", rep("\\checkmark", length(models))),
                           c("Donor FEs", rep("$\\times$", length(models))),
                           c("Number of Dyads", number_fe)),
          notes.append = TRUE,
          #notes.align = "m{45em}",
          font.size = "footnotesize",
          omit.stat = "ser",
          type = "text",
          out = "./tab_1_main_overall.tex")

# -- By Democracy Orientation --------------------

# RELATIONSHIP ORIENTATION
rel.0 <- felm(tot.log ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
rel.1 <- felm(dem.log ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))
rel.2 <- felm(econ.log ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))

# DONOR ORIENTATION
dnr.0 <- felm(tot.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
dnr.1 <- felm(dem.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))
dnr.2 <- felm(econ.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))

models_r <- list(rel.0, rel.1, rel.2)
models_d <- list(dnr.0, dnr.1, dnr.2)

number_fe_r <- NA
number_fe_d <- NA

for(i in 1:length(models_r)){
  number_fe_r[i] <- length(levels((models_r)[[i]]$fe$d_iso[[1]]))
  number_fe_d[i] <- length(levels((models_d)[[i]]$fe$dyad_id[[1]]))
}

stargazer(list(models_r, models_d), 
          title = "Marginal Effect of Law Passage by Democracy Orientation",
          label = "tab:me",
          dep.var.labels = c("Total Aid", "Democracy Aid", "Economic Aid", "Total Aid", "Democracy Aid", "Economic Aid"),
          covariate.labels = c("NGO Law Passage", "DO Relationship (DOR)", "NGO Law * DOR", "DO Donor (DOD)", "NGO Law * DOD"),
          add.lines = list(c("Year FEs", rep("\\checkmark", 6)), 
                           c("Dyad FEs", "$\\times$", "$\\times$", "$\\times$", "\\checkmark", "\\checkmark", "\\checkmark"),
                           c("Donor FEs", "\\checkmark", "\\checkmark", "\\checkmark", "$\\times$", "$\\times$", "$\\times$"),
                           c("No. of Dyad FE", "NA", "NA", "NA", number_fe_d),
                           c("No. of Donor FE", number_fe_r, "NA", "NA", "NA")),
          font.size = "footnotesize",
          omit.stat = "ser",
          type = "latex", 
          no.space = T,
          out = "./tab_2_main_interaction.tex")


# ----------------------------------------------------------------------------------
## Robustness Checks ---------------------------------------------------------------

# -- Include Control Dyads (staggered DiD) ---------------------------------------

did.0 <- felm(tot.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, df %>% filter(r_small == 0)) 
did.1 <- felm(dem.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, df %>% filter(r_small == 0))
did.2 <- felm(econ.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, df %>% filter(r_small == 0))

models <- list(did.0, did.1, did.2)

number_fe <- NA

for(i in 1:length(models)){
  number_fe[i] <- length(levels((models)[[i]]$fe$dyad_id[[1]]))
}

stargazer(models, 
          title = "Staggered DID: Inclusion of Never-Treated Dyads",
          label = "tab:controls_dnr_orient",
          dep.var.labels = c("Total Aid", "Democracy Aid", "Economic Aid"),
          covariate.labels = c("NGO Law Passage", "DOD", "NGO Law*DOD"),
          add.lines = list(c("Year FEs", rep("\\checkmark", length(models))), 
                           c("Dyad FEs", rep("\\checkmark", length(models))),
                           c("Donor FEs", rep("$\\times$", length(models))),
                           c("No. of Dyad FEs", number_fe),
                           c("No. of Donor FEs", rep("NA", length(models)))),
          font.size = "footnotesize",
          omit.stat = "ser",
          type = "text",
          out = "./tab_A6_staggered_did.tex")

# -- Add covariates ----------------------------------------------------------------

# Relationship Orientation
rel.0.cov <- felm(tot.log ~ dir.sig.event.l*DOR 
                  + tradeflow.ln.l + r_pop.ln.l + r_gdppc.ln.l + v2x_libdem.l + human_rights.l + v2csreprss.l + civil_war.l + dis_death.ihs.l
                  | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
rel.1.cov <- felm(dem.log ~ dir.sig.event.l*DOR 
                  + tradeflow.ln.l + r_pop.ln.l + r_gdppc.ln.l + v2x_libdem.l + human_rights.l + v2csreprss.l + civil_war.l + dis_death.ihs.l
                  | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))
rel.2.cov <- felm(econ.log ~ dir.sig.event.l*DOR 
                  + tradeflow.ln.l + r_pop.ln.l + r_gdppc.ln.l + v2x_libdem.l + human_rights.l + v2csreprss.l + civil_war.l + dis_death.ihs.l
                  | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))

# Donor Orientation
dnr.0.cov <- felm(tot.log ~ dir.sig.event.l*DOD 
                  + tradeflow.ln.l + r_pop.ln.l + r_gdppc.ln.l + v2x_libdem.l + human_rights.l + v2csreprss.l + civil_war.l + dis_death.ihs.l
                  | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
dnr.1.cov <- felm(dem.log ~ dir.sig.event.l*DOD 
                  + tradeflow.ln.l + r_pop.ln.l + r_gdppc.ln.l + v2x_libdem.l + human_rights.l + v2csreprss.l + civil_war.l + dis_death.ihs.l
                  | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))
dnr.2.cov <- felm(econ.log ~ dir.sig.event.l*DOD 
                  + tradeflow.ln.l + r_pop.ln.l + r_gdppc.ln.l + v2x_libdem.l + human_rights.l + v2csreprss.l + civil_war.l + dis_death.ihs.l
                  | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))

models_r <- list(rel.0.cov, rel.1.cov, rel.2.cov)
models_d <- list(dnr.0.cov, dnr.1.cov, dnr.2.cov)

number_fe_r <- NA
number_fe_d <- NA

for(i in 1:length(models_r)){
  number_fe_r[i] <- length(levels((models_r)[[i]]$fe$d_iso[[1]]))
  number_fe_d[i] <- length(levels((models_d)[[i]]$fe$dyad_id[[1]]))
}

stargazer(list(models_r, models_d), 
          title = "Inclusion of Covariates",
          label = "tab:me",
          dep.var.labels = c("Total Aid", "Democracy Aid", "Economic Aid", "Total Aid", "Democracy Aid", "Economic Aid"),
          covariate.labels = c("NGO Law Passage", "DO Relationship (DOR)", "DO Donor (DOD)", 
                               "Dyadic Trade", "Recipient Population", "Recipient GDP per capita", 
                               "Democracy", "Human Rights", "CSO Repression", "Civil War", "Disaster Deaths",
                               "NGO Law * DOR", "NGO Law * DOD"),
          add.lines = list(c("Year FEs", rep("\\checkmark", 6)), 
                           c("Dyad FEs", "$\\times$", "$\\times$", "$\\times$", "\\checkmark", "\\checkmark", "\\checkmark"),
                           c("Donor FEs", "\\checkmark", "\\checkmark", "\\checkmark", "$\\times$", "$\\times$", "$\\times$"),
                           c("No. of Dyad FE", "NA", "NA", "NA", number_fe_d),
                           c("No. of Donor FE", number_fe_r, "NA", "NA", "NA")),
          font.size = "footnotesize",
          omit.stat = "ser",
          type = "text",
          no.space = T,
          out = "./tab_A7_covariates.tex")





# -- Alternative coding scheme -----------------------------------------------------

# Relationship Orientation
rel.0.code <- felm(tot.log ~ dir.sig.event.l*DOR_v2 | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
rel.1.code <- felm(dem.log.v2 ~ dir.sig.event.l*DOR_v2 | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))
rel.2.code <- felm(econ.log.v2 ~ dir.sig.event.l*DOR_v2 | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))

# Donor Orientation
dnr.0.code <- felm(tot.log ~ dir.sig.event.l*DOD_v2 | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
dnr.1.code <- felm(dem.log.v2 ~ dir.sig.event.l*DOD_v2 | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))
dnr.2.code <- felm(econ.log.v2 ~ dir.sig.event.l*DOD_v2 | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))

models_r <- list(rel.0.code, rel.1.code, rel.2.code)
models_d <- list(dnr.0.code, dnr.1.code, dnr.2.code)

number_fe_r <- NA
number_fe_d <- NA

for(i in 1:length(models_r)){
  number_fe_r[i] <- length(levels((models_r)[[i]]$fe$d_iso[[1]]))
  number_fe_d[i] <- length(levels((models_d)[[i]]$fe$dyad_id[[1]]))
}

stargazer(list(models_r, models_d), 
          title = "Alternative Coding of Democracy Aid",
          label = "tab:alt_coding",
          dep.var.labels = c("Total Aid", "Democracy Aid", "Economic Aid", "Total Aid", "Democracy Aid", "Economic Aid"),
          covariate.labels = c("NGO Law Passage", "DO Relationship (DOR)", "NGO Law * DOR", "DO Donor (DOD)", "NGO Law * DOD"),
          add.lines = list(c("Year FEs", rep("\\checkmark", 6)), 
                           c("Dyad FEs", "$\\times$", "$\\times$", "$\\times$", "\\checkmark", "\\checkmark", "\\checkmark"),
                           c("Donor FEs", "\\checkmark", "\\checkmark", "\\checkmark", "$\\times$", "$\\times$", "$\\times$"),
                           c("No. of Dyad FE", "NA", "NA", "NA", number_fe_d),
                           c("No. of Donor FE", number_fe_r, "NA", "NA", "NA")),
          font.size = "footnotesize",
          omit.stat = "ser",
          type = "text",
          no.space = T,
          out = "./tab_A8_alternative_coding.tex")



# -- Continuous measure of democracy orientation -----------------------------------

event_hist = event_hist %>%
  mutate(DEM_PCT_100 = DEM_PCT * 100,
         DEM_PCT_DNR_100 = DEM_PCT_DNR * 100)

# Relationship Orientation
rel.0.cont <- felm(tot.log ~ dir.sig.event.l*DEM_PCT_100 | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
rel.1.cont <- felm(dem.log ~ dir.sig.event.l*DEM_PCT_100 | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))
rel.2.cont <- felm(econ.log ~ dir.sig.event.l*DEM_PCT_100 | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))

# Donor Orientation
dnr.0.cont <- felm(tot.log ~ dir.sig.event.l*DEM_PCT_DNR_100 | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
dnr.1.cont <- felm(dem.log ~ dir.sig.event.l*DEM_PCT_DNR_100 | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))
dnr.2.cont <- felm(econ.log ~ dir.sig.event.l*DEM_PCT_DNR_100 | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))


models_r <- list(rel.0.cont, rel.1.cont, rel.2.cont)
models_d <- list(dnr.0.cont, dnr.1.cont, dnr.2.cont)

number_fe_r <- NA
number_fe_d <- NA

for(i in 1:length(models_r)){
  number_fe_r[i] <- length(levels((models_r)[[i]]$fe$d_iso[[1]]))
  number_fe_d[i] <- length(levels((models_d)[[i]]$fe$dyad_id[[1]]))
}

stargazer(list(models_r, models_d), 
          title = "Continuous Measure of Democracy Orientation",
          label = "tab:cont",
          dep.var.labels = c("Total Aid", "Democracy Aid", "Economic Aid", "Total Aid", "Democracy Aid", "Economic Aid"),
          covariate.labels = c("NGO Law Passage", "DOR (continuous)", "NGO Law * DOR (continuous)", "DOD (continuous)", "NGO Law * DOD (continuous)"),
          add.lines = list(c("Year FEs", rep("\\checkmark", 6)), 
                           c("Dyad FEs", "$\\times$", "$\\times$", "$\\times$", "\\checkmark", "\\checkmark", "\\checkmark"),
                           c("Donor FEs", "\\checkmark", "\\checkmark", "\\checkmark", "$\\times$", "$\\times$", "$\\times$"),
                           c("No. of Dyad FE", "NA", "NA", "NA", number_fe_d),
                           c("No. of Donor FE", number_fe_r, "NA", "NA", "NA")),
          font.size = "footnotesize",
          omit.stat = "ser",
          type = "text",
          out = "./tab_A9_continuous.tex")

# interflex plots
d = as.data.frame(event_hist %>% filter(r_small == 0))
D = "dir.sig.event.l" 
Dlabel<-"Law Passage"
vcov<-"cluster"
cl<-"r_iso" 
Y="dem.log"
Ylabel<-"Democracy Aid"

# Relationship
X = "DEM_PCT_100" 
FE = c("year", "d_iso") #donor fixed effects
Xlabel<-"Relationship Democracy Orientation (%)"
rel.int = interflex("binning", Y=Y,D=D,X="DEM_PCT_100",data=d, FE=FE,
                      Ylabel=Ylabel, Xlabel=Xlabel, Dlabel=Dlabel, 
                      ylim = c(-3, 2.5), vcov.type = vcov, cl = cl,
                      nbins = 4, # correspond to 0-25th, 25th-50th, 50th - 75th, 75th+ percentile
                      bin.labs = F, theme.bw = T,
                      na.rm = T)
jpeg("./fig_A1a_interflex_rel.jpg", width = 6, height = 4, unit = "in", res = 300)
rel.int$figure
dev.off()

# Donor
X = "DEM_PCT_DNR_100" 
FE = c("year", "dyad_id") #donor fixed effects
Xlabel<-"Donor Democracy Orientation (%)"
dnr.int = interflex("binning", Y=Y,D=D,X=X,data=d, FE=FE,
                    Ylabel=Ylabel, Xlabel=Xlabel, Dlabel=Dlabel, 
                    ylim = c(-3, 2.5), vcov.type = vcov, cl = cl,
                    nbins = 4, # correspond to 0-25th, 25th-50th, 50th - 75th, 75th+ percentile
                    bin.labs = F, theme.bw = T,
                    na.rm = T)
jpeg("./fig_A1b_interflex_dnr.jpg", width = 6, height = 4, unit = "in", res = 300)
dnr.int$figure
dev.off()

rm(d, D, Dlabel, vcov, cl, Y, Ylabel, X, FE, Xlabel, dnr.int, rel.int)

# -- Drop United States --------------------------------------------------------------

# Relationship Orientation
rel.0.usa <- felm(tot.log ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0 & d_iso != "USA")) 
rel.1.usa <- felm(dem.log ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0 & d_iso != "USA"))
rel.2.usa <- felm(econ.log ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0 & d_iso != "USA"))

# Donor Orientation
dnr.0.usa <- felm(tot.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0 & d_iso != "USA")) 
dnr.1.usa <- felm(dem.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0 & d_iso != "USA"))
dnr.2.usa <- felm(econ.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0 & d_iso != "USA"))

models_r <- list(rel.0.usa, rel.1.usa, rel.2.usa)
models_d <- list(dnr.0.usa, dnr.1.usa, dnr.2.usa)

number_fe_r <- NA
number_fe_d <- NA

for(i in 1:length(models_r)){
  number_fe_r[i] <- length(levels((models_r)[[i]]$fe$d_iso[[1]]))
  number_fe_d[i] <- length(levels((models_d)[[i]]$fe$dyad_id[[1]]))
}

stargazer(list(models_r, models_d), 
          title = "Dropping Dyads with USA",
          label = "tab:usa",
          dep.var.labels = c("Total Aid", "Democracy Aid", "Economic Aid", "Total Aid", "Democracy Aid", "Economic Aid"),
          covariate.labels = c("NGO Law Passage", "DO Relationship (DOR)", "NGO Law * DOR", "DO Donor (DOD)", "NGO Law * DOD"),
          add.lines = list(c("Year FEs", rep("\\checkmark", length(models))), 
                           c("Dyad FEs", "$\\times$", "$\\times$", "$\\times$", "\\checkmark", "\\checkmark", "\\checkmark"),
                           c("Donor FEs", "\\checkmark", "\\checkmark", "\\checkmark", "$\\times$", "$\\times$", "$\\times$"),
                           c("No. of Dyad FE", "NA", "NA", "NA", number_fe_d),
                           c("No. of Donor FE", number_fe_r, "NA", "NA", "NA")),
          font.size = "footnotesize",
          omit.stat = "ser",
          type = "text",
          no.space = T,
          out = "./tab_A10_drop_USA.tex")

# -- Exclude Zero Values --------------------------------------------------------------

# Relationship Orientation
rel.0.nozero <- felm(tot.log ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0 & commit_tot > 0))
rel.1.nozero <- felm(dem.log ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0 & commit_dem > 0))
rel.2.nozero <- felm(econ.log ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0 & commit_econ > 0))

# Donor Orientation
dnr.0.nozero <- felm(tot.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0 & commit_tot > 0))
dnr.1.nozero <- felm(dem.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0 & commit_dem > 0))
dnr.2.nozero <- felm(econ.log ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0 & commit_econ > 0))

models <- list(rel.0.nozero, rel.1.nozero, rel.2.nozero, dnr.0.nozero, dnr.1.nozero, dnr.2.nozero)

number_fe_r <- length(levels((models)[[1]]$fe$d_iso[[1]]))
number_fe_d <- length(levels((models)[[4]]$fe$dyad_id[[1]]))

stargazer(models, 
          title = "Non-Zero Aid Flows Only",
          label = "tab:het_positive",
          dep.var.labels = c("Total Aid", "Democracy Aid", "Economic Aid", "Total Aid", "Democracy Aid", "Economic Aid"),
          covariate.labels = c("NGO Law Passage", "DO Relationship (DOR)", "NGO Law * DOR", "DO Donor (DOD)", "NGO Law * DOD"),
          add.lines = list(c("Year FEs", rep("\\checkmark", length(models))), 
                           c("Dyad FEs", "$\\times$", "$\\times$", "$\\times$", "\\checkmark", "\\checkmark", "\\checkmark"),
                           c("Donor FEs", "\\checkmark", "\\checkmark", "\\checkmark", "$\\times$", "$\\times$", "$\\times$"),
                           c("No. Dyad FEs","NA", "NA", "NA", number_fe_d, number_fe_d, number_fe_d),
                           c("No. Donor FEs", number_fe_r, number_fe_r, number_fe_r, "NA", "NA", "NA")),
          font.size = "footnotesize",
          omit.stat = "ser",
          no.space = T,
          type = "text",
          out = "./tab_A11_nonzero.tex")

# -- IHS-transformation of outcomes ------------------------------------------------

# Relationship Orientation
rel.0.ihs <- felm(tot.ihs ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
rel.1.ihs <- felm(dem.ihs ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))
rel.2.ihs <- felm(econ.ihs ~ dir.sig.event.l*DOR | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))

# Donor Orientation
dnr.0.ihs <- felm(tot.ihs ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
dnr.1.ihs <- felm(dem.ihs ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))
dnr.2.ihs <- felm(econ.ihs ~ dir.sig.event.l*DOD | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))

models_r <- list(rel.0.ihs, rel.1.ihs, rel.2.ihs)
models_d <- list(dnr.0.ihs, dnr.1.ihs, dnr.2.ihs)

number_fe_r <- NA
number_fe_d <- NA

for(i in 1:length(models_r)){
  number_fe_r[i] <- length(levels((models_r)[[i]]$fe$d_iso[[1]]))
  number_fe_d[i] <- length(levels((models_d)[[i]]$fe$dyad_id[[1]]))
}

stargazer(list(models_r, models_d), 
          title = "IHS-transformed Outcomes",
          label = "tab:ihs",
          dep.var.labels = c("Total Aid (ihs)", "Democracy Aid (ihs)", "Economic Aid (ihs)", "Total Aid (ihs)", "Democracy Aid (ihs)", "Economic Aid (ihs)"),
          covariate.labels = c("NGO Law Passage", "DO Relationship (DOR)", "NGO Law * DOR", "DO Donor (DOD)", "NGO Law * DOD"),
          add.lines = list(c("Year FEs", rep("\\checkmark", length(models))), 
                           c("Dyad FEs", "$\\times$", "$\\times$", "$\\times$", "\\checkmark", "\\checkmark", "\\checkmark"),
                           c("Donor FEs", "\\checkmark", "\\checkmark", "\\checkmark", "$\\times$", "$\\times$", "$\\times$"),
                           c("No. of Dyad FE", "NA", "NA", "NA", number_fe_d),
                           c("No. of Donor FE", number_fe_r, "NA", "NA", "NA")),
          font.size = "footnotesize",
          omit.stat = "ser",
          type = "text",
          no.space = T,
          out = "./tab_A12_ihs.tex")


## GENERALIZED SYNTHETIC CONTROL ###################################################
# ----------------------------------------------------------------------------------
## Run Analysis --------------------------------------------------------------------

# -- Prep data ---------------------------------------------------------------------
# Note: GSC takes a long time to run. If you want to load the final results, skip this
# section and jump to "Report Results" below

# partition data
dyads_t_dod = (df %>% filter(year == first.dir.sig & treated_dyad == 1 & DOD == 1))$dyad_name 
dyads_t_ndo = (df %>% filter(year == first.dir.sig & treated_dyad == 1 & DOD == 0))$dyad_name
dyads_c = unique((df %>% filter(treated_dyad == 0))$dyad_name)

# prep data
prep.df = df %>% 
  filter(first.dir.sig %in% 2009:2018 | is.na(first.dir.sig)) %>%  # only retain treated dyads with at least 5 pre-treatment (incl. passage yr) periods & 1 post-treatment period
  filter(r_small == 0) %>%                                         # remove small countries
  mutate(treat = ifelse(is.na(first.dir.sig), 0,                   # treatment = 0 in year of passage and all prior years (consistent with gsynth package)
                        ifelse(first.dir.sig < year, 1, 0)),
         status = ifelse(dyad_name %in% dyads_c, "control",
                         ifelse(dyad_name %in% dyads_t_dod, "treated dod",
                                ifelse(dyad_name %in% dyads_t_ndo, "treated non-dod", NA)))) %>%                          
  select(dyad_id, dyad_name, year, r_iso, treated_dyad, treat, status, DOD,   
         tot.log, dem.log, econ.log,
         r_gdppc.ln, tradeflow.ln, dis_death.ihs, v2x_libdem, human_rights, civil_war, v2csreprss) 

# dyads should have same number of observations within the panel (removes newly-formed countries, new donors, or recipients who became donors)

incomplete = (prep.df %>%                    # 114/1938 dyads have less than 15 years coverage
                group_by(dyad_name) %>% 
                summarize(n_years = n_distinct(year)) %>% 
                filter(n_years < 15))$dyad_name                           

prep.df = prep.df %>% filter(!dyad_name %in% incomplete) # filter out incomplete dyads

# dyad with missing values for all obs of a covariate will be removed by model - remove proactively

cov_complete = prep.df %>% group_by(dyad_name) %>% summarise_all(funs(sum(is.na(.))))
cov_complete = (cov_complete %>% filter(r_gdppc.ln == 15 | tradeflow.ln == 15 | dis_death.ihs == 15 | v2x_libdem == 15 | human_rights == 15 | civil_war == 15 | v2csreprss == 15))$dyad_name  
prep.df = prep.df %>% filter(!dyad_name %in% cov_complete)  # filter out dyads with NAs for every year

nrow(unique(prep.df %>% filter(status == "treated non-dod") %>% select(dyad_name)))
# -- Overall -----------------------------------------------------------------------
# Note that models take a while to run. 
# Total Aid
scm_tot =  gsynth(tot.log ~ treat + r_gdppc.ln + tradeflow.ln + dis_death.ihs + v2x_libdem + human_rights + civil_war + v2csreprss, 
                  na.rm = T, data = prep.df, index = c("dyad_name", "year"), min.T0 = 5,
                  se = T, inference = "parametric", r = c(0,3), CV = T,
                  force = "two-way", nboots = 1000, seed = 5195, EM = TRUE)

# Economic Aid
scm_econ =  gsynth(econ.log ~ treat + r_gdppc.ln + tradeflow.ln + dis_death.ihs + v2x_libdem + human_rights + civil_war + v2csreprss, 
                   na.rm = T, data = prep.df, index = c("dyad_name", "year"), min.T0 = 5,
                   se = T, inference = "parametric", r = c(0,3), CV = T,
                   force = "two-way", nboots = 1000, seed = 5195, EM = TRUE)

# Democracy Aid
scm_dem =  gsynth(dem.log ~ treat + r_gdppc.ln + tradeflow.ln + dis_death.ihs + v2x_libdem + human_rights + civil_war + v2csreprss, 
                  na.rm = T, data = prep.df, index = c("dyad_name", "year"), min.T0 = 5,
                  se = T, inference = "parametric", r = c(0,3), CV = T,
                  force = "two-way", nboots = 1000, seed = 5195, EM = TRUE)

# -- By Democracy Orientation ------------------------------------------------------

prep.df.dem = prep.df %>% filter(status %in% c("control", "treated dod"))
prep.df.non = prep.df %>% filter(status %in% c("control", "treated non-dod"))

# Total Aid
scm_tot_dod = gsynth(tot.log ~ treat + r_gdppc.ln + tradeflow.ln + dis_death.ihs + v2x_libdem + human_rights + civil_war + v2csreprss, 
                     na.rm = T, data = prep.df.dem, index = c("dyad_name", "year"), min.T0 = 5,
                     se = T, inference = "parametric", r = c(0,3), CV = T,
                     force = "two-way", nboots = 1000, seed = 5195, EM = TRUE)

scm_tot_non = gsynth(tot.log ~ treat + r_gdppc.ln + tradeflow.ln + dis_death.ihs + v2x_libdem + human_rights + civil_war + v2csreprss, 
                     na.rm = T, data = prep.df.non, index = c("dyad_name", "year"), min.T0 = 5,
                     se = T, inference = "parametric", r = c(0,3), CV = T,
                     force = "two-way", nboots = 1000, seed = 5195, EM = TRUE)

# Democracy Aid
scm_dem_dod = gsynth(dem.log ~ treat + r_gdppc.ln + tradeflow.ln + dis_death.ihs + v2x_libdem + human_rights + civil_war + v2csreprss, 
                     na.rm = T, data = prep.df.dem, index = c("dyad_name", "year"), min.T0 = 5,
                     se = T, inference = "parametric", r = c(0,3), CV = T,
                     force = "two-way", nboots = 1000, seed = 5195, EM = TRUE)

scm_dem_non = gsynth(dem.log ~ treat + r_gdppc.ln + tradeflow.ln + dis_death.ihs + v2x_libdem + human_rights + civil_war + v2csreprss, 
                     na.rm = T, data = prep.df.non, index = c("dyad_name", "year"), min.T0 = 5,
                     se = T, inference = "parametric", r = c(0,3), CV = T,
                     force = "two-way", nboots = 1000, seed = 5195, EM = TRUE)

# Economic Aid
scm_econ_dod = gsynth(econ.log ~ treat + r_gdppc.ln + tradeflow.ln + dis_death.ihs + v2x_libdem + human_rights + civil_war + v2csreprss, 
                     na.rm = T, data = prep.df.dem, index = c("dyad_name", "year"), min.T0 = 5,
                     se = T, inference = "parametric", r = c(0,3), CV = T,
                     force = "two-way", nboots = 1000, seed = 5195, EM = TRUE)

scm_econ_non = gsynth(econ.log ~ treat + r_gdppc.ln + tradeflow.ln + dis_death.ihs + v2x_libdem + human_rights + civil_war + v2csreprss, 
                     na.rm = T, data = prep.df.non, index = c("dyad_name", "year"), min.T0 = 5,
                     se = T, inference = "parametric", r = c(0,3), CV = T,
                     force = "two-way", nboots = 1000, seed = 5195, EM = TRUE)

# -- Save results ------------------------------------------------------------------
# GSC takes a long time to run. Save results for use next time.

save(scm_dem, scm_dem_dod, scm_dem_non, scm_econ, scm_econ_dod, scm_econ_non, scm_tot, scm_tot_dod, scm_tot_non, file = "./gsc_results.Rdata")

# ----------------------------------------------------------------------------------
## Report Results ------------------------------------------------------------------

# -- Prep --------------------------------------------------------------------------
# If loading results instead of running them, use command below 
load("./gsc_results.Rdata")


# Function to read gsc individual effects into data frame
read_ind = function(data_array){
  
  dim1_names <- dimnames(data_array)[[1]]
  dim2_names <- dimnames(data_array)[[2]]
  dim3_names <- dimnames(data_array)[[3]]
  
  # Convert the array to a data frame
  temp <- as.data.frame.table(data_array)
  names(temp) <- c("time", "metric", "dyad_name", "value")
  
  # Pivot
  temp_wide <- temp %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    mutate(d_iso = str_split_fixed(dyad_name, "-", 2)[,2])
  
  return(temp_wide)
}

# Reshape Data 
res_tot = as.data.frame(scm_tot$est.att) %>% rownames_to_column()
res_dem = as.data.frame(scm_dem$est.att) %>% rownames_to_column()
res_econ = as.data.frame(scm_econ$est.att) %>% rownames_to_column()

tot_dod = as.data.frame(scm_tot_dod$est.att) %>% rownames_to_column()
dem_dod = as.data.frame(scm_dem_dod$est.att) %>% rownames_to_column()
econ_dod = as.data.frame(scm_econ_dod$est.att) %>% rownames_to_column()

tot_non = as.data.frame(scm_tot_non$est.att) %>% rownames_to_column()
dem_non = as.data.frame(scm_dem_non$est.att) %>% rownames_to_column()
econ_non = as.data.frame(scm_econ_non$est.att) %>% rownames_to_column()


# -- Plot Effects: All Donors ------------------------------------------------------

jpeg("./fig_4a_tot_all.jpg", width = 6, height = 4, units = "in", res = 300)
ggplot(res_tot, aes(x = as.numeric(rowname), y = ATT)) +
  #geom_col(aes(y = (n.Treated / 1000), fill = "Treated Dyads"), fill = "darkblue", alpha = 0.05) + 
  geom_vline(aes(xintercept = 0), size = 3, color = "darkgrey", alpha = 0.5)+
  geom_hline(aes(yintercept = 0), size = 1, color = "darkgrey", alpha = 0.5)+
  geom_line(color = "darkblue", size = 1) + 
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2) +  
  labs(x = "Time Relative to Treatment", y = "ATT (Total Aid)") + # Scale ESS to fit on the same plot
  theme_minimal()+
  xlim(-4,6)+ylim(-5,2.2)+
  theme(panel.grid.minor.x = element_blank())
dev.off()

jpeg("./fig_4b_dem_all.jpg", width = 6, height = 4, units = "in", res = 300)
ggplot(res_dem, aes(x = as.numeric(rowname), y = ATT)) +
  #geom_col(aes(y = (n.Treated / 1000), fill = "Treated Dyads"), fill = "darkblue", alpha = 0.05) + 
  geom_vline(aes(xintercept = 0), size = 3, color = "darkgrey", alpha = 0.5)+
  geom_hline(aes(yintercept = 0), size = 1, color = "darkgrey", alpha = 0.5)+
  geom_line(color = "darkblue", size = 1) + 
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2) +  
  labs(x = "Time Relative to Treatment", y = "ATT (Democracy Aid)") + 
  theme_minimal() +
  xlim(-4,6)+ylim(-5,2.2)+
  theme(panel.grid.minor.x = element_blank())
dev.off()

jpeg("./fig_A2a_econ_all.jpg", width = 6, height = 4, units = "in", res = 300)
ggplot(res_econ, aes(x = as.numeric(rowname), y = ATT)) +
  #geom_col(aes(y = (n.Treated / 1000), fill = "Treated Dyads"), fill = "darkblue", alpha = 0.05) + 
  geom_vline(aes(xintercept = 0), size = 3, color = "darkgrey", alpha = 0.5)+
  geom_hline(aes(yintercept = 0), size = 1, color = "darkgrey", alpha = 0.5)+
  geom_line(color = "darkblue", size = 1) + 
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2) +  
  labs(x = "Time Relative to Treatment", y = "ATT (Economic Aid)") + 
  theme_minimal() +
  xlim(-4,6)+ylim(-5,2.2)+
  theme(panel.grid.minor.x = element_blank())
dev.off()

# -- Plot Effects: Democracy-Oriented Donors -----------------------------------------------

jpeg("./fig_NR_tot_dod.jpg", width = 6, height = 4, units = "in", res = 300)
ggplot(tot_dod, aes(x = as.numeric(rowname), y = ATT)) +
  #geom_col(aes(y = (n.Treated / 500), fill = "Treated Dyads"), fill = "darkblue", alpha = 0.05) + 
  geom_vline(aes(xintercept = 0), size = 3, color = "darkgrey", alpha = 0.5)+
  geom_hline(aes(yintercept = 0), size = 1, color = "darkgrey", alpha = 0.5)+
  geom_line(color = "darkblue", size = 1) + 
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2) + 
  labs(x = "Time Relative to Treatment", y = "ATT (Total Aid)") + 
  theme_minimal() +
  xlim(-4,6)+ylim(-5,2.2)+
  theme(panel.grid.minor.x = element_blank())
dev.off()

jpeg("./fig_5a_dem_dod.jpg", width = 6, height = 4, units = "in", res = 300)
ggplot(dem_dod, aes(x = as.numeric(rowname), y = ATT)) +
  #geom_col(aes(y = (n.Treated / 500), fill = "Treated Dyads"), fill = "darkblue", alpha = 0.05) + 
  geom_vline(aes(xintercept = 0), linewidth = 3, color = "darkgrey", alpha = 0.5)+
  geom_hline(aes(yintercept = 0), linewidth = 1, color = "darkgrey", alpha = 0.5)+
  geom_line(color = "darkblue", linewidth = 1) + 
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2) +  
  labs(x = "Time Relative to Treatment", y = "ATT (Democracy Aid)") + 
  theme_minimal() +
  xlim(-4,6)+ylim(-5,2.2)+
  theme(panel.grid.minor.x = element_blank())
dev.off()

jpeg("./fig_A2b_econ_dod.jpg", width = 6, height = 4, units = "in", res = 300)
ggplot(econ_dod, aes(x = as.numeric(rowname), y = ATT)) +
  #geom_col(aes(y = (n.Treated / 500), fill = "Treated Dyads"), fill = "darkblue", alpha = 0.05) + 
  geom_vline(aes(xintercept = 0), size = 3, color = "darkgrey", alpha = 0.5)+
  geom_hline(aes(yintercept = 0), size = 1, color = "darkgrey", alpha = 0.5)+
  geom_line(color = "darkblue", size = 1) + 
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2) +  
  labs(x = "Time Relative to Treatment", y = "ATT (Economic Aid)") + 
  theme_minimal() +
  xlim(-4,6)+ylim(-5,2.2)+
  theme(panel.grid.minor.x = element_blank())
dev.off()
# -- Plot Effects: Non Democracy-Oriented ------------------------------------------

jpeg("./fig_NR_tot_non.jpg", width = 6, height = 4, units = "in", res = 300)
ggplot(tot_non, aes(x = as.numeric(rowname), y = ATT)) +
  #geom_col(aes(y = (n.Treated / 500), fill = "Treated Dyads"), fill = "darkblue", alpha = 0.05) + 
  geom_vline(aes(xintercept = 0), size = 3, color = "darkgrey", alpha = 0.5)+
  geom_hline(aes(yintercept = 0), size = 1, color = "darkgrey", alpha = 0.5)+
  geom_line(color = "darkblue", size = 1) + 
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2) +  # Confidence interval shading
  labs(x = "Time Relative to Treatment", y = "ATT (Total Aid)") + 
  theme_minimal() +
  xlim(-4,6)+ylim(-5,2.2)+
  theme(panel.grid.minor.x = element_blank())
dev.off()

jpeg("./fig_5b_dem_non.jpg", width = 6, height = 4, units = "in", res = 300)
ggplot(dem_non, aes(x = as.numeric(rowname), y = ATT)) +
  #geom_col(aes(y = (n.Treated / 500), fill = "Treated Dyads"), fill = "darkblue", alpha = 0.05) + 
  geom_vline(aes(xintercept = 0), size = 3, color = "darkgrey", alpha = 0.5)+
  geom_hline(aes(yintercept = 0), size = 1, color = "darkgrey", alpha = 0.5)+
  geom_line(color = "darkblue", size = 1) + 
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2) +  # Confidence interval shading
  labs(x = "Time Relative to Treatment", y = "ATT (Democracy Aid)") + 
  theme_minimal() +
  xlim(-4,6)+ylim(-5,2.2)+
  theme(panel.grid.minor.x = element_blank())
dev.off()

jpeg("./fig_A2c_econ_non.jpg", width = 6, height = 4, units = "in", res = 300)
ggplot(econ_non, aes(x = as.numeric(rowname), y = ATT)) +
  #geom_col(aes(y = (n.Treated / 500), fill = "Treated Dyads"), fill = "darkblue", alpha = 0.05) + 
  geom_vline(aes(xintercept = 0), size = 3, color = "darkgrey", alpha = 0.5)+
  geom_hline(aes(yintercept = 0), size = 1, color = "darkgrey", alpha = 0.5)+
  geom_line(color = "darkblue", size = 1) + 
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2) +  # Confidence interval shading
  labs(x = "Time Relative to Treatment", y = "ATT (Economic Aid)") + 
  theme_minimal() +
  xlim(-4,6)+ylim(-5,2.2)+
  theme(panel.grid.minor.x = element_blank())
dev.off()
# -- Output Table ------------------------------------------------------------------

dem_dod_tab = dem_dod %>%  select(rowname, ATT, S.E., p.value, n.Treated) %>%
  mutate(ATT = round(ATT, 3), S.E. = round(S.E., 3), p.value = round(p.value, 3)) %>%
  dplyr::rename(Year = rowname, SE = S.E., p = p.value, "Treated Dyads" = n.Treated)

stargazer(dem_dod_tab, summary = F, rownames = F,
          type = "text",   title = "Democracy Aid: Democracy-Oriented Donors",
          out = "./tab_NR_gsc_dem_dod.tex")

dem_non_tab = dem_non %>%  select(rowname, ATT, S.E., p.value, n.Treated) %>%
  mutate(ATT = round(ATT, 3), S.E. = round(S.E., 3), p.value = round(p.value, 3)) %>%
  dplyr::rename(Year = rowname, SE = S.E., p = p.value, "Treated Dyads" = n.Treated)

stargazer(dem_non_tab, summary = F, rownames = F,
          type = "text",   title = "Democracy Aid: Non-Democracy-Oriented Donors",
          out = "./tab_NR_gsc_dem_non.tex")


## DISCUSSION ######################################################################
# -- Control for Democratic Backsliding --------------------------------------------

# RELATIONSHIP ORIENTATION
rel.0.dem <- felm(tot.log ~ dir.sig.event.l*DOR 
                  + v2xel_frefair.l + v2x_freexp_altinf.l + v2meharjrn.l + cpj.ihs.l + oppseat.l
                  | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))
rel.1.dem <- felm(dem.log ~ dir.sig.event.l*DOR 
                  + v2xel_frefair.l + v2x_freexp_altinf.l + v2meharjrn.l + cpj.ihs.l + oppseat.l
                  | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))
rel.2.dem <- felm(econ.log ~ dir.sig.event.l*DOR
                  + v2xel_frefair.l + v2x_freexp_altinf.l + v2meharjrn.l + cpj.ihs.l + oppseat.l
                  | year + d_iso | 0 | r_iso, event_hist %>% filter(r_small == 0))

# DONOR ORIENTATION
dnr.0.dem <- felm(tot.log ~ dir.sig.event.l*DOD 
              + v2xel_frefair.l + v2x_freexp_altinf.l + v2meharjrn.l + cpj.ihs.l + oppseat.l
              | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0)) 
dnr.1.dem <- felm(dem.log ~ dir.sig.event.l*DOD
              + v2xel_frefair.l + v2x_freexp_altinf.l + v2meharjrn.l + cpj.ihs.l + oppseat.l
              | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))
dnr.2.dem <- felm(econ.log ~ dir.sig.event.l*DOD 
              + v2xel_frefair.l + v2x_freexp_altinf.l + v2meharjrn.l + cpj.ihs.l + oppseat.l
              | year + dyad_id | 0 | r_iso, event_hist %>% filter(r_small == 0))

models_r <- list(rel.0.dem, rel.1.dem, rel.2.dem)
models_d <- list(dnr.0.dem, dnr.1.dem, dnr.2.dem)

number_fe_r <- NA
number_fe_d <- NA

for(i in 1:length(models_r)){
  number_fe_r[i] <- length(levels((models_r)[[i]]$fe$d_iso[[1]]))
  number_fe_d[i] <- length(levels((models_d)[[i]]$fe$dyad_id[[1]]))
}

stargazer(list(models_r, models_d), 
          title = "Controlling for Democratic Backsliding",
          label = "tab:dem_bck",
          dep.var.labels = c("Total Aid", "Democracy Aid", "Economic Aid", "Total Aid", "Democracy Aid", "Economic Aid"),
          covariate.labels = c("NGO Law Passage", "DO Relationship", "DO Donor",
                               "Clean Elections", "Freedom of Expression", "Media Harassment", "Attacks on Journalists", "Opposition Seat Share", 
                               "NGO Law * DOR", "NGO Law * DOD"),
          add.lines = list(c("Year FEs", rep("\\checkmark", 6)), 
                           c("Dyad FEs", "$\\times$", "$\\times$", "$\\times$", "\\checkmark", "\\checkmark", "\\checkmark"),
                           c("Donor FEs", "\\checkmark", "\\checkmark", "\\checkmark", "$\\times$", "$\\times$", "$\\times$"),
                           c("No. of Dyad FE", "NA", "NA", "NA", number_fe_d),
                           c("No. of Donor FE", number_fe_r, "NA", "NA", "NA")),
          font.size = "scriptsize",
          omit.stat = "ser",
          type = "text",
          no.space = T,
          out = "./tab_A13_backsliding.tex")

# -- Channel -----------------------------------------------------------------------

demch = df %>%
  group_by(year) %>%
  summarize(dem_total = sum(commit_dem, na.rm = T),
            across(c(dem_ngo,dem_r_public,dem_d_public, dem_o_public, dem_bypass, dem_none), ~sum(., na.rm = T)),
            across(c(dem_ngo,dem_r_public,dem_d_public, dem_o_public, dem_bypass, dem_none), ~ ./dem_total, .names = "p.{.col}")) 

# Figure A3 --
jpeg("./fig_A3_channel.jpg", units = "in", width = 8, height = 4, res = 300)
par(mar = c(5,5,2,9), xpd = F)
plot(x = demch$year, y = demch$p.dem_ngo, type = "l", col = "coral4", lwd = 2,
     #main = "Advocacy Aid Channels by year, 2005-2019",
     xlab = "Year",
     ylab = "Share of Democracy Aid through channel",
     ylim = c(0,.6), xlim = c(2005,2019))
lines(x = demch$year, y = demch$p.dem_bypass, type = "l", col = "brown1", lwd = 2)
lines(x = demch$year, y = demch$p.dem_r_public, type = "l", col = "blue4", lty = 2, lwd = 2)
lines(x = demch$year, y = demch$p.dem_d_public, type = "l", col = "cyan4", lty = 2, lwd = 2)
lines(x = demch$year, y = demch$p.dem_o_public, type = "l", col = "darkgrey", lty = 2,lwd = 2)
lines(x = demch$year, y = demch$p.dem_none, type = "l", col = "lightgrey", lty = 3, lwd = 2)
legend(x = 2020.25, y = .45,
       legend = c("NGOs", "Other Bypass", "Recipient Gov't", "Donor Gov't", "Unspecified Gov't", "Unclassified"),
       lty = c(1,1,2,2,2,3),
       lwd = c(2,2,2,2,2,2),
       col = c("coral4", "brown1", "blue4", "cyan4", "darkgrey", "lightgrey"),
       xpd = T,
       bty = "n",
       cex = .8)
dev.off()

# TWFE models ---

# Log values and filter out years with limited channel coverage
channel = event_hist %>% 
  mutate(across(c(dem_bypass, dem_ngo, dem_r_public, dem_d_public), ~log(.x + 1), .names = "{.col}.log")) %>%
  filter(year > 2009)

# Democracy-Oriented Donors
ch.0.dem = felm(dem_ngo.log ~ dir.sig.event.l | year + dyad_id | 0 |r_iso, channel %>% filter(r_small == 0 & DOR == 1))
ch.1.dem = felm(dem_bypass.log ~ dir.sig.event.l | year + dyad_id | 0 |r_iso, channel %>% filter(r_small == 0 & DOR == 1))
ch.2.dem = felm(dem_r_public.log ~ dir.sig.event.l | year + dyad_id | 0 |r_iso, channel %>% filter(r_small == 0 & DOR == 1))
ch.3.dem = felm(dem_d_public.log ~ dir.sig.event.l | year + dyad_id | 0 |r_iso, channel %>% filter(r_small == 0 & DOR == 1))

# Non-DO donors
ch.0.non = felm(dem_ngo.log ~ dir.sig.event.l | year + dyad_id | 0 |r_iso, channel %>% filter(r_small == 0 & DOR == 0))
ch.1.non = felm(dem_bypass.log ~ dir.sig.event.l | year + dyad_id | 0 |r_iso, channel %>% filter(r_small == 0 & DOR == 0))
ch.2.non = felm(dem_r_public.log ~ dir.sig.event.l | year + dyad_id | 0 |r_iso, channel %>% filter(r_small == 0 & DOR == 0))
ch.3.non = felm(dem_d_public.log ~ dir.sig.event.l | year + dyad_id | 0 |r_iso, channel %>% filter(r_small == 0 & DOR == 0))

channel_models = list(ch.0.dem, ch.1.dem, ch.2.dem, ch.3.dem,
              ch.0.non, ch.1.non, ch.2.non, ch.3.non)

for(i in 1:length(channel_models)){
  number_fe[i] <- length(levels((channel_models)[[i]]$fe$dyad_id[[1]]))
}

# Table A14 
stargazer(list(channel_models), 
          title = "Effect of NGO Laws on Democracy By Channel",
          label = "",
          dep.var.labels = c("NGOs", "Bypass", "Recipient Gov't", "Donor Gov't",  "NGOs", "Bypass", "Recipient Gov't", "Donor Gov't"),
          covariate.labels = c("NGO Law Passage"),
          add.lines = list(c("Year FEs", rep("\\checkmark", 8)), 
                           c("Dyad FEs", rep("\\checkmark", 8)),
                           c("No. Dyad FEs", number_fe)),
          font.size = "scriptsize",
          omit.stat = "ser",
          type = "text",
          no.space = T,
          out = "./tab_NR_channel.tex")

# Coefficient Plot

result_channel <- do.call(rbind, lapply(channel_models, function(mod) {
  ci <- confint(mod)  # 95% CI
  data.frame(
    dep_var = mod$lhs,
    coef = coef(mod)[1],
    ci_lower = ci[1, 1],
    ci_upper = ci[1, 2]
  )
}))

result_channel$type = c(rep("DO", 4), rep("Non-DO", 4))
result_channel$channel = rep(c("NGO", "Bypass", "Recipient Gov't", "Donor Gov't"),2)
result_channel$index = rep(1:4, 2)

jpeg("./fig_A4_channel_coefs.jpg", units = "in", width = 7, height = 5, res = 300)
ggplot(result_channel, aes(x = reorder(channel, desc(index)), y = coef))+
  geom_hline(yintercept = 0, color = "#2b2b2b")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 11, color = "#2b2b2b"),
        axis.text.x = element_text(size = 11, color = "#2b2b2b"),
        legend.position = "top",
        legend.background = element_rect(color = "#2b2b2b", size = 0.5),
        legend.text = element_text(size = 11, color = "#2b2b2b")) +
  geom_point(aes(color = type),position = position_dodge(width = -0.4), size = 3)+
  geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower, color = type),
                position = position_dodge(width = -0.4), 
                size = 1, width = 0)+
  scale_color_manual(values = c("DO" = "darkblue", "Non-DO" = "darkorange2"), name = "")+
  geom_hline(yintercept = 0)+
  coord_flip()  +
  xlab("")+ 
  ylab("Coefficient")
dev.off()


#Code for twists - Nour's contribution/part
#twist 1
library(ggplot2)
library(dplyr)
library(stargazer)
df_twist <- df %>% 
  arrange(dyad_id, year) %>%
  mutate(
    #forces the outcome to be 0 or 1
    event_int = as.integer(as.logical(dir.sig.event)), 
    prev_aid = lag(tot.log)
  ) %>%
  filter(!is.na(prev_aid)) #removing NAs
twist1_model <- glm(event_int ~ prev_aid + civil_war + cpj.ihs, 
                    data = df_twist, 
                    family = binomial(link = "logit"))
#table for coefficients
invisible(stargazer(twist1_model, 
                    type = "latex", 
                    title = "Twist 1: Logistic Regression Predicting NGO Law Enactment",
                    covariate.labels = c("Total Aid (Lagged)", "Civil War", "Press Freedom (CPJ)"),
                    dep.var.labels = "Probability of NGO Law",
                    out = "twist1_table.tex"))

#twist 2
#"commit_tot" is the aid amounts in df
df_twist2 <- df %>%
  mutate(
    #binary indicator: 1 if aid > 0, else 0
    gave_aid = as.integer(commit_tot > 0)
  )
twist2_model <- glm(gave_aid ~ dir.sig.event.l + DOD + civil_war, 
                    data = df_twist2, 
                    family = binomial(link = "logit"))

#generating the table
invisible(stargazer(twist2_model,  
                    type = "latex",
                    header = FALSE,
                    apply.coef = exp,
                    t.auto = FALSE, p.auto = FALSE,  
                    title = "Twist 2: Odds Ratio of Receiving Any Aid After NGO Law",
                    covariate.labels = c("NGO Law (Lagged)", "Donor Democracy (DOD)", "Civil War"),
                    dep.var.labels = "Probability of Receiving Aid",
                    notes = "Table reports Odds Ratios. Values < 1 indicate a decrease in likelihood.",
                    out = "twist2_odds_table.tex"))