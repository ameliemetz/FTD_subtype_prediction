# group comparisons
# demographics tables
# tables counting missing values / sample sizes

# 1, demographics table FTD participants
# 2, demographics table healthy controls
# 3, ANOVA post-hoc tests table (between FTD groups)
# 4, demographics table at UCSF site
# 5, count participant numbers
# 6, missing values table

library(dplyr)
library(readxl)
library(finalfit)
library(lubridate)
setwd("/data/dadmah/metame/NIFD_PLSpaper")

##################################################################################
# prep data

# read data
Info <- read.csv("NIFD_MRI_Measures.csv", stringsAsFactors=TRUE)

# remove odd participant (incoherent data)
Info <- Info[Info$visit_id != '1_S_0020_1', ]

# add PPVT composite score
PPVT_indices <- 35:38  
PPVTsum <- Info[, PPVT_indices]
sumPPVT <- rowSums(PPVTsum, na.rm = FALSE)
Info$sum_PPVT <- sumPPVT

# select behavioral indices and convert to a numeric matrix
ind_beh <- c(13, 14, 16, 18, 19, 20, 22, 24:26, 28:32, 252, 8:10) #maximal model
data_subset <- as.matrix(Info[, ind_beh])

# replace negative values with NaN
data_subset[data_subset < 0] <- NA
Info[, ind_beh] <- as.data.frame(data_subset)

# select baseline visit
Info <- Info %>%
  group_by(LONI_ID) %>%
  mutate(visit = row_number()) %>%
  ungroup()
Base <-  Info[Info$visit==1,]

# pick variables to include
variables = c('AGE','GENDER','EDUCATION', #demographics
              'CDR_BOX_SCORE','CDR_LANG','CDR_BEHAV', #CDR
              'MMSE_TOT', #MMSE
              'TRCOTOT','CORR30','CORR10','RECOG', #CVLT
              'DIGITFW','DIGITBW', #digit span
              'MTCORR','MTERROR', #modified trail making test
              'DCORR','ANCORR', #verbal fluency
              'BNTCORR', #boston naming test
              'sum_PPVT') #PPVT


##################################################################################
# 1
# demographics table FTD participants

# pick FTD participants
Base_FTD <- Base[Base$DX != "CON",]

Demo_Table <- Base_FTD %>% 
  summary_factorlist("DX", variables,
                     p=TRUE, na_include=TRUE)

# save as table
write.csv(Demo_Table,'NIFD_PLS_demotable.csv',row.names = FALSE, quote = FALSE)

##################################################################################
# 2
# demographics table with healthy controls

Demo_Table <- Base %>% 
  summary_factorlist("DX", variables,
                     p=TRUE, na_include=TRUE)

# save as table
write.csv(Demo_Table,'NIFD_PLS_demotable_controls.csv',row.names = FALSE, quote = FALSE)


##################################################################################
# 3
# post-hoc tests table (between FTD groups)

results <- data.frame(Variable = character(),
                      Comparison = character(),
                      `p-value` = numeric(),
                      stringsAsFactors = FALSE)

for (var in variables) {
  model <- aov(as.formula(paste(var, "~ DX")), data = Base_FTD)
  posthoc <- TukeyHSD(model)
  posthoc_results <- as.data.frame(posthoc$DX)
  comparisons <- rownames(posthoc_results)
  results <- rbind(results,
                   data.frame(Variable = rep(var, length(comparisons)),
                              Comparison = comparisons,
                              p = posthoc_results$`p adj`,
                              stringsAsFactors = FALSE))}
results$p <- round(results$p, 3)
results <- results %>%
  mutate(pval = case_when(p < 0.001 ~ '< 0.001*', p < 0.05 ~ paste0(p, "*"), TRUE ~ as.character(p)))

write.csv(results,'NIFD_demo_posthoc.csv',row.names = FALSE,quote = FALSE)


##################################################################################
# 4
# demographics table at UCSF site

# read data (saved in MATLAB, PLS code)
SITE <- read.csv("table_min_subjects.csv", stringsAsFactors=TRUE)

# combine DX and site column in one
SITE <- SITE %>%
  mutate(group = ifelse(DX == "BV", paste0("bvFTD_", SITE),
                        ifelse(DX == "SV", paste0("svPPA_", SITE),
                               ifelse(DX == "PNFA", paste0("nfvPPA_", SITE), 
                                      paste0("control_", SITE)))))

SITE_Table <- SITE %>% 
  summary_factorlist("group", variables,
                     p=TRUE, na_include=TRUE)

write.csv(SITE_Table,'NIFD_PLS_SITEtable.csv',row.names = FALSE, quote = FALSE)

cdr_lang_summary <- SITE %>%
  group_by(group) %>%
  summarise(
    CDR_LANG_mean = mean(CDR_LANG, na.rm = TRUE),
    CDR_LANG_sd = sd(CDR_LANG, na.rm = TRUE),
    .groups = 'drop'
  )

bv <- SITE[SITE$DX=='BV',]
bv %>% 
  summary_factorlist("SITE", variables,
                     p=TRUE, na_include=TRUE)


##################################################################################
# 5
# count participant numbers

# read full NIFD table
NIFD <- read_excel("NIFD_Clinical_Data_20200724_updated.xlsx")

# remove odd participant
NIFD[NIFD$LONI_ID == '1_S_0020' & NIFD$GENDER == 2, ] <- NA

# pick baseline visit
NIFD <- NIFD %>%
  group_by(LONI_ID) %>%
  mutate(visit = row_number()) %>%
  ungroup()
NIFD <-  NIFD[NIFD$visit==1,]

# print diagnostic groups and sizes
print('FTLDNI all participants')
table(NIFD$DX)

# exclude diagnoses we will not look at
NIFD <-  NIFD[NIFD$DX%in%c("BV","CON","PNFA","SV"),]

# check if all participants have MRI scans 
directory_path <- "/data/dadmah/metame/NIFD_MRI"
folder_names <- list.dirs(directory_path, full.names = FALSE, recursive = FALSE)
NIFD <-  NIFD[NIFD$LONI_ID %in% folder_names,]

# print number of participants that have MRI scans
print('FTLDNI MRI available')
table(NIFD$DX)

##################################################################################
# 6
# missing values table

# pick follow-up visit
# calculate exact age
# keep visit closest to one year

Info$DOB <- Info$DOB <- mdy(Info$DOB)
Info$DOB <- Info$DOB - years(100)
Info$CLINICAL_LINKDATE <- mdy(Info$CLINICAL_LINKDATE)
Info$age_days <- interval(Info$DOB, Info$CLINICAL_LINKDATE) %/% days(1)

Info <- Info %>%
  group_by(LONI_ID) %>%
  mutate(visit = row_number()) %>%
  ungroup()

ids <- unique(Info$LONI_ID)
Info$TimefromBaseline <- NA
Info$N_visits <- NA
for (i in seq_along(ids)) {
  ind <- which(Info$LONI_ID == ids[i])
  Info$TimefromBaseline[ind] <- Info$age_days[ind] - min(Info$age_days[ind])
  Info$N_visits[ind] <- length(ind)
}

Info$keep <- ifelse(Info$visit==1,1,
                    ifelse(Info$N_visits < 3 & Info$visit != 1, 2, 0))

for (i in seq_along(ids)) {
  ind <- which(Info$LONI_ID == ids[i])
  if (Info$N_visits[ind[1]] > 2) {
    closest_index <- ind[which.min(abs(Info$TimefromBaseline[ind] - 365))]
    #  Info$keep[Info$LONI_ID == ids[i]] <- 1
    Info$keep[closest_index] <- 2  
  }
}

# write a table counting numbers of available and missing values per group and visit

ind_beh <- c(12, 13, 14, 17:19, 21, 23:25, 27:31, 251, 7:9) #maximal model

NAs <- Info %>%
  group_by(DX, visit) %>%
  summarize(total_avail = sum(rowSums(pick(all_of(ind_beh)) %>% is.na()) == 0),
            total_visit =   n_distinct(LONI_ID),
            across(all_of(ind_beh), 
                   list(avail = ~ sum(!is.na(.)),
                        missing = ~ sum(is.na(.))), 
                   .names = '{col}_{fn}'), 
            .groups = 'drop') %>%
  mutate(total = ifelse(DX=='BV',70,ifelse(DX=="SV",36,ifelse(DX=='PNFA',30,133)))) %>%
  mutate(across(
    ends_with("_missing"),
    ~ round( . *100/ total_visit, 2 ),  
    .names = "{col}_NAperc"  
  ))

# we only use baseline and 1 year follow-up visit

NAs <- NAs[NAs$visit<3,]

# add percentage of missing values

base_names <- unique(sub("_missing$", "", grep("_missing$", names(NAs), value = TRUE)))

for (name in base_names) {
  NAs[[name]] <- paste0(NAs[[paste0(name, "_missing")]], 
                        " (", NAs[[paste0(name, "_missing_NAperc")]], "%)")
}

# write a nice table
missing <- NAs[base_names]
missing$DX <- NAs$DX
missing$visit <-  NAs$visit
missing <- missing %>%
  mutate(visit = ifelse(visit==1, 'baseline','1-year follow-up'))
missing$total <- NAs$total_visit
missing$total_avail <- NAs$total_avail
missing_transposed <- as.data.frame(t(missing))
colnames(missing_transposed) <- missing_transposed[20, ]
missing_transposed <- missing_transposed[-c(17:20), ]  

# save
write_csv(missing_transposed,'missingvalues_table.csv')