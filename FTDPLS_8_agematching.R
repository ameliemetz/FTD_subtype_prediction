# here we select an age-matched subgroup 
# see https://github.com/AlBRzez/brain_size_sex_brzezinskirittner_2024/blob/main/code/analyses/03_matching.r
# based on nfvPPA group (older than others)

library(lubridate)
library(dplyr)

setwd("/data/dadmah/metame/NIFD_PLSpaper")

# read data
Info <- read.csv("NIFD_MRI_Measures.csv", stringsAsFactors=TRUE)

# remove odd participant (incoherent data)
Info <- Info[Info$visit_id != '1_S_0020_1', ]

# add PPVT composite score
PPVT_indices <- 35:38  
PPVTsum <- Info[, PPVT_indices]
sumPPVT <- rowSums(PPVTsum, na.rm = FALSE)
Info$sum_PPVT <- sumPPVT

# save original table
Info_all <- Info

# select behavioral indices and convert to a numeric matrix
ind_beh <- c(13, 14, 16, 18, 19, 20, 22, 24:26, 28:32, 252, 8:10) #maximal model
#ind_beh <- c(13, 14, 16, 32, 8:10) #minimal model
data_subset <- as.matrix(Info[, ind_beh])

# replace negative values with NaN
data_subset[data_subset < 0] <- NA
Info[, ind_beh] <- as.data.frame(data_subset)

# remove participants with missing values in ind_beh columns
Info <- Info[rowSums(is.na(data_subset)) <= 1, ]

# calculate exact age
Info$DOB <- Info$DOB <- mdy(Info$DOB)
Info$DOB <- Info$DOB - years(100)
Info$CLINICAL_LINKDATE <- mdy(Info$CLINICAL_LINKDATE)
Info$age_months <- interval(Info$DOB, Info$CLINICAL_LINKDATE) %/% months(1)
Info$age_years <- interval(Info$DOB, Info$CLINICAL_LINKDATE) %/% years(1)

#pick baseline visit
Info <- Info %>%
  group_by(LONI_ID) %>%
  mutate(visit = row_number()) %>%
  ungroup()
Info <-  Info[Info$visit==1,]

# split data into diagnostic groups
nfv <- Info %>% filter(DX == "PNFA")  # nfvPPA
bv <- Info %>% filter(DX == "BV")  # bvFTD
sv <- Info %>% filter(DX == "SV")  # svPPA
con <- Info %>% filter(DX == "CON")  # healthy controls

# matching function
get_age_sample <- function(distance, nfv, bv, sv) {
  matc <- data.frame()  
  
  for (i in 1:nrow(nfv)) {
    # find age-matched participants from each group
    bv_match <- bv %>% filter(abs(nfv$age_months[i] - age_months) < (distance * nfv$age_months[i]))
    sv_match <- sv %>% filter(abs(nfv$age_months[i] - age_months) < (distance * nfv$age_months[i]))
    #con_match <- con %>% filter(abs(nfv$age_months[i] - age_months) < (distance * nfv$age_months[i]))
    
    # randomly select one participant from each matching group
    if (nrow(bv_match) > 0) bv_selected <- bv_match %>% slice_sample(n = 1) else next
    if (nrow(sv_match) > 0) sv_selected <- sv_match %>% slice_sample(n = 1) else next
    #if (nrow(con_match) > 0) con_selected <- con_match %>% slice_sample(n = 1) else next
    
    # combine matched participants 
    mat <- rbind(nfv[i, ], bv_selected, sv_selected)
    mat$round <- i  
    
    # add to final data frame
    matc <- rbind(matc, mat)
    
    # remove selected participants from their respective groups to prevent reuse
    bv <- bv %>% filter(!row_number() %in% which(bv$age_months == bv_selected$age_months))
    sv <- sv %>% filter(!row_number() %in% which(sv$age_months == sv_selected$age_months))
    #con <- con %>% filter(!row_number() %in% which(con$age_months == con_selected$age_months))
  }
  
  print(dim(matc))  
  return(matc)  
}

# call the function with a small distance (e.g., 0.002)
age_mat <- get_age_sample(.05, nfv, bv, sv)
age_mat %>%
  group_by(DX) %>%
  summarise(mean_age_years = sd(age_years, na.rm = TRUE))

# plot entire group versus age matched group

col = c("BV"="#1982C4","SV"="#8AC926","PNFA"="#FF595E","CON"="#6A4C93")

comp_plot <- function(df, title) {
  p1 <- 
    df |>
    
    ggplot(aes(x = age_years, color = DX)) +
    #facet_wrap(~DX, scales = "free") +
    # facet_grid(sex~name, scales = "free_x") +
    geom_histogram(alpha = .3, position = "stack", fill = NA) +
    
    labs(title = paste("age-matched", title),
         subtitle = paste(nrow(df), "subjects :
",nrow(df)/3, "per group")) +
    scale_color_manual(values = col) +
    ylim(0,15)+
    theme_light() +
    theme(
      legend.position = "none"
    )
  p2 <- 
    Info |>
    
    ggplot(aes(x = age_years, color = DX)) +
  #  facet_wrap(~DX, scales = "free") +
    # facet_grid(sex~name, scales = "free_x") +
    geom_histogram(alpha = .3, position = "stack", fill = NA) +
    
    labs(title = paste("all subjects", title),
         subtitle = paste(nrow(Info), "subjects: 
bvFTD =",nrow(Info[Info$DX=='BV',]),"svPPA =",nrow(Info[Info$DX=='SV',]),"nfvPPA =",nrow(Info[Info$DX=='PNFA',]))) +
    scale_color_manual(values = col,
                       aesthetics = c("color", "fill")) +
    ylim(0,15)+
    theme_light() +
    theme(
      legend.position = "none"
    )
  cowplot::plot_grid(p1, p2, nrow = 1)
}

comp_plot(age_mat, " ")

# save as table
matched <- age_mat$LONI_ID
Info_all <- Info_all[Info_all$LONI_ID %in% matched,]
write_csv(Info_all,'Info_agematched.csv')
