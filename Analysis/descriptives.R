library(matrixStats)
library(data.table)
library(readr)

# Save Location
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
load(paste0(dir_data_qd,"2_all_qd.RData"))
cohort_tmp <- setDT(new_data$w)
zip.tmp <- setDT(new_data$x)

mdat <- with(zip.tmp, data.frame(zip = zip, year = year, mw = 1 - (regionNORTHEAST + regionSOUTH + regionWEST),
                                 ne = regionNORTHEAST, south = regionSOUTH, west = regionWEST))

cohort <- merge(cohort_tmp, mdat, by = c("zip","year"), all.x = TRUE)
w <- cohort$time_count
zip <- merge(zip.tmp, aggregate(w, by = list(zip = cohort$zip, year = cohort$year), sum),
             by = c("zip", "year"), all.x = TRUE)
w2 <- zip$x

# ZIP Code Statistics

## write a function to add a row for continuous variables

add_cont <- function(x, w = rep(1, length(x)), nm_var, nm_level, ndig = 2) {
  
    n_x<-sum(!is.na(x))
    mean_x<-round(weightedMean(x = x, w = w, na.rm = T),ndig)
    sd_x<-round(weightedSd(x = x, w = w, na.rm=T), ndig)
    return(c(nm_var, nm_level, n_x, mean_x, sd_x))
  
}

## write a function to add a row for a level of a categorical variable
add_cat <- function(x, w = rep(1, length(x)), nm_var, nm_level, ndig = 2){

    n_x<-sum(!is.na(x))
    n_x1<-sum(w[x==1], na.rm=T)
    pct_x1<-round(100*n_x1/sum(w), ndig)
    return(c(nm_var, nm_level, n_x, n_x1, pct_x1))
  
}

## write a function to add a death-row
add_death <- function(x, w = rep(1, length(x)), nm_var){
  
  n_x<-sum(!is.na(x))
  n_x1<-sum(x, na.rm=T)
  pct_x1<-round(100*n_x1/sum(w), ndig)
  return(c(nm_var, nm_level = '', n_x, n_x1, pct_x1))
  
}

# Table 1 -----------------------------------------------------------------

## sex ##
table1<-c(add_cat(x=as.numeric(cohort[,female]==1),w = w, nm_var='Female', nm_level=''))

## age ##
table1<-rbind(table1, add_cat(x=as.numeric(cohort[,age_break]=="[65,75)"),
                              w = w, nm_var='Age', nm_level='65-74'))
table1<-rbind(table1, add_cat(x=as.numeric(cohort[,age_break]=="[75,85)"),
                              w = w, nm_var='Age', nm_level='75-84'))
table1<-rbind(table1, add_cat(x=as.numeric(cohort[,age_break]=="[85,95)"),
                              w = w, nm_var='Age', nm_level='85-94'))
table1<-rbind(table1, add_cat(x=as.numeric(cohort[,age_break]=="[95,125)"),
                              w = w, nm_var='Age', nm_level='95+'))

## race ##
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==1),
                               w = w, nm_var='Race',nm_level='White')))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==2),
                               w = w, nm_var='Race',nm_level='Black')))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==5),
                               w = w, nm_var='Race',nm_level='Hispanic')))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==4),
                               w = w, nm_var='Race',nm_level='Asian')))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==3),
                               w = w, nm_var='Race',nm_level='Other or Unknown')))

## dual eligible ##
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,dual]==1), 
                               w = w, nm_var='Medicaid Eligible',nm_level='')))

## region ##
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,ne]==1), 
                               w = w, nm_var='Region',nm_level='Northeast')))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,south]==1), 
                               w = w, nm_var='Region',nm_level='South')))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,mw]==1), 
                               w = w, nm_var='Region',nm_level='Midwest')))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,west]==1), 
                               w = w, nm_var='Region',nm_level='West')))

## deaths ##
table1<-rbind(table1,c(add_death(x=as.numeric(cohort[,dead]), w = w, nm_var = "Deaths")))
table1<-rbind(table1,c(add_death(x=as.numeric(cohort[,time_count]), w = w, nm_var = "Person-years")))

colnames(table1) <- c("Variable", "Level", "Rows", "N", "%")

write_csv(data.frame(table1), path = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/table1.csv')

# Table 2 -----------------------------------------------------------------

## demographics and health ##
table2 <- c(add_cont(x=as.numeric(zip[,pm25]),nm_var='PM2.5', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,mean_bmi]),
                                 w = w2, nm_var='Average BMI', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,smoke_rate]),
                                 w = w2, nm_var='Smoking Rate', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,hispanic]),
                                 w = w2, nm_var='Percent Hispanic', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,pct_blk]),
                                 w = w2, nm_var='Percent Black', nm_level=''))

## socioeconomic variables ##
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,popdensity]),
                                 w = w2, nm_var='Population Density', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,medhouseholdincome]),
                                 w = w2, nm_var='Median Household Income', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,medianhousevalue]),
                                 w = w2, nm_var='Median House Value', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,poverty]),
                                 w = w2, nm_var='Poverty', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,education]),
                                 w = w2, nm_var='Less than High-School', nm_level=''))

table2 <- rbind(table2, add_cont(x=as.numeric(zip[,pct_owner_occ]),
                                 w = w2, nm_var='Percenent Owner Occupied', nm_level=''))

## temperature/humidity ##
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,summer_tmmx]),
                                 w = w2, nm_var='Summer Temperature', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,winter_tmmx]),
                                 w = w2, nm_var='Winter Temperature', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,summer_rmax]),
                                 w = w2, nm_var='Summer Humidity', nm_level=''))
table2 <- rbind(table2, add_cont(x=as.numeric(zip[,winter_rmax]),
                                 w = w2, nm_var='Winter Humidity', nm_level=''))

colnames(table2) <- c("Variable", "Level", "Rows", "Mean", "SD")

write_csv(data.frame(table2), path = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/table2.csv')
