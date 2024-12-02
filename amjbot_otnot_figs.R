# amjbot_otnot_figs
## script to create figures for American Journal of Botany "On the Nature of Things" manuscript on photosynthetic acclimation

## load up packages
library(dplyr)
library(zoo)
library(lubridate)
library(plantecophys)
library(R.utils)
library(ggplot2)
library(grid)

## load up model functions (can be sourced from github.com/smithecophyslab/optimal_vcmax_R)
# source('../optimal_vcmax_R/calc_optimal_vcmax.R')
# sourceDirectory('../optimal_vcmax_R/functions')

## load NEON data (can be sourced from NEON data portal)
# HF_temp = read.csv('/Users/nicksmith/Documents/Research/Timescale_acclimation/Sensitivity/NEON_HarvardForest/NEON_temp-air-single/stackedFiles/SAAT_30min.csv')
# HF_par = read.csv('/Users/nicksmith/Documents/Research/Timescale_acclimation/Sensitivity/NEON_HarvardForest/NEON_par/stackedFiles/PARPAR_30min.csv')
# HF_rh = read.csv('/Users/nicksmith/Documents/Research/Timescale_acclimation/Sensitivity/NEON_HarvardForest/NEON_rel-humidity/stackedFiles/RH_30min.csv')
nrow(HF_temp)
nrow(HF_par)
nrow(HF_rh)

###############################################################################################
## 1. figure showing photosynthesis and costs of acclimating to different time scales
###############################################################################################

## join average data together
HF_temp_mean = HF_temp[, 5:7]
HF_par_mean = HF_par[, 5:7]
HF_rh_mean = HF_rh[, 5:7]
HF_mean = left_join(HF_temp_mean, HF_rh_mean)
HF_mean = left_join(HF_mean, HF_par_mean)

## calculate VPD from T and rh
vpd_from_t <- function(t, rh){ 
	
	svp = 610.7 * 10^(7.5*t/(237.3+t))
	vpd = ((100-rh) / 100)*svp
	vpd #Pa
	
}

HF_mean$VPD = vpd_from_t(HF_mean$tempSingleMean, HF_mean$RHMean) / 1000

## remove NAs
HF_mean <- subset(HF_mean, VPD!='NA')

## keep only daytime values
HF_mean_day = subset(HF_mean, PARMean > 0)

## aggregate by day
### separate out dates and times
### see: http://www.neonscience.org/dc-convert-date-time-POSIX-r
HF_mean_day$date = as.factor(as.Date(HF_mean_day$endDateTime)) 
### aggregate
HF_group_by_date = group_by(HF_mean_day[, 2:7], date)
HF_mean_date = summarize(HF_group_by_date, 
                         par = mean(PARMean, na.rm = T), 
                         temp = mean(tempSingleMean, na.rm = T), 
                         vpd = mean(VPD, na.rm = T))

## calculate running means
### previous 30 minutes: i = 1
### previous 3 months: i = 4320 (48 * 90 days)
MA <- list()
for (i in 1:90){
	
	temp_MA = rollapply(HF_mean_date$temp, i, mean, fill = NA, align ='right')
	par_MA = rollapply(HF_mean_date$par, i, mean, fill = NA, align ='right')
	vpd_MA = rollapply(HF_mean_date$vpd, i, mean, fill = NA, align ='right')
	
	temporary = cbind(temp_MA, par_MA, vpd_MA)
	colnames(temporary) = c(paste('temp_MA', i, sep = '_'), paste('par_MA', i, sep = '_'), paste('vpd_MA', i, sep = '_'))
	
	MA[[i]] = temporary
	
}

## calculate optimal traits throughout the summer (June to September)
days_interest = 515:636 # June through September
acclimation_values = list()
for (i in 1:90){
	
	vars = calc_optimal_vcmax(cao = 400, 
	                              tg_c = MA[[i]][days_interest, 1], 
	                              paro = MA[[i]][days_interest, 2], 
	                              z = 300, 
	                              vpdo = MA[[i]][days_interest, 3])
	
	vars_list = cbind(vars, paste(i - 1))
	
	acclimation_values[[i]] = vars_list
	
}

## calculate actual photosynthesis
HF_mean_date_2017 = HF_mean_date[days_interest, ]
HF_patm = (acclimation_values[[1]]$patm / 1000)[1]

photosynthesis_list = list()
for (i in 1:90){
	
	photosynthesis = Photosyn(VPD = HF_mean_date_2017$vpd, # vpd on day of interest
	                          Ca = 400, 
	                          PPFD = HF_mean_date_2017$par, # par on day of interest
	                          Tleaf = HF_mean_date_2017$temp, # temp on day of interest
	                          Patm = HF_patm, # patm on day of interest
	                          Jmax = acclimation_values[[i]]$jmax, # run through vcmax and jmax values
	                          Vcmax = acclimation_values[[i]]$vcmax)
	
	photosynthesis_list[[i]] = photosynthesis
	
}

## get seasonal photosynthesis total
photosynthesis_season = c()
for (i in 1:90){
	
	temp = mean(photosynthesis_list[[i]]$ALEAF)
	photosynthesis_season = c(photosynthesis_season, temp)
	
}

plot(photosynthesis_season ~ seq(1, 90, 1))

## calculate the cost of acclimation as the variability in Vcmax
vcmax_var = c()
for (i in 1:90){
	
	var = sd(acclimation_values[[i]]$vcmax)
	
	vcmax_var = c(vcmax_var, var)
	
}

plot(vcmax_var ~ seq(1, 90, 1))

## make a pretty plot with them together
timescale_dataframe <- data.frame(cbind(photosynthesis_season, vcmax_var))
timescale_plot <- ggplot(data = timescale_dataframe, aes(y = photosynthesis_season/photosynthesis_season[1], x = seq(1,90,1))) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(aes(color='blue'), linewidth = 2) +
  geom_line(aes(y=vcmax_var/vcmax_var[1], color='red'), linewidth = 2) +
  ylab('Relative seasonal value') +
  xlab('Acclimation timescale (days)') +
  scale_colour_manual(name = NULL, values =c('blue'='blue','red'='red'), labels = c('Photosynthesis','Cost')) +
  ylim(c(0.2, 1)) +
  xlim(c(0, 90)) +
  labs(tag = "A")

###############################################################################################
## 2. figure showing 1 week of diurnal photosynthesis of average and midday acclimation
###############################################################################################

## calculate max values for acclimation
HF_mean_day$date = as.factor(as.Date(HF_mean_day$endDateTime)) 
HF_group_by_date = group_by(HF_mean_day[, 2:7], date)
HF_max_date = summarize(HF_group_by_date, 
                         par = max(PARMean, na.rm = T), # change needed to get acclimation to max PAR
                         temp = mean(tempSingleMean, na.rm = T), 
                         vpd = mean(VPD, na.rm = T))

## calculate running means
### previous 30 minutes: i = 1
### previous 3 months: i = 4320 (48 * 90 days)
MA_max <- list()
for (i in 1:90){
  
  temp_MA = rollapply(HF_max_date$temp, i, mean, fill = NA, align ='right')
  par_MA = rollapply(HF_max_date$par, i, mean, fill = NA, align ='right')
  vpd_MA = rollapply(HF_max_date$vpd, i, mean, fill = NA, align ='right')
  
  temporary = cbind(temp_MA, par_MA, vpd_MA)
  colnames(temporary) = c(paste('temp_MA', i, sep = '_'), paste('par_MA', i, sep = '_'), paste('vpd_MA', i, sep = '_'))
  
  MA_max[[i]] = temporary
  
}

## calculate optimal traits throughout the summer (June to September)
days_interest = 515:636 # June through September
acclimation_values_max = list()
for (i in 1:90){
  
  vars = calc_optimal_vcmax(cao = 400, 
                            tg_c = MA[[i]][days_interest, 1], 
                            paro = MA_max[[i]][days_interest, 2], 
                            z = 300, 
                            vpdo = MA[[i]][days_interest, 3])
  
  vars_list = cbind(vars, paste(i - 1))
  
  acclimation_values_max[[i]] = vars_list
  
}

## make a dataframe with 30-min HF data that includes acclimated vcmax and jmax values
HF_mean$date = as.factor(as.Date(HF_mean$startDateTime)) 
HF_30min_data <- subset(HF_mean, as.Date(date) > as.Date("2017-05-31") & as.Date(date) < as.Date("2017-09-30"))
head(HF_30min_data)
tail(HF_30min_data)

##
HF_30min_data_group <- group_by(HF_30min_data, startDateTime, endDateTime, date)
HF_30min_data_summary <- summarise(HF_30min_data_group,
                                 tmp = mean(tempSingleMean, na.rm = T),
                                 rh = mean(RHMean, na.rm = T),
                                 par = mean(PARMean, na.rm = T),
                                 vpd = mean(VPD, na.rm = T))
head(HF_30min_data_summary)

## create data frame to merge in acclimated values
merge_doys <- data.frame(as.factor(as.Date(seq(as.Date("2017-06-01"), as.Date("2017-09-30"), 1))), days_interest)
colnames(merge_doys) <- c('date', 'doy')

## add doy to 30 min data
HF_30min_data_doy <- left_join(HF_30min_data, merge_doys)
head(HF_30min_data_doy)
tail(HF_30min_data_doy)
nrow(HF_30min_data_doy)

## boil it down to a week
HF_30min_week <- subset(HF_30min_data_doy, doy >545 & doy <553)
tail(HF_30min_week)
nrow(HF_30min_week)

HF_30min_week_group <- group_by(HF_30min_week, startDateTime, endDateTime, date, doy)
HF_30min_week_summary <- summarise(HF_30min_week_group,
                                   tmp = mean(tempSingleMean, na.rm = T),
                                   rh = mean(RHMean, na.rm = T),
                                   par = mean(PARMean, na.rm = T),
                                   vpd = mean(VPD, na.rm = T))
head(HF_30min_week_summary)
nrow(HF_30min_week_summary)

### merge in vcmax and jmax
doy = days_interest

acclimation_values_max_week <- cbind(acclimation_values_max[[7]], doy)
HF_30min_week_max <- left_join(HF_30min_week_summary, acclimation_values_max_week, by = c('doy'))
head(HF_30min_week_max)

acclimation_values_mean_week <- cbind(acclimation_values[[7]], doy)
HF_30min_week_mean <- left_join(HF_30min_week_summary, acclimation_values_mean_week, by = c('doy'))
head(HF_30min_week_mean)

### calculate photosynthesis
photosynthesis_max = Photosyn(VPD = HF_30min_week_max$vpd.x, # vpd on day of interest
                          Ca = 400, 
                          PPFD = HF_30min_week_max$par.x, # par on day of interest
                          Tleaf = HF_30min_week_max$tmp, # temp on day of interest
                          Patm = HF_patm, # patm on day of interest
                          Jmax = HF_30min_week_max$jmax, # run through vcmax and jmax values
                          Vcmax = HF_30min_week_max$vcmax)
plot(photosynthesis_max$ALEAF, type = 'l')

photosynthesis_mean = Photosyn(VPD = HF_30min_week_mean$vpd.x, # vpd on day of interest
                              Ca = 400, 
                              PPFD = HF_30min_week_mean$par.x, # par on day of interest
                              Tleaf = HF_30min_week_mean$tmp, # temp on day of interest
                              Patm = HF_patm, # patm on day of interest
                              Jmax = HF_30min_week_mean$jmax, # run through vcmean and jmean values
                              Vcmax = HF_30min_week_mean$vcmax)
plot(photosynthesis_mean$ALEAF, type = 'l')

## make pretty figure
diurnal_plot <- ggplot(data = photosynthesis_mean, aes(y=ALEAF, x = seq(1, 312, 1))) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(aes(color='black'), linewidth = 1.5, lty = 1) +
  geom_line(data = photosynthesis_max, aes(color='orange'), 
            linewidth = 1.5, lty = 1, alpha = 0.7) +
  ylab(expression('Photosynthesis (µmol m'^'2'*' s'^'-1'*')')) +
  xlab('Time (7 total days)') +
  scale_colour_manual(name = 'Acclimated conditions', 
                      values =c('black'='black','orange'='orange'), 
                      labels = c('Mean','Max')) +
  ylim(c(-10, 50)) +
  labs(tag = "B")

###############################################################################################
## 3. figure showing impact of altered soil resource acquisition costs on 
# stomatal conductance and photosynthetic N
###############################################################################################

## make a beta sensitivity simulation
beta_sensitivity <- calc_optimal_vcmax(beta = seq(10, 1000, 1))
beta_sensitivity$gs <- (beta_sensitivity$Ac/beta_sensitivity$cao) / (1 - beta_sensitivity$chi)

## make a vpd sensitivity simulation
vpd_sensitivity <- calc_optimal_vcmax(vpdo = seq(0.5, 5, 0.1))
vpd_sensitivity$gs <- (vpd_sensitivity$Ac/vpd_sensitivity$cao) / (1 - vpd_sensitivity$chi)

beta_gs_plot <- ggplot(data = beta_sensitivity, aes(y = gs/gs[991], x = beta)) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(color = 'blue', linewidth = 3, lty = 1) +
  ylab(expression('Relative stomatal conductance')) +
  xlab('Cost of acquiring nutrients relative to water (unitless)') +
  ylim(c(0, 1))

tiff(filename = "plots/beta_gs_plot.tiff", 
    width = 9, height = 9, units = 'in', res = 300)
grid.newpage()
grid.draw(beta_gs_plot)
dev.off()

vpd_gs_plot <- ggplot(data = vpd_sensitivity, aes(y = gs/gs[1], x = vpdo)) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(color = 'blue', linewidth = 3, lty = 1) +
  ylab(expression('Relative stomatal conductance')) +
  xlab('Atmospheric vapor pressure deficit (kPa)') +
  ylim(c(0, 1)) +
  xlim(c(0,5))


tiff(filename = "plots/vpd_gs_plot.tiff", 
     width = 9, height = 9, units = 'in', res = 300)
grid.newpage()
grid.draw(vpd_gs_plot)
dev.off()


beta_n_plot <- ggplot(data = beta_sensitivity, aes(y = nphoto/nphoto[1], x = beta)) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(color = 'darkorange', linewidth = 3, lty = 1) +
  ylab(expression('Relative photosynthetic nitrogen')) +
  xlab('Acq. cost of nutrients relative to water (unitless)') +
  ylim(c(0.9, 1))

tiff(filename = "plots/beta_n_plot.tiff", 
     width = 9, height = 9, units = 'in', res = 300)
grid.newpage()
grid.draw(beta_n_plot)
dev.off()

vpd_n_plot <- ggplot(data = vpd_sensitivity, aes(y = nphoto/nphoto[46], x = vpdo)) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(color = 'darkorange', linewidth = 3, lty = 1) +
  ylab(expression('Relative photosynthetic nitrogen')) +
  xlab('Atmospheric vapor pressure deficit (kPa)') +
  ylim(c(0.9, 1)) +
  xlim(c(0,5))

tiff(filename = "plots/vpd_n_plot.tiff", 
     width = 9, height = 9, units = 'in', res = 300)
grid.newpage()
grid.draw(vpd_n_plot)
dev.off()

###################################################################################################
## output plots
###################################################################################################

### Figure 2 (second two combined)
beta_gs_plot_g <- ggplotGrob(beta_gs_plot)
beta_n_plot_g <- ggplotGrob(beta_n_plot)
beta_gs_n_plot_g <- rbind(beta_gs_plot_g, 
                               beta_n_plot_g, 
                              size = "max")
vpd_gs_plot_g <- ggplotGrob(vpd_gs_plot)
vpd_n_plot_g <- ggplotGrob(vpd_n_plot)
vpd_gs_n_plot_g <- rbind(vpd_gs_plot_g, 
                          vpd_n_plot_g, 
                          size = "max")
beta_vpd_gs_n_plot_g <- cbind(beta_gs_n_plot_g,
                              vpd_gs_n_plot_g,
                              size = "max")

# tiff(filename = "/Users/nicksmith/Documents/Research/amjbot_otnot/revision/revision2/plots/Figure2.tiff", 
#      width = 18, height = 18, units = 'in', res = 300)
# grid.newpage()
# grid.draw(beta_vpd_gs_n_plot_g)
# dev.off()

timescale_plot_g <- ggplotGrob(timescale_plot)
diurnal_plot_g <- ggplotGrob(diurnal_plot)
timescale_diurnal_plot_g <- cbind(timescale_plot_g, 
                                  diurnal_plot_g, 
                                  size = "max")

# tiff(filename = "/Users/nicksmith/Documents/Research/amjbot_otnot/revision/revision2/plots/Figure1.tiff", 
#      width = 18, height = 8, units = 'in', res = 300)
# grid.newpage()
# grid.draw(timescale_diurnal_plot_g)
# dev.off()

#####################################
## additional calculations
#####################################
mean(HF_30min_week_max$vcmax) # average vcmax for simulations assuming acclimation to max light
mean(HF_30min_week_mean$vcmax) # average vcmax for simulations assuming acclimation to mean light









