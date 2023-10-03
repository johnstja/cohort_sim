##############################################
#### Identify tree establishment cohorts  ####
##############################################

#Author:  James Johnston

#Load some libraries
library(tidyverse)
library(KernSmooth)
library(ggplot2)

#Read the establishment data for the upper Middle Fork study sites (sites 1-16) and the MAR and TUM demonstration sites (sites 17-18)
establishment_data <- structure(list(site=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18), estab_year=c(1572,1571,1571,1573,1576,1575,1576,1578, 1583,1584,1594,1594,1595,1599,1601,1609,1724,1639,1649, 1668,1667,1670,1670,1679,1690,1698,1704,1774,1784,1784,1803,1614,1623,1648,1657,1657,1658,1682,1701,1707,1708,1710,1762,1637,1667,1668,1669,1671,1675,1686,1712,1761,1845,1851,1581,1590,1607,1602,1608,1608,1611,1614,1623,1630,1627,1644,1615,1635,1652,1668,1664,1690,1690,1690,
1694,1699,1716,1729,1411,1468,1471,1558,1570,1661,1675,1691,1764,1799,1800,1802,1818,1509,1519,1661,1679,1687,1733,1747,1747,1746,1753,1799,1810,1588,1624,1623,1630,1627,1633,1631,1634,1642,1681,1699,1616,1633,1669,1676,1689,1716,1727,1737,1761,1768,1860,1865,1248,1498,1512,1514,1513,1513,1515,1517,1521,1521,1528,1557,1562,1568,1598,1607,1743,1718,1723,1733,1731,1728,1731,1732,1736,1736,1737,1743,1744,1740,1744,1746,1748,1747,1753,1760,1754,1528,1550,1575,1578,1576,1601,1666,1669,1696,1788,1871,1920,1545,1597,1602,1610,1658,1662,1672,1672,1680,1840,1852,1530,1533,1534,1531,1564,1569,1570,1573,1569,
1570,1572,1574,1576,1668,1668,1715,1795,1539,1540,1541,1547,1548,1552,1565,1590,1660,1659,1659,1662,1866,1872,1878,1883,1889,1898,1908,1918,1923,1932,1941,1944,1950,1958,1962,1862,2010,2010,2011,2011,2011,2012,2012,2012,2012,2012,2013,2013,2013,2013)), row.names=c("15","4","6","5","1","13","12","11","8","7","16","17","10","3","2","9","14","18","27","28","20","22","184","24","21","23","19","26","182","183","25","30","29","33","34","40","32","37","31","36","35","39","38","42","43","47","41","46","48","49","44","181","180","45","50","56","58","59","57","54","61","51","53","52","60","55","64","63","69","62","190","67","70","66","68","191","71","65","84","78","76","77","83","72","73","74","80","82","81","75","79","85","88","87","89","95","90","86","94","91","92","93","188","96","100","104","105","99","102","103","189","101","97","98","107","201","106","196","193","198","109","108","202","200","197","199","119","124","111","118","121","117","122","113","114","123","120","110","194","115","112","116","195","137","140","136","138","144","139","142","130","131","129","134","143","132","125","126","128","141","133","127","135","153","187","148","152","151","146","147","145","185","150","186","149","156","160","205","206","154","159","155","158","157","204","203","165","163","166","210","192","208","209","167","207","161","162","212","211","213","214","215","164","174","170","171","173","169","172","168","179","177","178","176","175","1100","216","310","410","510","610","710","810","910","1010","1110","1210","1310","1410","1510","1610","1710","1810","1910","2010","217","221","231","241","251","261","271","281","291","301"),class="data.frame")

#Cohort detection function
cohort_function <- function(x, yearsvec, nosims, siglevel){
  estab_vec <- x[,yearsvec]
  estab_vec_pad <- c(min(estab_vec)-50, estab_vec, max(estab_vec)+50)
  estab <- if(max(estab_vec)-min(estab_vec) >= 150) 
    estab_vec else estab_vec_pad
  minest <- min(estab)
  maxest <- max(estab)
  notrees <- length(estab)
  nosims <- nosims
  unif_resampfunct <- function(){
    simestab <- sample(minest:maxest, size = notrees, replace = TRUE)
    return(simestab)
  }
  sim_data <- data.frame(replicate(n = nosims, 
    expr = unif_resampfunct()))
  names(sim_data) <- gsub(x = names(sim_data), pattern = "\\X", replacement = "sim") 
  all_data <- cbind(sim_data, estab)
  bw <- dpik(estab)
  all_density <- apply(all_data, 2, bkde, bandwidth=bw/2.75)  
  long_density <- do.call(rbind.data.frame, all_density)
  long_density$run <- rownames(long_density)
  long_density$group <- substr(long_density$run, 1, 3)
  sim_dens <- long_density[(long_density$group=="sim"),]
  est_dens <- long_density[(long_density$group=="est"),]
  sim_dens_sig <- quantile(sim_dens$y, siglevel)
  if(any(est_dens$y >= sim_dens_sig)){    
    estab_sig <- est_dens[(est_dens$y >= sim_dens_sig),]
    estab_sig$lower_bound <- rep(sim_dens_sig, nrow(estab_sig))
    estab_sig <- estab_sig[,c("x", "y", "lower_bound")]
    estab_sig$run_no <- as.numeric(sub(".*b.", "", 
        rownames(estab_sig)))
    rownames(estab_sig) <- 1:nrow(estab_sig)
    estab_sig$rows <- as.numeric(rownames(estab_sig))
    ints <- c(0, which(diff(estab_sig$run_no) != 1), length(estab_sig$run_no))
    estab_sig$cohort_no <- cut(estab_sig$rows, breaks=ints, labels=FALSE)
    estab_sig$rows <- NULL
  } else {
    estab_sig <- data.frame(x = 0, y = 0, lower_bound = 0, run_no = 0, cohort_no="NA")
  }
  return(list(crit_value=sim_dens_sig, sim_dens=sim_dens, est_dens=est_dens, estab_sig=estab_sig))
}

#Split the establishment data into a list split by site 
site_list <- split(establishment_data, establishment_data$site)

#Run the cohort detection function on every element of the list (run the function for each site).  "no-sims" sets the number of simulated establishment years (for each site) and siglevel sets the critical threshold for evaluating significant of actual establishment data.  1,000 simulations will yield a reasonable estimate of significant cohorts and will run on most computers in less than 20 seconds.  We ran 10,000 simulations per site when completing statistical and graphical analysis for our manuscript.  
set.seed(123)
results_list <- map(site_list, cohort_function, yearsvec="estab_year", nosims=1000, siglevel=.99)

#View the resulting list structure to get an idea how the function reports out results
data.tree::FromListSimple(results_list)

#Create a dataframe with critical values for cohort significance
crit_values <- unname(unlist(map(results_list, 1)))
crit_values_df <- data.frame(site=as.character(1:length(crit_values)), crit_value=crit_values) 

#Create data frame of real establishment determined to represent tree cohorts (significant given critical threshold). 
estab_sig <- map(results_list, 4)
estab_sig <- do.call(rbind.data.frame, estab_sig)
estab_sig$site <- substr(rownames(estab_sig), 1, 2)
estab_sig$site <- gsub("[[:punct:]]", "", estab_sig$site)

#Some necessary housekeeping:  Remove the critical values and significant cohorts data from the results list
results_list <- lapply(results_list, function(x) x[-1])
results_list <- lapply(results_list, function(x) x[-3])

#Create data frame with establishment and simulated establishment kernel density estimate curves
est_sim <- do.call(rbind, unlist(results_list, recursive = FALSE))
est_sim$site <- substr(rownames(est_sim), 1, 2)
est_sim$site <- gsub("[[:punct:]]", "", est_sim$site)
rownames(est_sim) <- NULL

#Split this dataframe into actual and simulated establishment
est_dens <- est_sim %>% filter(group=="est")
sim_dens <- est_sim %>% filter(group=="sim")

#Subset the demonstration site data (to graph separately) and then remove it
demo_estab_sig <- estab_sig %>% filter(site==17 | site==18)
estab_sig <- estab_sig %>% filter(!site %in% demo_estab_sig$site)
demo_sim_dens <- sim_dens %>% filter(site==17 | site==18)
sim_dens <- sim_dens %>% filter(!site %in% demo_sim_dens$site)
demo_est_dens <- est_dens %>% filter(site==17 | site==18)
est_dens <- est_dens %>% filter(!site %in% demo_est_dens$site)
demo_crit_values_df <- crit_values_df %>% filter(site==17 | site==18)
crit_values_df <- crit_values_df %>% filter(!site %in% demo_crit_values_df$site)

#Give the demonstration sites their appropriate names back
demo_estab_sig$site <- ifelse(demo_estab_sig$site==17, "Mar", "Tum")
demo_sim_dens$site <- ifelse(demo_sim_dens$site==17, "Mar", "Tum")
demo_est_dens$site <- ifelse(demo_est_dens$site==17, "Mar", "Tum")
demo_crit_values_df$site <- ifelse(demo_crit_values_df$site==17, "Mar", "Tum")

#Order the Middle Fork site factors to match the order in manuscript Figure 6
estab_sig$site <- factor(estab_sig$site, levels=c("12", "5", "11", "16", "1", "6", "8", "13", "9", "15", "3",  "7", "14", "2", "10", "4"))
sim_dens$site <- factor(sim_dens$site, levels=c("12", "5", "11", "16", "1", "6", "8", "13", "9", "15", "3",  "7", "14", "2", "10", "4"))
est_dens$site <- factor(est_dens$site, levels=c("12", "5", "11", "16", "1", "6", "8", "13", "9", "15", "3",  "7", "14", "2", "10", "4"))
crit_values_df$site <- factor(crit_values_df$site, levels=c("12", "5", "11", "16", "1", "6", "8", "13", "9", "15", "3",  "7", "14", "2", "10", "4"))

#Create graphical file for the Middle Fork sites output (but don't plot yet)
rig_sites_plot <- ggplot() +
  geom_ribbon(data=estab_sig, aes(x=x, ymax=y, ymin=lower_bound, group=cohort_no, fill=cohort_no), fill="darkolivegreen", alpha=0.5) +
  geom_path(data=sim_dens, aes(x=x, y=y), color="snow4", size=.025, alpha=.1) +
  geom_path(data=est_dens, aes(x=x, y=y), color="darkolivegreen", size=.65, alpha=1) +
  geom_hline(data=crit_values_df, aes(yintercept=crit_value), linetype='longdash', color = "black", size=.1) +
  scale_x_continuous("", breaks = seq(1450, 1850, 100), limits=c(1450, 1850)) + 
  scale_y_continuous("Density", breaks = c(0, .04)) + 
  facet_wrap(~ site, ncol = 2, strip.position="right") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(axis.ticks.x = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(2.5,0,0,0))) +
  theme(strip.background = element_blank()) +
  #theme(legend.position="none") +
  theme(panel.grid.minor = element_line(linewidth = 0.1), panel.grid.major = element_line(linewidth = .25)) +
  theme(panel.spacing.y = unit(.35, "lines")) 
#rig_sites_plot

#Create graphical file for the demonstration sites output (but don't plot yet)
demo_sites_plot <- ggplot() +
  geom_ribbon(data=demo_estab_sig, aes(x=x, ymax=y, ymin=lower_bound, group=cohort_no, fill=cohort_no), fill="darkolivegreen", alpha=0.5) +
  geom_path(data=demo_sim_dens, aes(x=x, y=y), color="snow4", size=.025, alpha=.1) +
  geom_path(data=demo_est_dens, aes(x=x, y=y), color="darkolivegreen", size=.65, alpha=1) +
  geom_hline(data=demo_crit_values_df, aes(yintercept=crit_value), linetype='longdash', color = "black", size=.1) +
  geom_vline(data=filter(demo_estab_sig, site=="Tum"), aes(xintercept=2009), size=.25, colour="darkred") + 
  scale_x_continuous("", breaks = seq(1850, 2020, 50), limits=c(1850, 2020)) + 
  scale_y_continuous("Density") + 
  facet_wrap(~ site, ncol = 2, strip.position="right") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(axis.ticks.x = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(2.5,0,0,0))) +
  theme(strip.background = element_blank()) +
  #theme(legend.position="none") +
  theme(panel.grid.minor = element_line(linewidth = 0.1), panel.grid.major = element_line(linewidth = .25)) +
  theme(panel.spacing.y = unit(.35, "lines")) 
#demo_sites_plot

#Tally establishment in the demonstration sites and rename the resulting file
estab_tally <- establishment_data %>%
  filter(site==17 | site==18) %>%
  group_by(site, estab_year) %>%
  count()
estab_tally$site <- ifelse(estab_tally$site==17, "Mar", "Tum")
estab_tally

#Histogram of the demonstration sites (don't plot yet)
agestructurehistogram <- ggplot() +
  geom_bar(data=estab_tally, aes(x=estab_year, y=n), stat = "identity", position = "stack", width=1, fill="darkolivegreen") +
  geom_vline(data=filter(estab_tally, site=="Tum"), aes(xintercept=2009), size=.25, colour="darkred") + 
  scale_x_continuous("", breaks = seq(1850, 2020, 50), limits=c(1850, 2020)) + 
  scale_y_continuous("Count", breaks = seq(0, 6, 3)) +
  facet_wrap(~ site, ncol = 2, strip.position="right") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(axis.ticks.x = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(2.5,0,0,0))) +
  theme(strip.background = element_blank()) +
  #theme(legend.position="none") +
  theme(panel.grid.minor = element_line(linewidth = 0.1), panel.grid.major = element_line(linewidth = .25)) +
  theme(panel.spacing.y = unit(.35, "lines")) 
#agestructurehistogram

#Combine all three graphs and print.  This is a fairly complex graphic that takes most computers at least 10-20 seconds to print, and it’ll require a large graphics window (i.e, 10 inches tall or taller).
library(ggpubr)
combine_graphs <- ggarrange(agestructurehistogram, demo_sites_plot, rig_sites_plot, heights = c(1, 1.23, 6), labels = c("", ""), ncol = 1, nrow = 3, align = "v")
combine_graphs
