
############################################################################
############################################################################
###                                                                      ###
###                           NITROGEN-STATS.R                           ###
###                                                                      ###
############################################################################
############################################################################

# This R script will analyse the nitrogen data generated from lab class week 2 and 5

##################################################################
##                            Set up                            ##
##################################################################

# Read in data files
germ_df <- read.csv("data/germination_data.csv", stringsAsFactors = TRUE)
length_df <- read.csv("data/length_data.csv", stringsAsFactors = TRUE)

# Load packages
library(car) # 'leveneTest' function
library(agricolae) # 'duncan.test' function


#################################################################
##                      Germination Data:                      ##
##                        One-way ANOVA                        ##
#################################################################

# Linear Model
germ_model <- lm(Germination~Treatment, data = germ_df)

# One-way ANOVA
germ_aov <- aov(germ_model)

# View p-values
summary(germ_aov)
#              Df Sum Sq Mean Sq F value   Pr(>F)    
# Treatment    9 1.3707  0.1523    15.7 8.36e-07 ***
# Residuals   18 0.1746  0.0097                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 1 observation deleted due to missingness


##################################################################
##                         Assumptions:                         ##
##                          Normality                           ##
##################################################################

# Plot normality of residuals; approx. linear; indicates normality
plot(germ_aov, 2)

# Extract the residuals
germ_residuals <- residuals(object = germ_aov)

# Run Shapiro-Wilk test; > 0.05 indicates normality
shapiro.test(x = germ_residuals)
# Shapiro-Wilk normality test
# 
# data:  germ_residuals
# W = 0.96502, p-value = 0.455

#################################################################
##                        Assumptions:                         ##
##              Equal variance / Homoscedasticity              ##
#################################################################


# Plot residuals vs Fitted; no wedge shape which is good; maybe homoscedastic.
plot(germ_aov, 1)

# Test for homogeniety of variance; > 0.05 indicates equal variance
leveneTest(germ_model, center = mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value Pr(>F)
# group  9  1.9205 0.1141
# 18     


#################################################################
##                      Germination Data:                      ##
##                      Post-hoc analysis                      ##
#################################################################


# Duncan's Multiple Range Test
germ_out <- duncan.test(germ_model,"Treatment")

# Visualise differences and plot significance
plot(germ_out, variation="SE",
     main = " ",
     xlab = "Treatment",
     ylab = "Germination Percentage (%)",
     las = 1,
     cex=1.5, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)

# > export > PDF > 18.75 x 10



#################################################################
#################################################################
#################################################################
#################################################################
##                          Root Data:                         ##
##                        One-way ANOVA                        ##
#################################################################

# Linear Model
len_model <- lm(Mean_Length~Treatment, data = length_df)

# One-way ANOVA
len_aov <- aov(len_model)

# View p-value
summary(len_aov)
#              Df Sum Sq Mean Sq F value  Pr(>F)   
# Treatment    9 248086   27565   5.275 0.00113 **
# Residuals   19  99291    5226                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##################################################################
##                         Assumptions:                         ##
##                          Normality                           ##
##################################################################

# Plot normality of residuals; approx. linear; indicates normality
plot(len_aov, 2)

# Extract the residuals
len_residuals <- residuals(object = len_aov)

# Run Shapiro-Wilk test; < 0.05 violates normality
shapiro.test(x = len_residuals)
# Shapiro-Wilk normality test
# 
# data:  len_residuals
# W = 0.80603, p-value = 0.0001073


#################################################################
##                        Assumptions:                         ##
##              Equal variance / Homoscedasticity              ##
#################################################################


# Plot residuals vs Fitted; no wedge shape which is good; maybe homoscedastic.
plot(len_aov, 1)

# Test for homogeniety of variance; < 0.05 violates equal variance
leveneTest(len_model, center = mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value    Pr(>F)    
# group  9  10.217 1.375e-05 ***
# 19                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



##################################################################
##                     Transform root data:                     ##
##                          Square root                         ##
##################################################################

# Square root transformations
length_df$Mean_Length_SqRoot ='^'(length_df$Mean_Length,1/2)

# Check new data values
head(length_df)


#################################################################
##                         Redo Stats:                         ##
##                        One-way ANOVA                        ##
#################################################################

# Linear Model
len_model_sqroot <- lm(Mean_Length_SqRoot~Treatment, data = length_df)

# One-way ANOVA
len_aov_sqroot <- aov(len_model_sqroot)

# View p-value
summary(len_aov_sqroot)
#              Df Sum Sq Mean Sq F value  Pr(>F)    
# Treatment    9  757.0   84.12   7.515 0.00012 ***
# Residuals   19  212.7   11.19                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#################################################################
##                      Redo Assumptions:                      ##
##                          Normality                          ##
#################################################################


# Plot normality of residuals; approx. linear; indicates normality
plot(len_aov_sqroot, 2)

# Extract the residuals
len_residuals_sqroot <- residuals(object = len_aov_sqroot)

# Run Shapiro-Wilk test; < 0.05 violates normality
shapiro.test(x = len_residuals_sqroot)
# Shapiro-Wilk normality test
# 
# data:  len_residuals_sqroot
# W = 0.83485, p-value = 0.000377


#################################################################
##                      Redo Assumptions:                      ##
##                        Equal variance                       ##
#################################################################



# Plot residuals vs Fitted; no wedge shape which is good; maybe homoscedastic.
plot(len_aov_sqroot, 1)

# Test for homogeniety of variance; < 0.05 violates equal variance
leveneTest(len_model_sqroot, center = mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value    Pr(>F)    
# group  9  9.1628 3.027e-05 ***
# 19                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


###########################################################################
###########################################################################
###                                                                     ###
###                     COMMENT ON TRANSFORMATIONS:                     ###
###                                                                     ###
###########################################################################
###########################################################################

# The square transformation did not resolve the violation of assumptions
# I will proceed with the untransformed root data for post-hoc analysis and visualisation
# The results will be intepreted with caution due to the violation of assumptions


#################################################################
##                          Root data:                         ##
##                      Post-hoc analysis                      ##
#################################################################


# Duncan's Multiple Range Test
len_out <- duncan.test(len_model,"Treatment")

# Visualise differences and plot significance
plot(len_out, variation="SE",
     main = " ",
     xlab = "Treatment",
     ylab = "Seedling root length (mm)",
     las = 1,
     cex=1.5, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)

# > export > PDF > 18.75 x 10
