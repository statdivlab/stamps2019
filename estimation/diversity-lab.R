## diversity-lab.R
## A script to introduce you to diversity estimation and comparison

## Before starting this tutorial, please run
devtools::install_github("adw96/breakaway")
library(breakaway)
# Let's test whether we can use the apples dataset in breakaway
# to confirm everything is loading.
data(apples) 
head(apples) # Something should output
devtools::install_github("adw96/DivNet")
library(DivNet)
divnet ## something should output

# If this doesn't work, run install.packages("devtools") and try again
# If this doesn't work, read the error message carefully 
# and try to debug it yourself (this is what you have to do at home -- so
# it's good to practice here)
# If this doesn't work, call a TA over

## Lab authors: Amy Willis, Bryan Martin, Pauline Trinh

## breakaway authors: Amy Willis, Kathryn Barger, John Bunge, 
##                        Bryan Martin, 2012+

## Start by checking what directory you are in with
getwd()

# You will need to change your working directory to what you want to work in,
# probably your STAMPS folder!
## Set it to your local copy: setwd([your directory name here])

library(tidyverse)

## Let's use data distributed from phyloseq
library(phyloseq)
data("GlobalPatterns")
GlobalPatterns 
## We see that GlobalPatterns is a phyloseq object containing an 
## otu table, sample data, taxonomy table, and a phylogenetic tree

## Have a look through GlobalPatterns samples
GlobalPatterns %>% sample_data


# Let's just look at the water samples of GlobalPatterns
# To speed things up, let's just aggregate taxa to the 
# order level. So we're going to estimate *order* level 
# diversity 

# This will take a while, but it will speed things up later
water <- GlobalPatterns %>%
  subset_samples(SampleType %in% c("Freshwater", 
                                   "Freshwater (creek)", 
                                   "Ocean", 
                                   "Sediment (estuary)")) %>%
  tax_glom("Order")

water

# phyloseq has some built in tools for exploring alpha
# diversity, but they're not great. They underestimate richness,
# underestimate uncertainty, and don't allow hypothesis testing

# Enter: breakaway :D
# breakaway was specifically designed for these tasks

## Let's look at the observed richness of the water samples
observed_c <- sample_richness(water)
summary(observed_c)
plot(observed_c, water, color = "SampleType")

# Hmmmm, but what if we observed these samples at different depth?
# Depth may be confounded with observed richness. 
# Let's practice ggplot and phyloseq to look at this

data.frame("observed_richness" = (observed_c %>% summary)$estimate,
           "depth" = phyloseq::sample_sums(water), # Easter egg! Phyloseq's function to get depth
           "type" = water %>% sample_data %>% get_variable("SampleType")) %>%
  ggplot(aes(x = depth, y = observed_richness, color = type)) +
  geom_point()

# So our lowest depth samples (Sediment) had lower richness
# compared to Freshwater (creek)!
# Not surprising!

# Recall that the best way to adjust for this is by estimating the 
# number of missing species using a species richness estimate.

ba <- breakaway(water)
ba 
plot(ba, water, color = "SampleType")

## Cool! How do these estimates work?

## Let's look at just one sample
tr <- water %>% subset_samples(X.SampleID == "TRRsed1")
tr

# what's the structure of this dataset?
fc <- tr %>% otu_table %>% make_frequency_count_table
# this is the frequency count table for this dataset
fc %>% head(10)
# So there are 11 singletons (Orders observed once) here

# Quiz: how many times was the most common Order observed?
# (Hint: what does tail() do?)
# (A: 8863 times)

# Let's fit the breakaway to this sample
ba_tr <- breakaway(fc)
ba_tr
# this is an alpha diversity estimate -- a special class for
# alpha diversity estimates

# breakaway picks lots of models and chooses the best
# Which model did it pick?
ba_tr$model

# Kemp -- that's a non-mixed Poisson model. Cool!

# Kemp models work by fitting a legit probabilistic model to
# a transformation of the data. We can plot the transformation and the fit:
ba_tr %>% plot
# (Not super important -- but if you want to know more, check out
# the paper: Willis & Bunge (2015), Biometrics)

# breakaway is easy to use! You can run it on one frequency table, but
# in general you'll run it on a phyloseq object, like here:
ba <- breakaway(water)
ba 

# you can plot using our default plotting function...
plot(ba, water, color = "SampleType")

# ... or you can take the estimates and turn them into 
# data frames so you can plot them yourself
# these are in the same order as the phyloseq samples:
summary(ba) %>%
  add_column("SampleNames" = water %>% otu_table %>% sample_names)


# Also, since the breakaway package implements lots of
# species richness estimates, you could choose a different one
# e.g., if you wanted something more stable
water %>%
  chao_bunge %>%
  plot(water, color = "SampleType")

# However, note that species richness is a challenging problem, 
# and that error bars will generally be large

# grrrrr don't use this one
water %>%
  chao1 %>%
  plot(water, color = "SampleType")
# (but note that these are the real error bars on this estimate)

# In many cases the absolute number of species
# isn't as interesting as comparing ecosystems. 
# Let's test the hypothesis that different types of water systems
# have the same microbial diversity

# betta() works like a regression model
# but it accounts for the uncertainty in estimating
# diversity


bt <- betta(summary(ba)$estimate,
            summary(ba)$error,
            make_design_matrix(water, "SampleType"))
bt$table

# betta() estimates that the mean Order-level
# diversity in Freshwater is 154 orders.
# It estimates that the diversity in
# creeks is significantly higher
# (on average 23 orders) while oceans
# have significantly lower diversity
# (on average, 16 orders). However, 
# estuaries do not have significantly different diversity 
# than fresh water sources.

# Note that these estimates account for
# different sequencing depths!
# breakaway estimates the number of missing species
# based on the sequence depth and
# number of rare taxa in the data

# The citation for betta() is
# Willis, Bunge & Whitman, 2016, JRSS-C

# Please use it! It's very important that
# you account for the error bars in
# diversity when doing hypothesis testing.
# t.test, lm(), aov(), anova() do not
# account for this -- betta does!

#### #### #### #### #### #### #### #### 
#### OTHER DIVERSITY INDICES
#### #### #### #### #### #### #### #### 

# Species richness counts all species equally.
# However, if a species is rare, you may think that
# it doesn't play a role in the community.
#  Another alpha diversity index, called the Shannon
# index, works similarly to species richness but it
# downweights the importance of rare taxa
# i.e. if a taxon is present but in low abundance, 
# it doesn't count for "1" taxon, but something less 
# (to reflect its low abundance)

# Since rare taxa may be dubious, the Shannon index
# is very popular in microbial ecology

# For the reasons discussed in the lecture, it's important to
# estimate the Shannon diversity using a valid
# estimator of the Shannon diversity.
# It's even more important to come up with
# the standard error, and use that
# standard error in testing using betta()

# Shockingly, until recently, there were no tools
# to estimate Shannon diversity
# in the presence of an ecological/microbial network!

# DivNet is a new tool that allows you to estimate this.
# It also adjusts for different sequencing depths
# so you don't have to throw away any data
# (you don't need to rarefy)!

# Check out our preprint: Willis & Martin (2018+), bioRxiv

# Let's load the package DivNet
library(DivNet)

# to test if DivNet loaded correctly, try
dv_water_testing <- divnet(water, ncores = 1, tuning = "test")
# If this command loads but other commands
# don't load, the problem is not DivNet,
# but a package that DivNet depends on

# DivNet is very flexible, but by default
# it estimates the
# microbial network and uses it only to
# adjust the standard errors on diversity estimates

# *If you looked at the package parallel*
# in the "more R" tutorials, 
# you should be able to run DivNet in parallel: 
dv_water <- divnet(water, ncores = 4)
# This might take a minute; if it errors,
# try running it in series (slower
# but more portable):
# dv_water <- divnet(water, ncores = 1)

# DivNet outputs a list of the estimates
# shannon, simpson (alpha diversity)
# bray-curtis, euclidean (beta diversity)
dv_water %>% names

# You can pull them out individually:
dv_water$shannon %>%
  summary %>%
  add_column("SampleNames" = water %>% otu_table %>% sample_names)

# or plot them:
plot(dv_water$shannon, 
     water, 
     col = "SampleType")

# Let's compare this to the naive approach of
# just "plugging in" the observed proportions
# to the Shannon diversity  formula

plot(water %>% sample_shannon, 
     water, 
     col = "SampleType") + ylim(0, 3.5)

# You will notice that the estimates are almost the same
# as previously, but the error bars differ significantly
# i.e., *there are error bars*

# The error bars were the same as previously
# because we haven't told DivNet anything about the experimental
# design. Here we observed
# 4 different water systems, so we will add this
# as a covariate
dv_water_st <- water %>%
  divnet(X = "SampleType", ncores = 8)

plot(dv_water_st$shannon, 
     water, 
     col = "SampleType")
# Now we see that organisms in the same group 
# are estimated to have the same diversity.
# We are now estimating the diversity 
# of this *type* of ecosystem. We're
# focusing on the ecosystem, not
# just the samples. If we want to 
# reproduce the results of our study,
# it's better to focus on the populations that
# the samples come from, not the samples themselves

# You can analyse sex, time, disease status in this way.
# It's actually ideal for longitudinal/epigenetic
# studies, because you can say something about the
# groups, not just the people/mice in your study

# Let's look at hypothesis testing for DivNet

testDiversity(dv_water_st, "shannon")
# So we have significantly lower Shannon diversity
# in creeks, oceans and estuaries than in
# freshwater. 

# This is just a wrapper for betta(), btw

# We can also do the same thing for Simpson, of course:
plot(dv_water_st$simpson, 
     water, 
     col = "SampleType")
testDiversity(dv_water_st, "simpson")

## To test hypotheses about beta diversity
# using DivNet, let's pull out our
# estimated distance matrix:
bc <- dv_water_st$`bray-curtis`

# You'll notice that all samples with the same
# SampleType at the same estimate. DivNet uses
# covariate information to share strength
# across samples and obtain an estimate
# about the beta diversity of the *ecosystems*
#   not the samples

# We can consider the  unique rows using
bc %>% unique

# Uniquely, DivNet also has variance estimates:
simplifyBeta(dv_water_st, water, "bray-curtis", "SampleType")
simplifyBeta(dv_water_st, water, "euclidean", "SampleType")

# You can plot this easily
simplifyBeta(dv_water_st, water, "bray-curtis", "SampleType") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est)) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")

# Finally, a teaser! Soon we will have
# support for phylodivnet(), to estimate UniFrac
# in the presence of a network. Subscribe
# to our github/Amy's Twitter feed
# for updates. 


## THANKS FOR STAYING THROUGH TO THE END! 
#
# I'm always happy to take questions
# either today, later, or after the course.
# 
# I'd also  love to hear if you
# hated it/loved it/worst day ever -- 
# please give me feedback to 
# better help you!! 
# 
# I want these tools to be used so if you
# have feedback on how to make them
# easier to use, please let me know.
#
# :) Amy
