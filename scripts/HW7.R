# HW 7 Script

library(tidyverse)
library(spatstat)

# Source the functions and then use them below
source('scripts/generate_clustered_pt_proc.R')

##### Q1: Simulate data #####

# Simulate some data using generate_clustered_pt_proc.R

ppp_obj <- generate_clustered_pt_proc()
plot_point_pattern(ppp_obj, title = 'Default settings')

# 1.1 Which parameter(s) control the spatial extent of clusters?
# Xmin, Xmax, Ymin, and Ymax (see code examples below)

# First try X range lower
generate_clustered_pt_proc(Xmin = -5, Xmax = 5) %>% 
  plot_point_pattern(title = 'X range: -5 to 5')
# Plot default in between to compare
plot_point_pattern(ppp_obj, title = 'X range: -50 to 50')
# Now X range higher
generate_clustered_pt_proc(Xmin = -500, Xmax = 500) %>% 
  plot_point_pattern(title = 'X range: -500 to 500')

# First try Y range lower
generate_clustered_pt_proc(Ymin = -5, Ymax = 5) %>% 
  plot_point_pattern(title = 'Y range: -5 to 5')
# Plot default in between to compare
plot_point_pattern(ppp_obj, title = 'Y range: -50 to 50')
# Now Y range higher
generate_clustered_pt_proc(Ymin = -500, Ymax = 500) %>% 
  plot_point_pattern(title = 'Y range: -500 to 500')

# 1.2 Which parameter(s) control the strength of clustering (i.e., density 
#   within clusters relative to outside of clusters)?

# `val.at.center` appears to control how tightly the points are 
# within a cluster

# First try `val.at.center` lower
generate_clustered_pt_proc(val.at.center = 0.1) %>% 
  plot_point_pattern(title = '`val.at.center` = 0.1')
# Plot default in between to compare
plot_point_pattern(ppp_obj, title = 'val.at.center` = 1')
# Now `val.at.center` higher
generate_clustered_pt_proc(val.at.center = 10) %>% 
  plot_point_pattern(title = '`val.at.center` = 10')

# `effect.range` also appears to do this, but has a nonlinear
# relationship because 1 is more spread than 10 and 100 is
# also more spread than 10.

# First try `effect.range` lower
generate_clustered_pt_proc(effect.range = 1) %>% 
  plot_point_pattern(title = '`effect.range` = 1')
# Plot default in between to compare
plot_point_pattern(ppp_obj, title = '`effect.range` = 10')
# Now `effect.range` higher
generate_clustered_pt_proc(effect.range = 100) %>% 
  plot_point_pattern(title = '`effect.range` = 100')

# Similar to `val.at.center`, `background` appears to control
# how tightly (or not) the points are packed within a cluster.

# First try `background` lower
generate_clustered_pt_proc(background = 0.0005) %>% 
  plot_point_pattern(title = '`background` = 0.0005')
# Plot default in between to compare
plot_point_pattern(ppp_obj, title = '`background` = 0.001')
# Now `background` higher
generate_clustered_pt_proc(background = 0.1) %>% 
  plot_point_pattern(title = '`background` = 0.1')
generate_clustered_pt_proc(background = 1) %>% 
  plot_point_pattern(title = '`background` = 1')

# 1.3 Generate point pattern data from a complete spatial randomness (CSR) 
#   process and a clustered process and paste the two plots below.

# Read that CSR is really a Poisson here: https://inside.mines.edu/~jdzimmer/tutorials/Section2.html

csr_ppp <- generate_csr_pt_proc(seed=19)
clust_ppp <- generate_clustered_pt_proc(seed=19)

plot_point_pattern(csr_ppp, random = TRUE, title = 'Complete spatial randomness')
plot_point_pattern(clust_ppp, title = 'Clustered')

##### Q2: Test significance from CSR #####

# Use the quadrat test to determine whether each of these plots differs 
# significantly from CSR. You can either code this yourself or, if that seems 
# daunting, use the quadrat.test() function in the spatstat library.  Report the 
# Chi-square statistic and p value for each plot above.

# ?quadrat.test
csr_quad_test <- quadrat.test(csr_ppp)
clust_quad_test <- quadrat.test(clust_ppp)

csr_quad_test$statistic
clust_quad_test$statistic

csr_quad_test$p.value
clust_quad_test$p.value

##### Q3: Ripley's K plot #####

# Describe the degree of clustering at different spatial scales using a Ripley's 
# K plot.  Either code it yourself using eq. 2.8 from Fortin and Dale or use the 
# Kest() function in the spatstat library and the envelope() function to 
# generate an envelope for the null expectation for K for CSR data. Paste the 
# plot below.

# Calculate Ripley's K
clust_ripK <- Kest(clust_ppp, correction = "Ripley")
plot(clust_ripK)

# Add the envelope (confidence interval)
clust_env <- envelope(clust_ppp, Kest)
plot(clust_env, add=TRUE)

# In my simulated clusters, they are all statistically spatially clustered
# no matter the spatial scale (all fall above the CSR K line). See this resource:
# https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/GUID-110520A9-402D-4C17-8486-A7EC0F827D83-web.png

##### Q4: Generate spatial point process data #####

# Can you generate spatial point process data that are clustered at smaller 
# spatial scales but random at larger scales?  Paste a plot of the spatial point 
# pattern and a plot of Ripley's K below.

test_ripleysK <- function(vac = 1, er = 10, title = "Ripley's K for default") {
  csr_ppp2 <- generate_csr_pt_proc(seed=19)
  clust_ppp2 <- generate_clustered_pt_proc(val.at.center = vac, effect.range = er, seed=19)
  
  plot_point_pattern(csr_ppp2, random = TRUE, title = 'Complete spatial randomness')
  plot_point_pattern(clust_ppp2, title = 'Clustered')
  
  # Calculate & plot Ripley's K
  clust_ripK2 <- Kest(clust_ppp2, correction = "Ripley")
  clust_env2 <- envelope(clust_ppp2, Kest)
  
  plot(clust_ripK2, main = title)
  plot(clust_env2, add=TRUE)
}

# Look for pattern where the simulated data line is above the CSR line & 
# envelope at low values of `r` (left side of x-axis), and below that line & 
# envelope at high values of `r` (right side of x-axis)

# Change `effect.range`? No, that just seems to make it converge to the CSR
test_ripleysK() # This is the default for comparison
test_ripleysK(er = 1, title = "effect.range = 1") 
test_ripleysK(er = 100, title = "effect.range = 100") 

# Change `val.at.center`? # Also converging on CSR
test_ripleysK() # This is the default for comparison
test_ripleysK(vac = 0.1, title = "val.at.center = 0.1") 
test_ripleysK(vac = 10, title = "val.at.center = 10") 
test_ripleysK(vac = 0.01, title = "val.at.center = 0.01") 
test_ripleysK(vac = 0.001, title = "val.at.center = 0.001") 
test_ripleysK(vac = 0.0001, title = "val.at.center = 0.0001") 

# Is the problem that the data being generated is "thinned"?
