##############################################################################################################
# RURO-STYLE INFO-PLOTS FROM DATA : by Mahima Agarwal (https://www.linkedin.com/in/mahima-agarwal-b9856562/ ) and Martin Hilbert @UCDavis.edu  #
##############################################################################################################
#  This code is a R implementation to produce the information theoretic measures and plots analyzed in:
    ##  Crutchfield & Feldman, (2003). Regularities Unseen, Randomness Observed ["RURO"]: Levels of entropy convergence. Chaos: An Interdisciplinary Journal of Nonlinear Science, 13(1), 25-54.
#  The code allows to make RURO-style plots from time-series data  #
#  The logic of the calculations is adapted from the Python library https://github.com/dit/dit ; Python package for information theory. http://docs.dit.io, 
#  We are thankful for the help and guidance of dit's author, Ryan James (https://github.com/Autoplectic)

# You can run this complete script by running source
# in the R console (Remember to change the working directory!)
##          OR          ##
# Run line by line in R-studio by pressing CTRL+ENTER to execute a line

###############################################################################
# Imports

# Note: You may need to install these packages first.
#   Do this using the install.packages('...') command 

library(zoo)
library(grid)
library(gridBase)
library(pander)

###############################################################################


###############################################################################
# Counts and Entropy

# This function calculates counts of elements in rolling windows of length l in 
# a time series.
# Inputs: (1) ts: Time Series vector
#         (2) l: Length of rolling window
# Output: (1): counts: Dataframe with 2 columns: keys: Unique combinations, 
#                                                counts: counts in the rolling window application
counts.from.data<- function(ts, l)
{
  ts.rolled <- as.data.frame(rollapplyr(ts, l, c))         # Return elements of rolling windows
  ts.rolled.string<- do.call(paste, c(ts.rolled, sep="-"))  # Join the elements of the rolling window into a string
  # Count number of times each combination occurs using aggregate
  counts = aggregate(rep(1, length(ts.rolled.string)), by=list(ts.rolled.string), sum)   
  colnames(counts)<- c('keys', 'values')
  return(counts)
}


# Entropy estimator
# This function calculates entropy from counts data
# Input: (1) counts: Dataframe with counts data. Must have a values column.
# Output: (1) entropy: The calculated entropy value
H.0<- function(counts)
{
  N = sum(counts$values)  # Total number of elements = sum of all counts
  k = counts$values   # Vector of values
  # Calculate entropy by applying (ki/N * log2(ki/N)) to each element in k, and then 
  # taking a sum over all these values
  entropy = -sum(k/N * log2(k/N))
  return(entropy)
}
##############################################################################


###############################################################################
# Plotting functions

# Create Plots
# Given a time series, this function plots RURO style plots for the data.
# Inputs: (1) ts: Time series vector
#         (2) max.l: Maximum word length (Plots will be for 0: max.l)
#         (3) counts.pref (Default: 60): For calculating the 'worst case average word samples'
#               This aims at having some x number of each of the words that are created by the sliding window.
#               It follows the rule: WordLength = LOG2(((sampleSize) / x) - alphabet ).
#         (4) entropy.func (Default: H.0): Function to calculate the chosen entropy metric 
#         (5) digits (Default=6): Precision to be used in the table
#         (6) legend (Default='ON'): Include label inserts in plots ('ON'/'OFF')
make.ruro.plots<- function(ts, max.l, counts.pref=60, entropy.func=H.0, digits=6, legend='ON')
{
  A = length(unique(ts))    # Length of unique elements in A: alphabet
  ls.l <- 0:max.l         # Construct a range from 0 to max.l for word lengths to use
  bes <- rep(0, max.l+1)   # Initialize a vector of length ls.l with 0s
  # Calculate entropy for each values of l
  for (l in ls.l)       
  {
    # Replace bes[l+1] with the calculated entropy for counts with l elements in sliding window
    bes[l+1] = entropy.func(suppressWarnings(counts.from.data(ts, l)))
  }
  hmuL = c(NA, diff(bes))      # Calculate lagged differences from time series
  hmu = hmuL[length(hmuL)]     # Last value in hmuL
  
  Es = c(NA, diff(bes) - hmu)     # Subtract hmu from each lagged difference
  E = bes[length(bes)] - max.l * hmu
  
  Ts = (E + ls.l*hmu - bes)
  T = sum(Ts)
  
  Gs = c(rep(NA, 2), -diff(bes, differences = 2))     # Calculate 2nd order lagged differences
  G = sum(Gs, na.rm = TRUE)         # Sum while ignoring NAs
  
  reliable.l = log2(length(ts)/counts.pref)/log2(A)     # worst case average word samples
  atol = 1e-2          # this is real data here, people!
  rtol = 1e-2
  # Check the number of values equal to 0, within tolerance measures
  R = sum(1 - (abs(E + hmu*ls.l - bes - 0) <= (atol + rtol * abs(0))))
  if (R >= max.l)
  {
    R = NA
  }

  # Make tables:
  # Format all numbers in the table with n digits
  output.table = cbind(ls.l, formatC(Ts, format='f', digits=digits, flag='0') ,
                       formatC(Es, format='f', digits=digits, flag='0'),
                       formatC(Gs, format='f', digits=digits, flag='0'),
                       formatC(hmuL, format='f', digits=digits, flag='0'))

  colnames(output.table) = c('word length', 'T', 'E', 'G', 'hmuL')
  # Add totals row
  output.table = rbind(output.table, c('Totals', formatC(sum(Ts, na.rm = TRUE), format='f', digits=digits, flag='0'),
                                       formatC(sum(Es, na.rm = TRUE), format='f', digits=digits, flag='0'),
                                       formatC(sum(Gs, na.rm = TRUE), format='f', digits=digits, flag='0'),
                                       c('-')))
  # Print formatted table to console
  pandoc.table(output.table, style='grid', plain.ascii=TRUE)
  
  #Create plots
  # Set plot parameters
  png('RURO_Style_TS_Plots.png', width=3600, height=3600, res=300)  # Plot will be created in the 
                                                                              # working directory
  par(mfrow = c(3,1))  # Plot 3 subplots in single plot, by row
  
  # Add chart 1: H vs L
  par(mar=c(1, 4, 3, 2) + 0.1)  # Set chart margins
  # Plot first line in chart, without any axes
  plot(x=ls.l, y=(E + ls.l*hmu), type='l', col = 'blue', lty=2, xaxs="i", xlab='',
       ylab = '', axes=FALSE)
  # Add second line
  lines(x=ls.l, y=bes, col = 'blue', lty=1, lwd=1.5)
  # Fill area between the two lines
  polygon(c(ls.l, rev(ls.l)), c(bes, rev(E + ls.l*hmu)),
          col = rgb(0,0,1,0.5), border = NA)
  plot.area = par('usr')          # Get plot area coordinates
  rect(plot.area[1], plot.area[3], plot.area[2], plot.area[4], border = 'black', lwd=2)    # Box around plot area
  ell = '\u2113'  # Unicode for l
  # Add y-axis label with formatted expression
  title(ylab=bquote('Entropy' ~ italic(H)~italic((.(ell)))), line=1.7, cex.lab=2)
  
  if (legend == 'ON')
  {
    legend.position = c(0.5, round(max(bes, 2)))     # Can modify legend position for better 
    # suited coordinates 
    legend(legend.position[1], legend.position[2], col = c('blue', 'blue', rgb(0, 0, 1, 0.5)), cex=2,
           legend=c(
             as.expression(
               bquote(bold(E) + .(ell) * h[mu] ~ '=' ~.(round(E, 4)) + .(ell) %.% .(round(hmu, 4)))
             ),
             as.expression(bquote(italic(H)~italic((.(ell))))),
             as.expression(bquote(italic(T) ~ '=' ~ .(round(T, 4))))),
           lty=c(2, 1, NA), merge=T, pch = c(NA, NA, 15), pt.cex = c(1,1,5)) 
  }
  axis(1, at=ls.l, labels=FALSE)     # Add x-axis without tick-labels
  axis(1, at=c(0), cex=1.5)          # Add tick-label for x=0
  axis(2, cex=1.5)                   # Add y-axis
  abline(v=reliable.l, lty=3)        # Add vertical line of "worst case average word samples"
  
  # Add chart 2: hmu vs L
  par(mar=c(1, 4, 1, 2) + 0.1)   # keep same figure margins, but reduce space at top
  plot(x=c(0, ls.l[-1]), c(NA, diff(bes)), 'l', col = 'blue', lty=1, xaxs="i", xlab='', ylab = '',
      xaxt='n', lwd=1.5, yaxt='n', axes=FALSE)             # Plot first line, but ensure that x-axis starts from 0
                                                            # This is important to ensure that the charts line up
  lines(ls.l[-1], y =rep(hmu, max.l), col='blue', lty=2)          # Add second line
  polygon(c(ls.l[-1], rev(ls.l[-1])), c(rep(hmu,max.l), rev(diff(bes))), 
          col = rgb(0, 1, 1, 0.5), border = NA)        # Fill area between lines
  plot.area = par('usr')     # Get plot coordinates
  rect(1,plot.area[3], plot.area[2], plot.area[4], border='black', lwd=2)  # Draw box starting from x=1 instead of x=0

  if (legend == 'ON')
  {
    legend.position = c(max.l - max.l/4, max(diff(bes)))    # Can modify legend position for better 
                                                               # suited coordinates 
    legend(legend.position[1], legend.position[2] , col = c('blue', 'blue', rgb(0, 1, 1, 0.5)), cex= 2,
           legend=c(
             as.expression(bquote(italic(h[mu]) ~ ~italic((.(ell))))),
             as.expression(bquote(italic(h[mu]))),
             as.expression(bquote(bold(E) ~ '=' ~ .(round(E, 4))))),
           lty=c(1, 2, NA), merge=T, pch = c(NA, NA, 15), pt.cex = c(1,1,5))
  }
  axis(1, at=ls.l[-1], labels=FALSE)     # Add x-axis starting from x=1
  axis(1, at=c(1), cex=1.5)           # Add x=1 label
  axis(2, pos=1, cex=1.5)        # Add y-axis at x=1
  abline(v=reliable.l, lty=3)        # Add vertical line of "worst case average word samples"
  
  # Identify location for y-axis title:
  x.1 = grconvertX(1, 'user', 'inches')   # Convert x=1 chart coordinate into inches 
  # Identify x and y-coordinates for y-axis labels
  x.label = grconvertX(x.1-0.37, 'inches', 'user')  # 0.37 inches before x=1, converted to plot coordinates
  y.label = (plot.area[4] + plot.area[3]) / 2       # y-axis midpoint in plot coordinates
  # Add y-axis title
  text(x.label, y.label, bquote('Entropy rate' ~ italic(h[mu])~italic((.(ell)))), cex = 2, srt=90)
  
  # Add chart 3: Delta^2 H(l)
  par(mar=c(4, 4, 1, 2) + 0.1)    # Use same margins as chart 2 except more space at the bottom, for x-axis label
  plot(x=ls.l, c(NA, NA, diff(bes, difference=2)), 'l', col = 'blue', lwd=1.5,  axes = FALSE,     # Plot first line, but ensure that x-axis starts from 0
       xaxs="i", xlab='', ylab = '', xaxt='n', ylim = c(min(diff(bes, difference=2)), 0))        # This is important to ensure that the charts line up
  lines(x=c(2, max.l), y= c(0, 0), col='blue', lty=2)   # Add second line
  # Fill area between lines
  polygon(c(ls.l[-c(1,2)], rev(ls.l[-c(1,2)])), c(diff(bes, difference=2), rev(rep(0,length(ls.l) - 2))),
          col = rgb(1, 0, 1, 0.5), border = NA)
  plot.area = par('usr')  # Get plot coordinates
  rect(2,plot.area[3], plot.area[2], plot.area[4], border='black', lwd=2)    # Draw box starting from x=2
  # Add x-axis title
  title(xlab=bquote('Word length' ~ italic(.(ell))), line=3, cex.lab=2, adj= ((max.l-2)/2 +2)/max.l)
  if (legend == 'ON')
  {
    legend.position = c(max.l - max.l/4, (par('usr')[3])/2)   # Can modify legend position for better 
    # suited coordinates 
    legend(legend.position[1], legend.position[2], col = c('blue', rgb(1, 0, 1, 0.5)), cex=2,
           legend=c(
             as.expression(bquote(Delta^2 ~ italic(H) ~ italic((.(ell))))),
             as.expression(bquote(G ~ '=' ~ .(round(G, 4))))),
           lty=c(1, NA), merge=T, pch = c(NA, 15), pt.cex = c(1, 5)) 
  }
  axis(1, at=c(2:max.l), cex=1.5)  # Add x-axis with tick labels, string from x=2
  axis(2, pos=2, cex=1.5)    # Add y-axis
  abline(v=reliable.l, lty=3)      # Add vertical line of "worst case average word samples"
  # Identify location for y-axis title:
  x.2 = grconvertX(2, 'user', 'inches')   # Convert x=2 to inches
  # Get title x and y coordinates
  x.label = grconvertX(x.2-0.37, 'inches', 'user')    # x=2 - 0.37 inches to plot coordinates
  y.label = (plot.area[4] + plot.area[3]) / 2        # midpoint of y-axis
  # Add y-axis title expression
  text(x.label, y.label, bquote(Delta^2 ~ italic(H) ~ italic((.(ell)))), cex = 2, srt=90)
  dev.off() # Close plotting device
}

###############################################################################
# Import time series
# set working directory 
setwd("C:/Users/hilbert/...")

# load csv file 
# in ('...') name of csv-file; in $'..' name of column (with special characters replaced with .)
GitHub.ts = read.csv('OSM_Russia.csv')$A.37
###############################################################################


###############################################################################
# Testing
# And now we can actually run the code with the simple line below!
# Adjust the first number to set the 'maximal word length' to assess the measures (aka number of consecutive letter in the time sequence)
# Adjust the second number (count.pref) to create a line in the graph that shows at what point each world should (probablisitcally) appear x number of times (to see sampling reliability)
# Adjust the third number to specify how many digits you want to be printed for the calculated information measures

make.ruro.plots(GitHub.ts, 6, counts.pref=10, entropy.func=H.0, digits=4, legend='ON')

# Find the output plots in your working directory at ./RURO_Style_TS_Plots.png