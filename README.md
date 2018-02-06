# AttributedSBM
#Objective: This function is for fitting the attributed stochastic block model to a network with multiple continuous attributes.
#Last Edited: February 6, 2017
#For bugs: Please contact NatalieStanley1318@gmail.com
#Implementation details: Tested in R version 3.4.3
#Dependencies:igraph and mvtnorm

#Inputs:
  #Network: An NxN adjaceny matrix. It can be sparse
  #Attribute Mat: NxP vector of attributes, where p is the number of attributes you have
 #Prob: An indicator for wheether your attributes represent a probability of being in each of p communities. Use 0 if your     attributes are not probabilities. 

#Outputs: A list object with entries
  $Comm: Node to community assignment
  $SBMProb: SBM probability parameters (between community probability matrix)
  $Mean: The mean of the Gaussian describing each community. This is a list object with each entry corresponding to the community of the index.
  $Cov: The covariance matrix describing each community. This is a list object with each entry corresponding to the community of the index.
  
  
  Example:
  
  library('igraph')
  library('mvtnorm')
  source('FitAttribute.R')
  Out=FitAttribute(MyNet,MyAttribute,0)
