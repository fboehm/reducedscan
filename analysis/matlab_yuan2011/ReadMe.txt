1) Description:
f.m is a matlab funtion for generating MCMC samples. 


Arguments in function f.m: 
datafile    a data matrix,  containing the variables in the model, including the parent information for RIX individuals.

             data format:
             1st column: sample ID
             2nd column: trait value
             3rd and 4th columns: parental RI lines
             rest of columns: genotype data

drop    Number of steps dropped from the chain to allow for a burn-in phase 

delta   An optional delta value used in the prior distributions of variance parameters, 
        as described in ter.Braak (2005).

thin    Thinning facor. 

hmax    Number of MCMC steps, i.e. length of the Markov chain. 


Value

a) a 2*num.marker matrix of posterior means of all markers. The odd rows are for additive effects, and the even rows 
are for dominance effects. 
b) two figures, one for the additive effect and another one for the dominance effect.
In each figure, x-axis refers the marker positions (if marker positions are given) or  marker numbers (if marker 
positions are not given); and y-axis refers the posterior means.







2) example matlab codes:
%% load data file, which is datamarker50.txt in this example
load datamarker50.txt 

%% call function f.m
f(datamarker50)

