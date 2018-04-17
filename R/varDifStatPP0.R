`varDifStatPP0` <-
function(pn, p0 = NULL, funcProfP0 = NULL, simplify=T) {
#
# Description:  Generic function to compute the variance of the fit statistic
#               sqrt(n)d(Pn,P0) (or equivalently the variance of 
#               sqrt(n){d(Pn,P0) - d(P,P0)}), where:
#               d is the square Euclidean distance, Pn is a sample GO functional
#               profile associated to the expanded profile 'pn' (based
#               on a sample of n genes) and P0 corresponds to a "model" profile
#               represented either by 'funcProfP0' or by 'p0' (then, in an
#               "expanded" profile form)
# 
# Usage:        varDifStatPP0(pn, p0 = NULL, funcProfP0 = NULL, simplify=T)
#     
# Arguments:
# pn: 
# p0:
# funcProfP0:      
# simplify: should the result be simplified, if possible? See the 'Value' section
#
# Details:
#
#
# Value:    A list or a vector of variance values. Its length equals
#           max(ncol(pn),ncol(p0),ncol(funcProfP0)). When simplifiy==T, the
#           result is a vector, otherwise it is a list.
#
# References: citar papers en construccio
#
# See Also:
#
# Examples:
#
    UseMethod("varDifStatPP0")
}

