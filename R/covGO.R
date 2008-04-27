`covGO` <-
function(pn, simplify=T) {
#
# Description:  Generic function to compute the covariance matrix associated to
#               sample (contracted) GO profile.
# 
# Usage:        covGO(pn, simplify=T)
#     
# Arguments:
# pn:       an object of class ExpandedGOProfile or an object that may be
#           coerced to this class
# simplify: should the result be simplified, if possible? See the 'Value' section
#
# Details:
#
#
# Value:    A list of covariance matrices. The length of the list equals the
#           number of GO profiles represented in the 'pn' argument. When this
#           length is one and simplifiy==T, the results consists in a single
#           matrix object, not a list of matrix objects.
#
# References: citar papers en construcció
#
# See Also:
#
# Examples:
#
#
    UseMethod("covGO")
}

