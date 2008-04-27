`vecProfile` <-
function (exp.Prof){
    vecProf <- exp.Prof$Freq[-length(exp.Prof$Freq)]
    names(vecProf)<-exp.Prof$res[-length(exp.Prof$res)]
    attr(vecProf, "n")<-sum(vecProf)
    vecProf <-vecProf/attr(vecProf, "n")
}

