`belongs.to` <-
function(element, setsList) {
    sapply(setsList, has.element, element)
}

