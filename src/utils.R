# simple function
# if x is NOT null evaluate expression y otherwise do nothing
ifnotnull <- function(x, y){
  if(!is.null(x)){
    return(eval(quote(y)))
  }
  else(return(x))
}

# another simple function (used in verify below)
check_dims <- function (x, d, par) 
{
  if (is.vector(x) | is.factor(x)) {
    dc <- length(x)
  }
  else {
    dc <- unname(dim(x))
  }
  d <- unname(as.integer(d))
  if (!identical(dc, d)) {
    s1 <- paste(d, collapse = " x ")
    s2 <- paste(dc, collapse = " x")
    s <- paste("Parameter", par, "should have dimension", 
               s1, "instead found dimensions", s2)
    stop(s)
  }
}