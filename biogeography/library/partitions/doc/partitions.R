### R code from vignette source 'partitions.Rnw'

###################################################
### code chunk number 1: partitions.Rnw:122-122
###################################################



###################################################
### code chunk number 2: partitions.Rnw:123-124
###################################################
require(partitions)


###################################################
### code chunk number 3: partitions.Rnw:132-133
###################################################
P(100)


###################################################
### code chunk number 4: partitions.Rnw:137-138
###################################################
diffparts(10)


###################################################
### code chunk number 5: partitions.Rnw:170-176
###################################################
f <- function(x){prod(factorial(x),factorial(tabulate(x)))}
prob <- function(a,n){
  jj <- restrictedparts(n,a,include.zero=FALSE)
  N <- factorial(a)*factorial(n)*sum(1/apply(jj,2,f))
  return(N/a^n)
}


