* output.do
set more off
capture log close

sjlog using output1, replace
sysuse auto
regress mpg weight
sjlog close, replace

sort weight
predict yhat
set scheme sj
scatter mpg yhat weight, c(. l) s(x i)
graph export output1.eps, replace

exit
