* output.do
version 7

sjlog using output1, replace
use auto, clear
regress mpg weight
sjlog close, replace

mat b = e(b)
cap noi outtable using table1, mat(b) caption(Estimated Parameters) replace

sort weight
predict yhat
graph mpg yhat weight, c(.l) s(xi) xsize(3) ysize(2.5) sav(output1, replace)
translate output1.gph output1.eps, replace

exit
