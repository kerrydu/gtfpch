capture log close
sjlog using ex.describe, replace

use example.dta
describe 

sjlog close, replace


capture log close
sjlog using ex.teddf, replace

teddf K L= Y: CO2, dmu(Province) time(year) sav(ex.teddf.result,replace)

sjlog close, replace


capture log close
sjlog using ex.teddf.direction, replace

gen gK=0
gen gL=0
gen gY=Y
gen gCO2=-CO2
teddf K L= Y: CO2, dmu(Province) time(year) gx(gK gL) gy(gY) gb(gCO2) sav(ex.teddf.direction.result,replace)

sjlog close, replace


capture log close
sjlog using ex.teddf.nonr, replace

teddf K L= Y: CO2, dmu(Province) time(year) nonr sav(ex.teddf.nonr.result,replace)

sjlog close, replace


capture log close
sjlog using ex.teddf.nonr.weight, replace

mat wmatrix=(0.5,0.5,1,1)
teddf K L= Y: CO2, dmu(Province) time(year) nonr wmat(wmatrix) sav(ex.teddf.nonr.weight.result,replace)

sjlog close, replace


capture log close
sjlog using ex.gtfpch, replace

egen id=group(Province)
xtset id year
gtfpch K L= Y: CO2, dmu(Province) global sav(ex.gtfpch.result,replace)

sjlog close, replace


capture log close
sjlog using ex.gtfpch.nonr, replace

gtfpch K L= Y: CO2, dmu( Province ) nonr  global sav(ex.gtfpch.nonr.result,replace)

sjlog close, replace