capture log close
sjlog using ex0, replace

use example.dta

describe 

sjlog close, replace


capture log close
sjlog using ex1, replace

teddf K L= Y: CO2, dmu( Province ) time(year) sav(ex1result,replace)

sjlog close, replace


capture log close
sjlog using ex2, replace

teddf K L= Y: CO2, dmu( Province ) time(year) nonr sav(ex2result,replace)

sjlog close, replace


capture log close
sjlog using ex3, replace

egen id=group(Province)
xtset id year
gtfpch K L= Y: CO2, dmu( Province ) global sav(ex3result,replace)

sjlog close, replace


capture log close
sjlog using ex4, replace

gtfpch K L= Y: CO2, dmu( Province ) nonr  global sav(ex4result,replace)

sjlog close, replace