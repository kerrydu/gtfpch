*! version 2.1, 7 Aug 2020
*! version 2.0
* By Kerry Du & Daoping Wang, 30 Jul 2020 
**
* 
capture program drop teddf
program define teddf, rclass
    version 16

    gettoken word 0 : 0, parse(" =:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
		local invars `invars' `word'
		gettoken word 0 : 0, parse("=:,")
    }
    unab invars : `invars'

	gettoken word 0 : 0, parse(" =:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
		local gopvars `gopvars' `word'
		gettoken word 0 : 0, parse(" =:,")
    }
    unab gopvars : `gopvars'
	
	
    syntax varlist [if] [in], Dmu(varname) [Time(varname) gx(varlist) gy(varlist) gb(varlist)  ///
	                                        BIennial  SEQuential GLObal VRS  NONRadial  TONE ///
										   Wmat(string) SAVing(string)  WINdow(numlist intege max=1 >=1)    ///
										   maxiter(numlist integer >0 max=1) tol(numlist max=1)]
    local bopvars `varlist'
	if `"`tone'"'!=""{
		sbmeff `invars' = `gopvars' : `bopvars', dmu(`dmu') time(`time')  `vrs' ///
		     `biennial' `sequential' `global'  window(`window') saving(`saving') maxiter(`maxiter') tol(`tol')
		exit
	}

										   
	marksample touse 
	markout `touse' `invars' `gopvars' `gx' `gy' `gb'
	
	local invars: list uniq invars
	local gopvars: list uniq gopvars
	local bopvars: list uniq bopvars
	
	local ninp: word count `invars'
    local ngo: word count `gopvars'
    local nbo: word count `bopvars'
	local nvar=`ninp'+`ngo'+`nbo'
	
	confirm numeric var `invars' `gopvars' `bopvars'

	local comvars: list invars & gopvars 
	if !(`"`comvars'"'==""){
		disp as error "`comvars' should not be specified as input and desriable output simultaneously."
		error 498
	}
	
	local comvars: list invars & bopvars
	if !(`"`comvars'"'==""){
		disp as error "`comvars' should not be specified as input and undesriable output simultaneously."
		error 498
	}	
	
	local comvars: list gopvars & bopvars
	if !(`"`comvars'"'==""){
		disp as error "`comvars' should not be specified as desriable and undesriable outputs simultaneously."
		error 498
	}	
		


	if "`nonradial'"==""{
		if `"`wmat'"'!=""{
			disp as red "Warning: wmat() should be only used when nonradial is specified."
			error 498
		}
	
	}
	else{
		tempname weightvec
		
		if `"`wmat'"'!=""{
			confirm matrix `wmat'
			local ncol=colsof(`wmat')
			if `ncol'!=`nvar'{
			    dis as error `"# of column of `wmat' != # of input-output variables"'
				exit 498
			}
			mat `weightvec'=`wmat'
			
			
		   }
		else{
			
			mat  `weightvec'=J(1,`nvar',1)
				
		}

		forval i = 1/`=colsof(`weightvec')' {
             local wvalues `wvalues'  `=`weightvec'[1,`i']'
           }
           	di 
            disp " The weight vector is (`wvalues')"
            //mat list `weightvec',noblank  noheader		
	}
	
   local techtype "contemporaneous production technology"
   //local techtype "contemporaneous"
   if `"`time'"'==""& "`sequential'`biennial'`window'"==""{
   		local techtype "global production technology"
   		local global global
   }

   if "`global'"==""  {

    	if "`time'"==""{
		   disp as error "For biennial/window/sequential/ model, time() should be specified."
		   error 498
		}

    }	
	
	if "`maxiter'"==""{
		local maxiter=-1
	}
	if "`tol'"==""{
		local tol=-1
	}	
	  

   if "`global'"!=""{
     if "`sequential'"!=""{
     
       disp as error "global and sequential cannot be specified together."
       error 498     
     
     }
     
     if "`window'"!=""{
     
       disp as error "global and window() cannot be specified together."
       error 498     
     
     }  

     if "`biennial'"!=""{
     
       disp as error "global and biennial cannot be specified together."
       error 498     
     
     }           
     
     local techtype "global production technology"
  
  } 
  
   

   if "`sequential'"!=""{
 
     if "`window'"!=""{
     
       disp as error "sequential and window() cannot be specified together."
       error 498     
     
     }  

     if "`biennial'"!=""{
     
       disp as error "sequential and biennial cannot be specified together."
       error 498     
     
     }            
     
     local techtype "sequential production technology"
  
  } 
    
 
     if "`window'"!=""{

       if "`biennial'"!=""{
       
         disp as error "biennial and window() cannot be specified together."
         error 498     
       
       }        
     
         local techtype "window(`window') production technology"   
     
     }  

       if "`biennial'"!=""{
       
        local techtype "biennial production technology"     
       
       }




	
	preserve
	

	
	
	if "`gx'"!=""{
		local ngx: word count `gx'
		if `ngx'!=`ninp'{
		    disp as error "# of input variables != # of variables specified in gx()."
		    error 498
		}
		local gmat `gmat' `gx'
		local gmatname `gmatname' `gx'
	
	}
	else{
		local invarscopy `invars' 
		forv k=1/`ninp'{
			gettoken word invarscopy:invarscopy
			//di "`word'"
		    tempvar gx_`k'
			qui gen `gx_`k''=-`word'
			local gmat `gmat' `gx_`k''
			local gmatname `gmatname' -`word'
		}
	
	}
	
	if "`gy'"!=""{
		local ngy: word count `gy'
		if `ngy'!=`ngo'{
		    disp as error "# of desriable output variables != # of variables specified in gy()."
		    error 498
		}
		local gmat `gmat' `gy'
		local gmatname `gmatname' `gy'
	
	}
	else{
	    local gopvarscopy `gopvars'
		forv k=1/`ngo'{
		    gettoken word gopvarscopy:gopvarscopy
		    tempvar gy_`k'
			qui gen `gy_`k''=`word'
			local gmat `gmat' `gy_`k''
			local gmatname `gmatname' `word'
		}
	
	}
		
	if "`gb'"!=""{
		local ngb: word count `gb'
		if `ngb'!=`nbo'{
		    disp as error "# of undesriable output variables != # of variables specified in gb()."
		    error 498
		}
		local gmat `gmat' `gb'
		local gmatname `gmatname' `gb'
	
	}
	else{
	    local bopvarscopy `bopvars'
		forv k=1/`nbo'{
		    gettoken word bopvarscopy:bopvarscopy
		    tempvar gb_`k'
			qui gen `gb_`k''=-`word'
			local gmat `gmat' `gb_`k''
			local gmatname `gmatname' -`word'

		}
	
	}	

        di 
        di " The diectional vector is (`gmatname')"
		
	local rts=0
	if "`vrs'"!=""{
	   local rts=1
	}
	
	

	qui keep   `invars' `gopvars' `bopvars' `dmu' `time' `gmat' `touse'
	qui gen Row=_n
	qui keep if `touse'
	label var Row "Row #"
	qui gen double Dval=.
	label var Dval "Value of NDDFs: `techtype'"
	

	
    tempvar tvar dmu2
	

	
	if  `"`time'"'!="" {
	    qui egen `tvar'=group(`time')
	    qui egen `dmu2'=group(`dmu')
	}
	else{
	    qui gen `tvar'=1
		qui gen `dmu2'=_n
	}
	
	sort Row `dmu2' `tvar' 


	tempvar touse2  touse3 touse4
    qui gen  byte `touse2'=1 
    markout `touse2' `invars' `opvars' `badvars'	
    qui gen byte `touse3'=0	
    qui gen byte `touse4'=`touse'	
	
	qui mata mata mlib index	
	local data `invars' `gopvars' `bopvars'
	if "`nonradial'"!=""{

		foreach v in `invars' `gopvars' `bopvars'{
		   qui gen double B_`v'=.
		   label var B_`v' `"beta:`v'"'
		   local slackvars `slackvars' B_`v'
		
		}

		if "`global'"!=""{

			mata: _tenddf("`data'",`ninp',`ngo',"`touse'", "`touse2'","`gmat'","`weightvec'","Dval `slackvars'",`rts',`maxiter',`tol')
		}


		qui su `tvar'
		local tmax=r(max)

		if "`sequential'"!=""{


	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'<=`t' & `touse2'==1) 

	 		   mata: _tenddf("`data'",`ninp',`ngo',"`touse4'", "`touse3'","`gmat'","`weightvec'","Dval `slackvars'",`rts',`maxiter',`tol')       	
	        }


		}	

		if "`biennial'"!=""{

	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'>=`t' & `tvar'<=`t'+1 & `touse2'==1) 
	 		   mata: _tenddf("`data'",`ninp',`ngo',"`touse4'", "`touse3'","`gmat'","`weightvec'","Dval `slackvars'",`rts',`maxiter',`tol')       	
	        }


		}


		if "`window'"!=""{
			local band=(`window'-1)/2
	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'>=`t'-`band' & `tvar'<=`t'+`band' & `touse2'==1) 
	 		   mata: _tenddf("`data'",`ninp',`ngo',"`touse4'", "`touse3'","`gmat'","`weightvec'","Dval `slackvars'",`rts',`maxiter',`tol')       	
	        }

		}


		if "`global'`sequential'`biennial'`window'"==""&"`time'"!=""{


	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'==`t' & `touse2'==1) 
	 		   mata: _tenddf("`data'",`ninp',`ngo',"`touse4'", "`touse3'","`gmat'","`weightvec'","Dval `slackvars'",`rts',`maxiter',`tol')       	
	        }


		}		


	}
	else{


		if "`global'"!=""{

			mata: _ddf("`data'",`ninp',`ngo',"`touse'", "`touse2'","`gmat'","Dval",`rts',`maxiter',`tol')
		}


		qui su `tvar'
		local tmax=r(max)

		if "`sequential'"!=""{


	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'<=`t' & `touse2'==1) 
	 		   mata: _ddf("`data'",`ninp',`ngo',"`touse4'", "`touse3'","`gmat'","Dval",`rts',`maxiter',`tol')       	
	        }


		}	

		if "`biennial'"!=""{

	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')	        	
	           qui replace `touse3'=(`tvar'>=`t' & `tvar'<=`t'+1 & `touse2'==1) 
	 		   mata: _ddf("`data'",`ninp',`ngo',"`touse4'", "`touse3'","`gmat'","Dval",`rts',`maxiter',`tol')       	
	        }


		}


		if "`window'"!=""{
			local band=(`window'-1)/2
	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')	        	
	           qui replace `touse3'=(`tvar'>=`t'-`band' & `tvar'<=`t'+`band' & `touse2'==1)
	 		   mata: _ddf("`data'",`ninp',`ngo',"`touse4'", "`touse3'","`gmat'","Dval",`rts',`maxiter',`tol')       	
	        }

		}

		if "`global'`sequential'`biennial'`window'"==""&"`time'"!=""{


	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')	        	
	           qui replace `touse3'=(`tvar'==`t' & `touse2'==1) 
	 		   mata: _ddf("`data'",`ninp',`ngo',"`touse4'", "`touse3'","`gmat'","Dval",`rts',`maxiter',`tol')   
	        }


		}	


	}

	
    qui keep if `touse'
	local NDDF=cond("`nonradial'"!="","Non-raidal","")

	format Dval  `slackvars' %9.4f
	order Row `dmu' `time' Dval  `slackvars'
	keep  Row `dmu' `time' Dval  `slackvars'
	
	disp _n(2) " `NDDF' Directional Distance Function Results:"
	disp "    (Row: Row # in the original data; Dval: Estimated value of DDF.)"

	list, sep(0) 
	di "Note: missing value indicates infeasible problem."
	
	if `"`saving'"'!=""{
	  save `saving'
	  gettoken filenames saving:saving, parse(",")
	  local filenames `filenames'.dta
	  disp _n `"Estimated Results are saved in `filenames'."'
	}

	return local file `filenames'

	
	restore 
	
	end
	
	
///////////////////////////////////////////////////////	


