*! version 2.0.1
* 2020-8-15, remove sup opt
*! version 2.0
* By Kerry Du & Daoping Wang, 8 Aug 2020 
capture program drop sbmeff
program define sbmeff, rclass
    version 16
	
    if strpos(`"`0'"',":")==0{  // Tone(2001) SBM without undesriable outputs
	   sbmeffnu2 `0'
	   exit
	   
	}	

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
	
	
    *syntax varlist [if] [in], dmu(varname) [Time(varname) SEQuential  VRS  SAVing(string) maxiter(numlist integer >0 max=1) tol(numlist max=1)]
	
	    syntax varlist [if] [in], Dmu(varname) [Time(varname)  rf(varname) ///
	                                        BIennial  SEQuential GLObal VRS  ///
										   SAVing(string)  WINdow(numlist intege max=1 >=1)    ///
										   maxiter(numlist integer >0 max=1) tol(numlist max=1) ]

										   
										   
	preserve
	marksample touse 
	markout `touse' `invars' `gopvars' 

	local bopvars `varlist'
	
	local invars: list uniq invars
	local gopvars: list uniq gopvars
	local bopvars: list uniq bopvars
	
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
		
	
	local ninp: word count `invars'
    local ngo: word count `gopvars'
    local nbo: word count `bopvars'
	
	local rts=0
	if "`vrs'"!=""{
	   local rts=1
	}
	
									   									   							 
	
	if "`maxiter'"==""{
		local maxiter=-1
	}
	if "`tol'"==""{
		local tol=-1
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
	
	
	
	qui mata mata mlib index
	
	
	qui keep  `invars' `gopvars' `bopvars' `dmu' `time' `touse' `rf'
	qui gen Row=_n
	label var Row "Row #"
	qui keep if `touse'	
	qui gen double TE=.
	label var TE "Technical Efficiency"
	
	foreach v in `invars' `gopvars' `bopvars'{
	   qui gen double S_`v'=.
	   label var S_`v' `"Slack:`v'"'
	   local slackvars `slackvars' S_`v'
	
	}
	
	

	
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
    if "`rf'"!=""{
    	qui replace `touse2'=(`rf'!=0) if !missing(`rf')
    }
    markout `touse2' `invars' `opvars' `badvars' `rf'

    qui gen byte `touse3'=0	
    qui gen byte `touse4'=`touse'	
	

	local data `invars' `gopvars' `bopvars'



		if "`global'"!=""{

			mata: _sbmeff("`data'",`ninp',`ngo',"`touse'", "`touse2'","TE `slackvars'",`rts',`maxiter',`tol')
		}


		qui su `tvar'
		local tmax=r(max)

		if "`sequential'"!=""{


	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'<=`t' & `touse2'==1) 

	 		   mata: _sbmeff("`data'",`ninp',`ngo',"`touse4'", "`touse3'","TE `slackvars'",`rts',`maxiter',`tol')       	
	        }


		}	

		if "`biennial'"!=""{

	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'>=`t' & `tvar'<=`t'+1 & `touse2'==1) 
	 		   mata: _sbmeff("`data'",`ninp',`ngo',"`touse4'", "`touse3'","Dval `slackvars'",`rts',`maxiter',`tol')       	
	        }


		}


		if "`window'"!=""{
			local band=(`window'-1)/2
	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'>=`t'-`band' & `tvar'<=`t'+`band' & `touse2'==1) 
	 		   mata: _sbmeff("`data'",`ninp',`ngo',"`touse4'", "`touse3'","TE `slackvars'",`rts',`maxiter',`tol')       	
	        }

		}


		if "`global'`sequential'`biennial'`window'"==""&"`time'"!=""{


	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'==`t' & `touse2'==1) 
	 		   mata: _sbmeff("`data'",`ninp',`ngo',"`touse4'", "`touse3'","TE `slackvars'",`rts',`maxiter',`tol')       	
	        }


		}		


	
	order Row `dmu' `time' TE `slackvars'
	keep Row `dmu' `time' TE `slackvars'
	
	disp _n(2) " SBM Efficiency Results:"
	disp "    (Row: Row # in the original data; TE: Efficiency Score;  S_X: Slack of X)"
	//disp "      S_X : Slack of X"
	*format TE `slackvars' %9.4f
	list Row `dmu' `time' TE `slackvars', sep(0) 

	//disp _n
	if `"`saving'"'!=""{
	  save `saving'
	  gettoken filenames saving:saving, parse(",")
	  local filenames `filenames'.dta
	  disp _n `"Estimated Results are saved in `filenames'."'
	}
	
	return local file `filenames'
	restore 
	
	end
	
///////////////////////////////////////////////////////////////////////////////	
capture program drop sbmeffnu2
program define sbmeffnu2, rclass
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


   
	
	
    *syntax varlist [if] [in], dmu(varname) [Time(varname) SEQuential  VRS  SAVing(string) maxiter(numlist integer >0 max=1) tol(numlist max=1)]
	
	    syntax varlist [if] [in], Dmu(varname) [Time(varname)  rf(varname) ///
	                                        BIennial  SEQuential GLObal VRS  ///
										   SAVing(string)  WINdow(numlist intege max=1 >=1)    ///
										   maxiter(numlist integer >0 max=1) tol(numlist max=1) ]

	local gopvars `varlist'		
	unab gopvars : `gopvars'
										   
	preserve
	marksample touse 
	markout `touse' `invars' `gopvars' 

	local invars: list uniq invars
	local gopvars: list uniq gopvars

	
	confirm numeric var `invars' `gopvars' 
	
	local comvars: list invars & gopvars 
	if !(`"`comvars'"'==""){
		disp as error "`comvars' should not be specified as input and desriable output simultaneously."
		error 498
	}
	
	
	local ninp: word count `invars'
    local ngo: word count `gopvars'
   
	
	local rts=0
	if "`vrs'"!=""{
	   local rts=1
	}
	
									   									   							 
	
	if "`maxiter'"==""{
		local maxiter=-1
	}
	if "`tol'"==""{
		local tol=-1
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
	
	
	
	qui mata mata mlib index
	
	
	qui keep  `invars' `gopvars'  `dmu' `time' `touse' `rf'
	qui gen Row=_n
	label var Row "Row #"
	qui keep if `touse'	
	qui gen double TE=.
	label var TE "Technical Efficiency"
	
	foreach v in `invars' `gopvars' {
	   qui gen double S_`v'=.
	   label var S_`v' `"Slack:`v'"'
	   local slackvars `slackvars' S_`v'
	
	}
	
	
	
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
    
    if "`rf'"!=""{
    	qui replace `touse2'=(`rf'!=0) if !missing(`rf')
    }
    markout `touse2' `invars' `opvars' `rf'
   * markout `touse2' `invars' `opvars' 
    qui gen byte `touse3'=0	
    qui gen byte `touse4'=`touse'	
	

	local data `invars' `gopvars'



		if "`global'"!=""{

			mata: _sbm("`data'",`ninp',`ngo',"`touse'", "`touse2'","TE `slackvars'",`rts',`maxiter',`tol')
		}


		qui su `tvar'
		local tmax=r(max)

		if "`sequential'"!=""{


	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'<=`t' & `touse2'==1) 

	 		   mata: _sbm("`data'",`ninp',`ngo',"`touse4'", "`touse3'","TE `slackvars'",`rts',`maxiter',`tol')       	
	        }


		}	

		if "`biennial'"!=""{

	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'>=`t' & `tvar'<=`t'+1 & `touse2'==1) 
	 		   mata: _sbm("`data'",`ninp',`ngo',"`touse4'", "`touse3'","Dval `slackvars'",`rts',`maxiter',`tol')       	
	        }


		}


		if "`window'"!=""{
			local band=(`window'-1)/2
	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'>=`t'-`band' & `tvar'<=`t'+`band' & `touse2'==1) 
	 		   mata: _sbm("`data'",`ninp',`ngo',"`touse4'", "`touse3'","TE `slackvars'",`rts',`maxiter',`tol')       	
	        }

		}


		if "`global'`sequential'`biennial'`window'"==""&"`time'"!=""{


	        forv t=1/`tmax'{
	           qui replace `touse4'=(`tvar'==`t' & `touse')
	           qui replace `touse3'=(`tvar'==`t' & `touse2'==1) 
	 		   mata: _sbm("`data'",`ninp',`ngo',"`touse4'", "`touse3'","TE `slackvars'",`rts',`maxiter',`tol')       	
	        }


		}		


	
	order Row `dmu' `time' TE `slackvars'
	keep Row `dmu' `time' TE `slackvars'
	*format TE `slackvars' %9.4f
	disp _n(2) " SBM Efficiency Results:"
	disp "    (Row: Row # in the original data; TE: Efficiency Score;  S_X: Slack of X)"
	//disp "      S_X : Slack of X"
	list Row `dmu' `time' TE `slackvars', sep(0) 

	//disp _n
	if `"`saving'"'!=""{
	  save `saving'
	  gettoken filenames saving:saving, parse(",")
	  local filenames `filenames'.dta
	  disp _n `"Estimated Results are saved in `filenames'."'
	}
	
	return local file `filenames'
	restore 
	
	end	
	
	

