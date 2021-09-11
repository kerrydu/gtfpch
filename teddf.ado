*! version 4.0.1, 11 Sep 2021
*! version 3.0.1, 20 June 2021
*! version 3.0, 20 Apr 2021
* use frame to store results
*! version 3.0, 6 Apr 2021
* add subsample bootstrap
*! version 2.1, 7 Aug 2020
*! version 2.0
* By Kerry Du, Daoping Wang, Ning Zhang, 30 Jul 2020 
**
* 
capture program drop teddf
program define teddf, rclass
    version 16


    **********************************************
    *******Check whether a new version is available
    // global c_m_d_0 is used in the updatecmd.ado
    global c_m_d_0 `0'  
    local pkg gtfpch // updatecmd should be replaced with your package name

    //install the updatecmd package if it is missing
    cap which updatecmd
    if _rc{
      cap net install updatecmd, from("https://github.com/kerrydu/kgitee/raw/master/") replace
      if _rc{
         cap net install updatecmd, from("https://gitee.com/kerrydu/kgitee/raw/master/") replace
         if _rc global up_grade_`pkg' "updatecmd_is_missing"
       }
      
    }
    local checkupdate 0
    //the first run of the command defines global up_grade_`pkg'
    if "${up_grade_`pkg'}"==""{ 
		local checkupdate 1
        updatecmd teddf, from("https://raw.githubusercontent.com/kerrydu/gtfpch/master/") ///
		froma("https://gitee.com/kerrydu/gtfpch/raw/master/") pkg(`pkg')       
    } 
	if `checkupdate'  exit
    ********************************************

	_get_version gtfpch
	_compile_mata, package(gtfpch) version(`package_version') verbose 

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
	
	
    syntax varlist [if] [in], Dmu(varname) [rf(varname) Time(varname) gx(varlist) gy(varlist) gb(varlist)  frame(name) ///
	                                        BIennial  SEQuential GLObal VRS  NONRadial  brep(integer 0) alpha(real 0.7) ///
										   Wmat(string) SAVing(string)  WINdow(numlist integer max=1 >=1)   level(real 95) ///
										   maxiter(numlist integer >0 max=1) tol(numlist max=1) NODOTS  noPRINT]
    local bopvars `varlist'

	if "`nodots'"!=""{
		local nodots = 1
	}
	else{
		local nodots = 0
	}
    
	if (`brep'<100 & `brep'>0) {
		di "{p 1 1 7}{error}# of  bootstrap replications must be at least 100."
		exit 198
	}
	if (`brep'<200 & `brep'>0) {
		di as red "Warning: Statistical inference may be unreliable for small number of bootstrap replications."
	}
	
	if (`level' < 10 |  `level' > 99.99) {
		di "{p 1 1 7}{error}level() must be between 10 and 99.99 inclusive{p_end}"
		exit 198
	}

	if (`alpha' < 0.5 |  `alpha' > 1) {
		di "{p 1 1 7}{error}alpha() must be between 0.5 and 1{p_end}"
		exit 198
	}

	if `"`frame'"'!=""{
    	confirm new frame `frame'
	}


    // copy data to ddfResults frame & saved resulst 
	qui pwf
    local currentframe = r(currentframe)

    *frame copy `currentframe' ddfResults
    *cwf ddfResults

	marksample touse 
	markout `touse' `invars' `gopvars' `gx' `gy' `gb'
	
	local invars: list uniq invars
	local gopvars: list uniq gopvars
	local bopvars: list uniq bopvars
	
	local ninp: word count `invars'
    local ngo: word count `gopvars'
    local nbo: word count `bopvars'
	local nvar=`ninp'+`ngo'+`nbo'
	
    local band = 0
	confirm numeric var `invars' `gopvars' `bopvars' `rf' `gx' `gy' `gb'

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
			disp as error "wmat() should be only used when nonradial is specified."
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
			forv i=1/`ncol'{
				local wival=`wmat'[1,`i']
				if `wival'<0{
					di as error `"The element of matrix `wmat' should not be less than 0."'
					exit 498
				}
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
			local rweightvec  (`wvalues')
            //mat list `weightvec',noblank  noheader		
	}
	
   local techtype "contemporaneous production technology"
   local ttype contemporaneous
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
     local ttype global
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
     local ttype sequential
  
  } 
    
 
     if "`window'"!=""{

       if "`biennial'"!=""{
       
         disp as error "biennial and window() cannot be specified together."
         error 498     
       
       }        
         local band = `window'
         local techtype "window(`window') production technology"   
         local ttype window
     
     }  

       if "`biennial'"!=""{
       
        local techtype "biennial production technology" 
        local ttype  biennial   
       
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
			*qui gen `gx_`k''=-`word'
			qui gen `gx_`k''=0
			local gmat `gmat' `gx_`k''
			*local gmatname `gmatname' -`word'
			local gmatname `gmatname' 0
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
		
	local rstype=0
	if "`vrs'"!=""{
	   local rstype=1
	}
	
	

	qui keep   `invars' `gopvars' `bopvars' `dmu' `time' `gmat' `touse' `rf' `gx' `gy' `gb'
	qui gen Row=_n
	
	label var Row "Row # in the original data"
	
    tempvar tvar dmu2
	
    qui egen `dmu2'=group(`dmu')
	
	if  `"`time'"'!="" {
	    qui egen `tvar'=group(`time')
	    
	}
	else{
	    qui gen `tvar'=1
	}
	
	sort `dmu2' `tvar'  // important for mm_sample() in the context of panel data


	tempvar touse2  
    qui gen byte `touse2'=1 
    if "`rf'"!=""{
    	qui replace `touse2'=(`rf'!=0) if !missing(`rf')
    }
    markout `touse2' `invars' `opvars' `badvars' `rf'
    *count if `touse2'	
	
	qui mata mata mlib index


	local data `invars' `gopvars' `bopvars'


    tempname data0 dataref gmat0

    qui keep if `touse'

    qui putmata `data0'   = (`dmu2' `tvar' `data') if `touse', replace
    qui putmata `gmat0'   = (`gmat') if `touse', replace
    qui putmata `dataref' = (`dmu2' `tvar' `data') if `touse2', replace


    qui gen double Dv=.
    label var Dv "value of distance function"
	local outputvar Dv 
    if `brep'>0{
        qui gen double Dv_bc=.
        label var Dv_bc "the bias-corrected Dv"
        qui gen double Dv_lower=.
        label var Dv_lower "the lower-bound estimate of Dv"
        qui gen double Dv_upper=.
        label var Dv_upper "the upper-bound estimate of Dv" 
        qui gen double bootvar=.
        label var bootvar "the bootstrap variance estimate" 
		  
    }

	

	if "`nonradial'"!=""{
		local Dvbeta Dv
		foreach v in `invars' `gopvars' `bopvars'{
		   qui gen double B_`v'=.
		   label var B_`v' `"beta:`v'"'
		   //local adjvars `adjvars' B_`v'
		   local Dvbeta `Dvbeta' B_`v'
		}

		 

		tempname wmat0
		mata: `wmat0'=st_matrix("`weightvec'")

		if `brep'==0{
			tempname DV
			mata: `DV'=nddf`ttype'(`data0',`dataref',`ninp',`ngo',`wmat0',`gmat0',`band',`rstype',`maxiter',`tol')
			qui getmata (`Dvbeta')= `DV', replace
			local outputvar `Dvbeta' 
		}
		else{
			mata: nddf`ttype'bt(`data0',`dataref',`ninp',`ngo',`wmat0',`gmat0',`band',`alpha',`level',`brep',`rstype',`maxiter',`tol',`nodots',"`Dvbeta'")
			local outputvar `Dvbeta' Dv_bc  Dv_lower  Dv_upper bootvar 
		}		


    }
	else{
		local outputvar Dv  
		if `brep'==0{
			tempname DV
			mata: `DV'=ddf`ttype'(`data0',`dataref',`ninp',`ngo',`gmat0',`band',`rstype',`maxiter',`tol')
			qui getmata Dv = `DV', replace
		}
		else{
			mata: ddf`ttype'bt(`data0',`dataref',`ninp',`ngo',`gmat0',`band',`alpha',`level',`brep',`rstype',`maxiter',`tol',`nodots')
			local outputvar Dv Dv_bc  Dv_lower  Dv_upper bootvar 
		}

	}

    local NDDF=cond("`nonradial'"!="","Non-raidal","")
 	format `outputvar'   %9.4f
	order Row `dmu' `time' `outputvar'
	keep  Row `dmu' `time' `outputvar'

    //di "`print'"
    if "`print'" != "noprint" {
        disp _n(2) " `NDDF' Directional Distance Function Results:"
        disp "    (Row: Row # in the original data; Dv: Estimated value of `NDDF' DDF.)"
        list, sep(0) 
        di "Note: Missing value indicates infeasible problem."
    }
	
	if `"`saving'"'!=""{
	  save `saving'
	  gettoken filenames saving:saving, parse(",")
	  local filenames `filenames'.dta
	  disp _n `"Estimated Results are saved in `filenames'."'
      *disp _n `"Estimated Results are also saved in ddfResults frame temporarily."'
      *cwf `currentframe'
      *frame drop ddfResults
	}

	if `"`frame'"'!=""{
		frame copy `currentframe' `frame'
		disp _n `"Estimated Results are also saved in `frame' frame temporarily."'
	}


	return local file `filenames'      
	return local frame `frame'
	return local gvec  (`gmatname')
	return local weight `rweightvec'

	restore 
	
	end
	
