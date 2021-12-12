*! version 4.1, 21 Oct 2021
*! version 4.022, 18 Sep 2021
*! version 2.02, 20 June 2021
*! version 2.01, 20 Apr 2021
* use frame to store results
*! version 1.2
* 19 Jul 2020
* added hybrid Oriented
*! version 1.1
* 30 Apr 2020
* add biennial option

*! version 1.0
* By Kerry Du,Daoping Wang, Ning Zhang, 3 Dec 2019 
**
* 
capture program drop gtfpch
program define gtfpch, rclass prop(xt)
    version 16

   * recompile the mlib for different Stata version
	  _get_version gtfpch
	  _compile_mata, package(gtfpch) version(`package_version') verbose 
    qui pwf
    local currentframe = r(currentframe)

    qui mata mata mlib index
    _xt, trequired 
    local id=r(ivar)
    local time=r(tvar)

   * extract input and desriable output variables
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
	
    syntax varlist [if] [in],    [dmu(varname) gx(varlist) gy(varlist) gb(varlist)  ///
	                              BIennial  SEQuential GLOBAL FGNZ RD  LUENberger ort(string) ///
								  WINdow(numlist intege max=1 >=1) SAVing(string)   Wmat(string) NONRadial ///
								  maxiter(numlist integer >0 max=1) tol(numlist max=1) frame(name) noPRINT NOCHeck]


   * check whether the new version is available
   if "`nocheck'"==""{ //default
      if "$gtfpchcheck"==""   gtfpchcheckupdate gtfpch 
   }
   
    preserve

    sort `id' `time'
    if `"`frame'"'!=""{
       confirm new frame `frame'
    }    
    
    marksample touse

    local bopvars `varlist'
    local ninp: word count `invars'
    local ngo:  word count `gopvars'
    local nbo:  word count `bopvars'
    local nvar: word count `invars' `gopvars' `bopvars'
    qui keep `invars' `gopvars' `bopvars' `touse' `id' `time' `dmu'  `gx' `gy' `gb'

    qui gen Row=_n
    label var Row "Row # in the original dataset"
    if "`nonradial'"!=""{
     local luenberger luenberger
	 *if "`ort'" =="" local ort h

     }

    if ("`ort'" == "") local ort = "out"
    else if ("`ort'" == "hybrid" | "`ort'" == "HYBRID" | "`ort'" == "h" | "`ort'" == "H") {
         local ort="h"
    }
    else {
        local ort = upper("`ort'")
        if ("`ort'" == "I" | "`ort'" == "IN" | "`ort'" == "INPUT") {
            local ort = "i"
        }
        else if ("`ort'" == "O" | "`ort'" == "OUT" | "`ort'" == "OUTPUT") {
            local ort = "out"
        }
        else {
            di as err "option ort allows for case-insensitive " _c
            di as err "(i|in|input|o|out|output|h|hybrid) or nothing."
            exit 198
        }
    }

  if  "`luenberger'"=="" & "`ort'"=="h"{
      di as err "ort(hybrid) can not be used for Malmquist-Luenberger productivity index."
      exit 198
  }


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
  
    if  `"`ort'"'!="out"{ // for ort(i) and ort(h)
      local invarscopy `invars'
      forv k=1/`ninp'{
        gettoken word invarscopy:invarscopy
        tempvar gx_`k'
        qui gen `gx_`k''=-`word'  // gx=-X for input orientied direction
        local gmat `gmat' `gx_`k''
        local gx `gx' `gx_`k''
        local gmatname `gmatname' -`word'
      }   
    
    }
    else{
      forv k=1/`ninp'{
        tempvar gx_`k'
        qui gen `gx_`k''=0  // gx=0 for output orientied direction
        local gmat `gmat' `gx_`k''
        local gx `gx' `gx_`k''
        local gmatname `gmatname' 0 
      }
    
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
  
    if  `"`ort'"'!="i"{ //for ort(out) and ort(h)
      local gopvarscopy `gopvars'
      forv k=1/`ngo'{
        gettoken word gopvarscopy:gopvarscopy
        tempvar gy_`k'
        qui gen `gy_`k''=`word'
        local gmat `gmat' `gy_`k''
        local gy `gy' `gy_`k''
        local gmatname `gmatname' `word'       
      }   
    
    }
    else{
      forv k=1/`ngo'{
        tempvar gy_`k'
        qui gen `gy_`k''=0
        local gmat `gmat' `gy_`k''
        local gy `gy' `gy_`k''
        local gmatname `gmatname' 0        
      }       
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
  
    if `"`ort'"'!="i"{
  
    local bopvarscopy `bopvars'
    forv k=1/`nbo'{
      gettoken word bopvarscopy:bopvarscopy
      tempvar gb_`k'
      qui gen `gb_`k''=-`word'
      local gmat `gmat' `gb_`k''
      local gb `gb' `gb_`k''
      local gmatname `gmatname' -`word'     
    } 
  
    }
    else{
      forv k=1/`nbo'{
        tempvar gb_`k'
        qui gen `gb_`k''=0
        local gmat `gmat' `gb_`k''
        local gb `gb' `gb_`k''
        local gmatname `gmatname'  0         
      }     
    
    
    }

  
  } 



    tempname weightvec wcheck
    if "`nonradial'"==""{
        if `"`wmat'"'!=""{
            disp as red "Warning: wmat() is only used when nonradial is specified."
        }
        
    }
    else{
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
              di as error `"The element of matrix `wmat' should not be negative."'
              exit 498
            }
          }
           
        //mat `weightvec'=`wmat'
        if "`ort'"=="out"{
            mat `weightvec'=`wmat'[1,1..`ninp']
            mata: st_numscalar("`wcheck'",sum(st_matrix("`weightvec'")))
            if `wcheck'!=0{
                disp as error "For output orientied Luenberger productivity indicator, the weight components for input reduction should be set to 0."
                exit 498
            }
        }
        else if "`ort'"=="i"{
            mat `weightvec'=`wmat'[1,(`ninp'+1)..(`ninp'+`nbo'+`ngo')]
            mata: st_numscalar("`wcheck'",sum(st_matrix("`weightvec'")))
            if `wcheck'!=0{
                disp as error "For input orientied Luenberger productivity indicator, the weight components for output reduction should be set to 0."
                exit 498
            }
            
        }       
         mat `weightvec'=`wmat'
        }
        else{
			
           if "`ort'"=="out"{       
              mat  `weightvec'=(J(1,`ninp',0),J(1,`ngo',1),J(1,`nbo',1))           
           }
           else if "`ort'"=="i"{
              mat  `weightvec'=(J(1,`ninp',1),J(1,`ngo'+`nbo',0))               
           }
           else{
              mat  `weightvec'=(J(1,`ninp',1),J(1,`ngo'+`nbo',1)) 
           }

        
        }

        forval i = 1/`=colsof(`weightvec')' {
             local wvalues `wvalues'  `=`weightvec'[1,`i']'
           }
        di 
        disp " The weight vector is (`wvalues')"
        local rweightname (`wvalues')

    }   


************************************** CRS case*************************************

    if "`nonradial'"==""{

        if "`luenberger'"==""{
            _malmqluen `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') gx(`gx') ///
                               gy(`gy') gb(`gb')  ort(`ort')   ///
                              `biennial' `global' `sequential' window(`window') ///
                               maxiter(`maxiter') tol(`tol')
            local indexname "Malmquist-Luenberger Productivity Index"

        }
        else{

            _luen_ddf `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') gx(`gx') ///
                               gy(`gy') gb(`gb')  ort(`ort')   ///
                              `biennial' `global' `sequential' window(`window') ///
                               maxiter(`maxiter') tol(`tol')
            local indexname "Luenberger Productivity Index (based on DDF)"
        }
    }
    else{

             _luen_nddf `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') gx(`gx') ///
                               gy(`gy') gb(`gb')  ort(`ort')  wmat(`weightvec') ///
                              `biennial' `global' `sequential' window(`window') ///
                               maxiter(`maxiter') tol(`tol')       
            local indexname "Luenberger Productivity Index (based on nonradial DDF)"

    }
            local resvars `r(rvars)'
*****************************************************************************************
  if "`rd'"==""&"`fgnz'"==""{

        di 
        di " The directional vector is (`gmatname')"
        format `resvars' %9.4f
        qui keep if  `touse'
        qui cap bys `id' (`time'): gen Pdwise=`time'[_n-1]+"~"+`time' if _n>1
        qui cap bys `id' (`time'): gen Pdwise=string(`time'[_n-1])+"~"+string(`time') if _n>1
        label var Pdwise "Period wise"       
        order Row `dmu' `id' Pdwise  `resvars' 
        qui keep  Row `dmu' `id' Pdwise  `resvars'    
        qui keep if !missing(Pdwise)

        if "`noprint'" == ""{
          disp _n(2) " Total Factor Productivity Change:`indexname'"
          disp "    (Row: Row # in the original data; Pdwise: periodwise)"
          list Row `dmu' `id'  Pdwise  `resvars', sep(0) 
          di "Note: missing value indicates infeasible problem."
        }


        if `"`saving'"'!=""{
          save `saving'
          gettoken filenames saving:saving, parse(",")
          local filenames `filenames'.dta
          disp _n `"Estimated Results are saved in `filenames'."'
       
        }

        if `"`frame'"'!=""{
          disp _n `"Estimated Results are saved in `frame' frame."'
          frame copy `currentframe'  `frame'    
        }           
        
        return local file `filenames'  
        return local frame `frame'
        *return local techtype `techtype'     
        return local weight  `rweightname'
        return local gvec  `gmatname'    
        restore

  }


  **************************** estiamte VRS and decompose scale effciency ************************
  else{


   foreach v of local resvars{
    
      rename `v' `v'_crs
    
    }

   if "`nonradial'"==""{

        if "`luenberger'"==""{
            _malmqluen `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') gx(`gx') ///
                               gy(`gy') gb(`gb')  ort(`ort')  vrs ///
                              `biennial' `global' `sequential' window(`window') ///
                               maxiter(`maxiter') tol(`tol')
            //local indexname "Malmquist-Luenberger Productivity Index"

        }
        else{

            _luen_ddf `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') gx(`gx') ///
                               gy(`gy') gb(`gb')  ort(`ort')  vrs  ///
                             `biennial' `global' `sequential' window(`window') ///
                               maxiter(`maxiter') tol(`tol')
            //local indexname "Luenberger Productivity Index (base on DDF)"
        }
    }
    else{

             _luen_nddf `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') gx(`gx') ///
                               gy(`gy') gb(`gb')  ort(`ort')  wmat(`weightvec') ///
                              `biennial' `global' `sequential' window(`window') vrs ///
                               maxiter(`maxiter') tol(`tol')       
            //local indexname "Luenberger Productivity Index (base on nonradial DDF)"

    }
    local resvars `r(rvars)'

    if "`luenberger'"==""{
      if "`rd'"!=""{
      
        qui gen SECH=TFPCH_crs/TFPCH
        label var SECH "Scale efficiecny change"
        qui replace TFPCH=TFPCH_crs
        local resvars `resvars' SECH
      
      }
      else{
        qui gen SECH=TECH_crs/TECH
        label var SECH "Scale efficiecny change"      
        qui replace TFPCH=TFPCH_crs
		qui replace TECCH=TECCH_crs 
        local resvars `resvars' SECH
      }

    }
    else{
      if "`rd'"!=""{
        
        //su TFPCH_crs TFPCH
        qui gen SECH=TFPCH_crs-TFPCH
        label var SECH "Scale efficiecny change"
        qui replace TFPCH=TFPCH_crs
        local resvars `resvars' SECH
      
      }
      else{
        //su TECH_crs TECH
        qui gen SECH=TECH_crs-TECH
        label var SECH "Scale efficiecny change"      
        qui replace TFPCH=TFPCH_crs
		
		qui replace TECCH=TECCH_crs 

        local resvars `resvars' SECH
      }

    }

    
    di 
    di " The directional vector is (`gmatname')"    
    format `resvars' %9.4f
	
    qui keep if `touse'
    qui cap bys `id' (`time'): gen Pdwise=`time'[_n-1]+"~"+`time' if _n>1
    qui cap bys `id' (`time'): gen Pdwise=string(`time'[_n-1])+"~"+string(`time') if _n>1
    label var Pdwise "Period wise"
          
    order Row `dmu' `id' Pdwise  `resvars' 
    qui keep if !missing(Pdwise) & `touse'
    qui keep  Row `dmu' `id' Pdwise  `resvars' 

    if "`noprint'"==""{
      disp _n(2) " Total Factor Productivity Change:`indexname'"
      disp "    (Row: Row # in the original data; Pdwise: periodwise)"

      list Row `dmu' `id'  Pdwise  `resvars' , sep(0) 
      di "Note: missing value indicates infeasible problem."          
    }

    if `"`saving'"'!=""{
      save `saving'
      gettoken filenames saving:saving, parse(",")
      local filenames `filenames'.dta
      disp _n `"Estimated Results are saved in `filenames'."'

    } 
        if `"`frame'"'!=""{
          disp _n `"Estimated Results are saved in `frame' frame."'
          frame copy `currentframe'  `frame'    
        }           
        
        return local file `filenames'  
        return local frame `frame'   
        *return local techtype `techtype'     
        return local weight  `rweightname'
        return local gvec  `gmatname'          

        restore



  }
        


end


***********************************************************
/////////////////////subprograms///////////////////////////
***********************************************************
cap program drop _malmqluen
program define _malmqluen,rclass
    version 16
  
    gettoken word 0 : 0, parse("=:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
    local invars `invars' `word'
    gettoken word 0 : 0, parse("=:,")
    }
    unab invars : `invars'

  gettoken word 0 : 0, parse("=:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
    local gopvars `gopvars' `word'
    gettoken word 0 : 0, parse(" =:,")
    }
    unab gopvars : `gopvars'

    syntax varlist [if] [in], id(varname) time(varname)  [gx(varlist)    ///
                       BIennial  gy(varlist) gb(varlist) VRS ort(string) ///
                       GLOBAL SEQuential WINdow(numlist intege max=1 >=1) ///
                  maxiter(numlist integer >0 max=1) tol(numlist max=1 >0)   ]
                 
  marksample touse 
  local bopvars `varlist'
  
  local ninp: word count `invars'
  local ngo:  word count `gopvars'
  local nbo:  word count `bopvars'
  local nvar: word count `invars' `gopvars' `bopvars'

  confirm numeric var `invars' `gopvars' `bopvars'
  
  local techtype "contemporaneous"
   

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
     
     local techtype "global"
  
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
     
     local techtype "sequential"
  
  } 
    
 
     if "`window'"!=""{

       if "`biennial'"!=""{
       
         disp as error "biennial and window() cannot be specified together."
         error 498     
       
       }        
     
         local techtype "window"   
     
     }  

       if "`biennial'"!=""{
       
        local techtype "biennial"     
       
       }  




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


     
  if "`maxiter'"==""{
    local maxiter=-1
  }
  if "`tol'"==""{
    local tol=-1
  } 
  
   
  tempvar period dmu
  
  qui egen `period'=group(`time')
  qui egen `dmu'=group(`id')  

  
    qui su  `period'
    local tmax=r(max)

    tempvar flag temp DD D21 D12
    
    qui gen `DD'=.
    qui gen `D21'=.
    qui gen `D12'=.
    qui gen `temp'=.
    qui gen `flag'=0
  
    local gmat `gx' `gy' `gb'
    sort `period' `dmu'
  
  if `"`techtype'"'=="contemporaneous"{
  
      qui{
	  	*D11
        forv t=1/`tmax'{ 
            qui replace `flag'= (`period'==`t')
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }
		*D21
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'==`t'+1) 
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`D21') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }  
		*D12
        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1)
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }       

    }
  
  
  }
  
  
    if `"`techtype'"'=="sequential"{
  
		*D11
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t')
            _ddf if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')

        }
		*D21
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'<=`t'+1) 
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`D21') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }  
		*D12
        forv t=2/`tmax'{
            qui replace `flag'= (`period'<=`t'-1) 
            _ddf if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }       

    
  
  
  }
  
  
  
     if `"`techtype'"'=="window"{
      
       local band=`window'
   
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t'+`band' & `period'>=`t'-`band') 
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'= (`period'<=`t'+1+`band' &  `period'>=`t'-`band'+1)
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`D21') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'<=`t'-1+`band' & `period'>=`t'-1-`band') 
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }       

    
  
  
  }



  if `"`techtype'"'=="biennial" {
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t' | `period'==`t'+1)
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }     

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1 | `period'==`t')
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')

        }       

    }
  
  
  }

 

  
 
  if `"`techtype'"'=="global"{

      qui replace `flag'=1
      _ddf  if `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol') 
      qui bys `dmu' (`period'): gen TFPCH=`temp'/`temp'[_n-1] 
      label var TFPCH "Total factor productivity change"
      cap drop `temp'   
    
    sort `period' `dmu'
    forv t=1/`tmax'{
      qui replace `flag'=(`period'==`t')
      _ddf  if `touse' & `period'==`t', rflag(`flag') gen(`DD') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
    }
    
    qui bys `dmu' (`period'): gen TECH=`DD'/`DD'[_n-1]  
  
	qui bys `dmu' (`period'): gen TECCH=TFPCH/TECH
  
    label var TECH  "Technical efficiency change" 
	label var TECCH "Technological change"
   
	local resvars TFPCH TECH  TECCH
    
    
  }

  else if `"`techtype'"'=="biennial"{

    qui bys `dmu' (`period'): gen TFPCH=`D12'/`DD'[_n-1] if _n>1
    label var TFPCH "Total factor productivity change"
    cap drop `temp' 

    sort `period' `dmu'
    forv t=1/`tmax'{
      qui replace `flag'=(`period'==`t')
      _ddf  if `touse' & `period'==`t', rflag(`flag') gen(`DD') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')

    }

    qui bys `dmu' (`period'): gen TECH=`DD'/`DD'[_n-1]  
    qui bys `dmu' (`period'): gen TECCH=TFPCH/TECH      
  
    label var TECH  "Technical efficiency change" 
    label var TECCH "Techological change"
    local resvars TFPCH TECH  TECCH           

  }

  else{
  
      //su `DD' `D12' `D21'
    qui {
      sort `dmu' `period'
      bys `dmu' (`period'): gen TECH=`DD'/`DD'[_n-1]
      bys `dmu' (`period'): gen TECCH=sqrt(`D12'/`DD'*`DD'[_n-1]/`D21'[_n-1])
      gen TFPCH= TECH*TECCH
      local resvars TFPCH TECH TECCH
        label var TFPCH "Total factor productivity change"
        label var TECH  "Technical efficiency change" 
        label var TECCH "Techological change"       
    } 
  
  }
  
  return local rvars "`resvars'"
  


end  

///////////////////////////////////////////////////////////////////////////
capture program drop _ddf
program define _ddf
    version 16

    syntax [if] [in], gen(string) INvars(varlist) OPvars(varlist) BADvars(varlist) [GVec(varlist) rflag(varname) ort(string) VRS maxiter(numlist) tol(numlist)]
        

    marksample touse 
    markout `touse' `invars' `opvars' `badvars' `gvec'
    
    tempvar touse2
    local rflag=cond("`rflag'"=="","","if `rflag'")
    mark `touse2' `rflag' // `rflag' might be empty
    
    markout `touse2' `invars' `opvars' `badvars'

     local data `invars' `opvars' `badvars'
     local num1: word count `invars'   
     local num2: word count `opvars'
      
    if "`vrs'"!=""{
      local rts=1
    }
    else{
      local rts=0
    
    }
    
    mata: _ddf("`data'",`num1',`num2',"`touse'", "`touse2'","`gvec'","`gen'",`rts',`maxiter',`tol')
	
	//di "`ort'"

    if "`ort'" =="out"{
      
      qui replace `gen'=1/(1+`gen') if `touse'
      }
    else{
      qui replace `gen'=(1-`gen') if `touse'
     }
     

end 

//////////////////////////////////////////////////
cap program drop _luen_ddf
program define _luen_ddf,rclass
    version 16
  

    gettoken word 0 : 0, parse("=:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
    local invars `invars' `word'
    gettoken word 0 : 0, parse("=:,")
    }
    unab invars : `invars'

  gettoken word 0 : 0, parse("=:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
    local gopvars `gopvars' `word'
    gettoken word 0 : 0, parse(" =:,")
    }
    unab gopvars : `gopvars'

********************************************************************************* 
  
    *local ninv: word count `invars'
    *local ngov: word count `gopvars'
    syntax varlist [if] [in], id(varname) time(varname)  [gx(varlist) ///
                   BIennial   gy(varlist) gb(varlist) VRS ort(string) ///
                             GLOBAL SEQuential WINdow(numlist intege max=1 >=1) ///
                 maxiter(numlist integer >0 max=1) tol(numlist max=1 >0)]
                 
  marksample touse 
  local bopvars `varlist'
  local ninp: word count `invars'
  local ngo: word count `gopvars'
  local nbo: word count `bopvars'
  local nvar: word count `invars' `gopvars' `bopvars' 
  confirm numeric var `invars' `gopvars' `bopvars'
 
 local techtype "contemporaneous"
   

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
     
     local techtype "global"
  
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
     
     local techtype "sequential"
  
  } 
    
 
     if "`window'"!=""{

       if "`biennial'"!=""{
       
         disp as error "biennial and window() cannot be specified together."
         error 498     
       
       }        
     
         local techtype "window"   
     
     }  

       if "`biennial'"!=""{
       
        local techtype "biennial"     
       
       }    
  
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


     
  if "`maxiter'"==""{
    local maxiter=-1
  }
  if "`tol'"==""{
    local tol=-1
  } 
  
   
    tempvar period dmu
  
  qui egen `period'=group(`time')
  qui egen `dmu'=group(`id')  

  
    qui su  `period'
    local tmax=r(max)

    tempvar flag temp DD D21 D12
    
    qui gen `DD'=.
    qui gen `D21'=.
    qui gen `D12'=.
    qui gen `temp'=.
    local gmat `gx' `gy' `gb'
  
    qui gen `flag'=0
  
    sort `period' `dmu'
  
  if `"`techtype'"'=="contemporaneous"{
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t')
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
  
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'==`t'+1) 
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`D21') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
  
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1)
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
    
        }       

    }
  
  
  }


  if `"`techtype'"'=="biennial"{
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t' | `period'==`t'+1)
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
   
        }    

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1 | `period'==`t')
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
       
        }       

    }
  
  
  }

  
  
    if `"`techtype'"'=="sequential"{
  
    
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t')
            _ddf_luen if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'<=`t'+1) 
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`D21') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        }  

        forv t=2/`tmax'{
            qui replace `flag'= (`period'<=`t'-1) 
            _ddf_luen if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')

        }       

    
  
  
  }
  
  
  
     if `"`techtype'"'=="window"{
    *local band=(`window'-1)/2
      local band=`window'
   
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t'+`band' & `period'>=`t'-`band') 
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
  
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'= (`period'<=`t'+1+`band' &  `period'>=`t'-`band'+1)
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`D21') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
    
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'<=`t'-1+`band' & `period'>=`t'-1-`band') 
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            
        }       

    
  
  
  }
 

  
 
  if `"`techtype'"'=="global"{

      qui replace `flag'=1
    _ddf_luen  if `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        
        qui bys `dmu' (`period'): gen TFPCH=`temp'[_n-1]-`temp' 
    label var TFPCH "Total factor productivity change"
    cap drop `temp'   
    
    sort `period' `dmu'
    forv t=1/`tmax'{
      qui replace `flag'=(`period'==`t')
      _ddf_luen  if `touse' & `period'==`t', rflag(`flag') gen(`DD') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
   
    }
    
    qui bys `dmu' (`period'): gen TECH=`DD'[_n-1]-`DD'
    *qui bys `dmu' (`period'): gen BPC=TFPCH-TECH 
	qui bys `dmu' (`period'): gen TECCH=TFPCH-TECH
  
    label var TECH  "Technical efficiency change" 
    *label var BPC "Best practice gap change"
	label var TECCH "Technological change"
    *local resvars TFPCH TECH  BPC
	local resvars TFPCH TECH  TECCH
    
    
  }
  else if `"`techtype'"'=="biennial"{

     qui bys `dmu' (`period'): gen TFPCH=`DD'[_n-1]-`D12' if _n>1
     label var TFPCH "Total factor productivity change"
     qui cap drop `temp'
     sort `period' `dmu'
    forv t=1/`tmax'{
      qui replace `flag'=(`period'==`t')
      _ddf_luen  if `touse' & `period'==`t', rflag(`flag') gen(`DD') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
     
    }
    
    qui bys `dmu' (`period'): gen TECH=`DD'[_n-1]-`DD'
    qui bys `dmu' (`period'): gen TECCH=TFPCH-TECH      
  
    label var TECH  "Technical efficiency change" 
    label var TECCH "Techological change" 
    local resvars TFPCH TECH  TECCH     


 }

  else{
  
      //su `DD' `D12' `D21'
    qui {
      sort `dmu' `period'
      bys `dmu' (`period'): gen TECH=`DD'[_n-1]-`DD'
      bys `dmu' (`period'): gen TECCH=0.5*(`DD'-`D12'+`D21'[_n-1]-`DD'[_n-1])
      gen TFPCH= TECH+TECCH
      local resvars TFPCH TECH TECCH
        label var TFPCH "Total factor productivity change"
        label var TECH  "Technical efficiency change" 
        label var TECCH "Techological change"       
    } 
  
  }
  
  return local rvars "`resvars'"
  


end  



///////////////////////////////////////////////////////////////
capture program drop _ddf_luen
program define _ddf_luen
    version 16

    syntax [if] [in], gen(string) INvars(varlist) OPvars(varlist) BADvars(varlist) [GVec(varlist) rflag(varname)  VRS maxiter(numlist) tol(numlist)]
        

    marksample touse 
    markout `touse' `invars' `opvars' `badvars' `gvec'
    local rflag=cond("`rflag'"=="","","if `rflag'") //fix missing cond 18-9-2021
    tempvar touse2
    mark `touse2' `rflag' // `rflag' might be empty
    markout `touse2' `invars' `opvars' `badvars'
 

    local data `invars' `opvars' `badvars'
    local num1: word count `invars'   
    local num2: word count `opvars'
    
  
    if "`vrs'"!=""{
      local rts=1
    }
    else{
      local rts=0
    
    }
    
    mata: _ddf("`data'",`num1',`num2',"`touse'", "`touse2'","`gvec'","`gen'",`rts',`maxiter',`tol')

     

end 

////////////////////////////////////////////////////////////////////////
cap program drop _luen_nddf
program define _luen_nddf,rclass
    version 16
  

    gettoken word 0 : 0, parse("=:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
    local invars `invars' `word'
    gettoken word 0 : 0, parse("=:,")
    }
    unab invars : `invars'

  gettoken word 0 : 0, parse("=:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
    local gopvars `gopvars' `word'
    gettoken word 0 : 0, parse(" =:,")
    }
    unab gopvars : `gopvars'

  
    *local ninv: word count `invars'
    *local ngov: word count `gopvars'
    syntax varlist [if] [in], id(varname) time(varname) wmat(string)  [ gx(varlist) ///
                           BIennial  gy(varlist) gb(varlist) VRS ort(string) ///
                             GLOBAL SEQuential WINdow(numlist intege max=1 >=1) ///
                 maxiter(numlist integer >0 max=1) tol(numlist max=1 >0)]
                 
  marksample touse 
  local bopvars `varlist'

  local ninp: word count `invars'
  local ngo: word count `gopvars'
  local nbo: word count `bopvars'
  local nvar: word count `invars' `gopvars' `bopvars' 
  confirm numeric var `invars' `gopvars' `bopvars'
 local techtype "contemporaneous"
   

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
     
     local techtype "global"
  
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
     
     local techtype "sequential"
  
  } 
    
 
     if "`window'"!=""{

       if "`biennial'"!=""{
       
         disp as error "biennial and window() cannot be specified together."
         error 498     
       
       }        
     
         local techtype "window"   
     
     }  

       if "`biennial'"!=""{
       
        local techtype "biennial"     
       
       }   
  
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

/*
  tempname weightvec wcheck
  if "`nonradial'"==""{
    if `"`wmat'"'!=""{
      disp as red "Warning: wmat() is only used when nonradial is specified."
    }
    
  }
  else{
    if `"`wmat'"'!=""{  
      confirm matrix `wmat'
      local ncol=colsof(`wmat')
      if `ncol'!=`nvar'{
        dis as error `"# of column of `wmat' != # of input-output variables"'
        exit 498
      }
      
    //mat `weightvec'=`wmat'
    if "`ort'"=="out"{
      mat `weightvec'=`wmat'[1,1..`ninp']
      mata: st_numscalar("`wcheck'",sum(st_matrix("`weightvec'")))
      if `wcheck'!=0{
        disp as error "For output orientied Luenberger productivity index, the weight components for input reduction should be set to 0."
        exit 498
      }
    }
    else{
      mat `weightvec'=`wmat'[1,(`ninp'+1)..(`ninp'+`nbo'+`ngo')]
      mata: st_numscalar("`wcheck'",sum(st_matrix("`weightvec'")))
      if `wcheck'!=0{
        disp as error "For input orientied Luenberger productivity index, the weight components for output reduction should be set to 0."
        exit 498
      }
      
    }   
     mat `weightvec'=`wmat'
    }
    else{
       if "`ort'"=="out"{   
          mat  `weightvec'=(J(1,`ninp',0),J(1,`ngo',1)*(1/2/`ngo'),J(1,`nbo',1)*(1/2/`nbo'))       
       }
       else{
          mat  `weightvec'=(J(1,`ninp',1/`ninp'),J(1,`ngo'+`nbo',0))         
       }
    
    }

  } 

*/
     
  if "`maxiter'"==""{
    local maxiter=-1
  }
  if "`tol'"==""{
    local tol=-1
  } 
  
   
    tempvar period dmu
  
  qui egen `period'=group(`time')
  qui egen `dmu'=group(`id')  

  
    qui su  `period'
    local tmax=r(max)

    tempvar flag temp DD D21 D12
    
    qui gen `DD'=.
    qui gen `D21'=.
    qui gen `D12'=.
    qui gen `temp'=.
    qui gen `flag'=0
  
    local gmat `gx' `gy' `gb'
    sort `period' `dmu'


  if `"`techtype'"'=="contemporaneous"{
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t')

            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
     
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'==`t'+1) 
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`D21') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
          
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1)
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')

        }       

    }
  
  
  }


  if `"`techtype'"'=="biennial"{
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t' | `period'==`t'+1)

             _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
      
        }    

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1 | `period'==`t')
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
    
        }       

    }
  
  
  }


  
  
  
    if `"`techtype'"'=="sequential"{
  
    
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t')
            _nddf if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
          
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'<=`t'+1) 
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`D21') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
          
        }  

        forv t=2/`tmax'{
            qui replace `flag'= (`period'<=`t'-1) 
            _nddf if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
           
        }       

    
  
  
  }
  
  
  
     if `"`techtype'"'=="window"{
    local band=`window'
   
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t'+`band' & `period'>=`t'-`band') 
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`DD') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'= (`period'<=`t'+1+`band' &  `period'>=`t'-`band'+1)
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`D21') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'<=`t'-1+`band' & `period'>=`t'-1-`band') 
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`D12') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
         
        }       

    
  
  
  }
 

  
 
  if `"`techtype'"'=="global"{

      qui replace `flag'=1
    _nddf  if `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
        
        qui bys `dmu' (`period'): gen TFPCH=`temp'[_n-1]-`temp' 
    label var TFPCH "Total factor productivity change"
    cap drop `temp'   
    
    sort `period' `dmu'
    forv t=1/`tmax'{
      qui replace `flag'=(`period'==`t')
      _nddf  if `touse' & `period'==`t', rflag(`flag') gen(`DD') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
   
    }
    
    qui bys `dmu' (`period'): gen TECH=`DD'[_n-1]-`DD'
    *qui bys `dmu' (`period'): gen BPC=TFPCH-TECH 
	qui bys `dmu' (`period'): gen TECCH=TFPCH-TECH
  
    label var TECH  "Technical efficiency change" 
    *label var BPC "Best practice gap change"
	label var TECCH "Technological Change"
    *local resvars TFPCH TECH  BPC
	local resvars TFPCH TECH  TECCH
    
    
  }
  else  if `"`techtype'"'=="biennial"{

    qui bys `dmu' (`period'): gen TFPCH=`DD'[_n-1]-`D12' 
    label var TFPCH "Total factor productivity change"
    cap drop `temp'   
    
    sort `period' `dmu'
    forv t=1/`tmax'{
      qui replace `flag'=(`period'==`t')
      _nddf  if `touse' & `period'==`t', rflag(`flag') gen(`DD') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
    
    }
    
    qui bys `dmu' (`period'): gen TECH=`DD'[_n-1]-`DD'
    qui bys `dmu' (`period'): gen TECCH=TFPCH-TECH      
  
    label var TECH  "Technical efficiency change" 
    label var TECCH "Techological change"
    local resvars TFPCH TECH  TECCH
    
    
  }

  else{
  
      //su `DD' `D12' `D21'
    qui {
      sort `dmu' `period'
      bys `dmu' (`period'): gen TECH=`DD'[_n-1]-`DD'
      bys `dmu' (`period'): gen TECCH=0.5*(`DD'-`D12'+`D21'[_n-1]-`DD'[_n-1])
      gen TFPCH= TECH+TECCH
      local resvars TFPCH TECH TECCH
        label var TFPCH "Total factor productivity change"
        label var TECH  "Technical efficiency change" 
        label var TECCH "Techological change"       
    } 
  
  }
  
  return local rvars "`resvars'"

end  

//////////////////////////////////////////////////////////////
capture program drop _nddf
program define _nddf
    version 16

    syntax [if] [in], gen(string) INvars(varlist) OPvars(varlist) ///
	                  BADvars(varlist) wmat(string) GVec(varlist) ///
					  [rflag(varname)  VRS maxiter(numlist) tol(numlist)]
        

    marksample touse 
    markout `touse' `invars' `opvars' `badvars' `gvec'
    
    tempvar touse2
     local rflag=cond("`rflag'"=="",""," if `rflag'")
    mark `touse2' `rflag' // `rflag' might be empty
    markout `touse2' `invars' `opvars' `badvars'
   

    local data `invars' `opvars' `badvars'
    local num1: word count `invars'   
    local num2: word count `opvars'
     
    if "`vrs'"!=""{
      local rts=1
    }
    else{
      local rts=0
    
    }
    
    mata: _nddf("`data'",`num1',`num2',"`touse'", "`touse2'","`gvec'","`wmat'","`gen'",`rts',`maxiter',`tol')
     

end 



*******************************************
///////////////////END//////////////////////
********************************************
