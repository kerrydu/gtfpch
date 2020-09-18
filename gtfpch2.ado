*! version 2.0
* 30 Apr 2020
* add sbmdir option

*! version 1.1
* 30 Apr 2020
* add biennial option

*! version 1.0
* By Kerry Du, 3 Dec 2019 
**
* 
capture program drop gtfpch2
program define gtfpch2, rclass prop(xt)
    version 16

    qui mata mata mlib index
    _xt, trequired 
    local id=r(ivar)
    local time=r(tvar)

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
	                          SBMdir  BIennial  SEQuential GLOBAL FGNZ RD  LUENberger ort(string) ///
										   WINdow(numlist intege max=1 >=1) SAVing(string)   Wmat(string) NONRadial       ///
										   maxiter(numlist integer >0 max=1) tol(numlist max=1)]


    preserve
    marksample touse
    local bopvars `varlist'
    local ninp: word count `invars'
    local ngo: word count `gopvars'
    local nbo: word count `bopvars'
    local nvar: word count `invars' `gopvars' `bopvars'
    qui keep `invars' `gopvars' `bopvars' `touse' `id' `time' `dmu'
    qui gen Row=_n
    label var Row "Row # in the original dataset"   

    if `"`sbmdir'"'!=""{

        tfp_sbm `varlist' if `touse', dmu(`dmu') gx(`gx') gy(`gy') gb(`gb') `biennial' `sequential'  ///
                                  `global' `fgnz' `rd'  window(`window') wmat(`wmat') maxiter(`maxiter')  ///
                                  tol(`tol')  saving(`saving')

        exit


    }







    if "`nonradial'"!=""{
     local luenberger luenberger

     }

    if ("`ort'" == "") local ort = "out"
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
            di as err "(i|in|input|o|out|output) or nothing."
            exit 198
        }
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
  
    if  `"`ort'"'!="out"{
      local invarscopy `invars'
      forv k=1/`ninp'{
        gettoken word invarscopy:invarscopy
        tempvar gx_`k'
        qui gen `gx_`k''=-`word'
        local gmat `gmat' `gx_`k''
        local gmatname `gmatname' -`word'
      }   
    
    }
    else{
      forv k=1/`ninp'{
        tempvar gx_`k'
        qui gen `gx_`k''=0
        local gmat `gmat' `gx_`k''
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
  
    if  `"`ort'"'=="out"{
      local gopvarscopy `gopvars'
      forv k=1/`ngo'{
        gettoken word gopvarscopy:gopvarscopy
        tempvar gy_`k'
        qui gen `gy_`k''=`word'
        local gmat `gmat' `gy_`k''
        local gmatname `gmatname' `word'       
      }   
    
    }
    else{
      forv k=1/`ngo'{
        tempvar gy_`k'
        qui gen `gy_`k''=0
        local gmat `gmat' `gy_`k''
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
  
    if `"`ort'"'=="out"{
  
      local bopvarscopy `bopvars'
    forv k=1/`nbo'{
        gettoken word bopvarscopy:bopvarscopy
        tempvar gb_`k'
      qui gen `gb_`k''=-`word'
      local gmat `gmat' `gb_`k''
      local gmatname `gmatname' -`word'     
    } 
  
    }
    else{
      forv k=1/`nbo'{
        tempvar gb_`k'
        qui gen `gb_`k''=0
        local gmat `gmat' `gb_`k''
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
            local indexname "Luenberger Productivity Index (base on DDF)"
        }
    }
    else{

            matrix colnames `weightvec' =`invars' `gopvars' `bopvars' 
            matrix rownames `weightvec' = "W"
            di 
            disp " The weight vector is:"
            mat list `weightvec',noblank  noheader

             _luen_nddf `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') gx(`gx') ///
                               gy(`gy') gb(`gb')  ort(`ort')  wmat(`weightvec') ///
                              `biennial' `global' `sequential' window(`window') ///
                               maxiter(`maxiter') tol(`tol')       
            local indexname "Luenberger Productivity Index (base on nonrial DDF)"

    }
            local resvars `r(rvars)'

  if "`rd'"==""&"`fgnz'"==""{

        di 
        di " The diectional vector is (`gmatname')"
        format `resvars' %9.4f
        qui keep if  `touse'
        qui cap bys `id' (`time'): gen Pdwise=`time'[_n-1]+"~"+`time' if _n>1
        qui cap bys `id' (`time'): gen Pdwise=string(`time'[_n-1])+"~"+string(`time') if _n>1
        label var Pdwise "Period wise"       
        order Row `dmu' `id' Pdwise  `resvars' 
        qui keep  Row `dmu' `id' Pdwise  `resvars'    
        qui keep if !missing(Pdwise)
        disp _n(2) " Total Factor Productivity Change:`indexname'"
        disp "    (Row: Row # in the original data; Pdwise: periodwise)"

        list Row `dmu' `id'  Pdwise  `resvars', sep(0) 
        di "Note: missing value indicates infeasible problem."

        if `"`saving'"'!=""{
          save `saving'
          gettoken filenames saving:saving, parse(",")
          local filenames `filenames'.dta
          disp _n `"Estimated Results are saved in `filenames'."'
        }   
        

        return local file `filenames'       

        restore

  }
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
            //local indexname "Luenberger Productivity Index (base on nonrial DDF)"

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
        if "`global'"!=""{
          qui replace BPC=BPC_crs
        }
        else{
          qui replace TECCH=TECCH_crs   
        }

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
        if "`global'"!=""{
          qui replace BPC=BPC_crs
        }
        else{
          qui replace TECCH=TECCH_crs   
        }

        local resvars `resvars' SECH
      }

    }

    
        di 
        di " The diectional vector is (`gmatname')"    
      format `resvars' %9.4f

    qui keep if `touse'
    qui cap bys `id' (`time'): gen Pdwise=`time'[_n-1]+"~"+`time' if _n>1
    qui cap bys `id' (`time'): gen Pdwise=string(`time'[_n-1])+"~"+string(`time') if _n>1
    label var Pdwise "Period wise"
          
    order Row `dmu' `id' Pdwise  `resvars' 
    qui keep if !missing(Pdwise) & `touse'
    qui keep  Row `dmu' `id' Pdwise  `resvars' 
  
    disp _n(2) " Total Factor Productivity Change:`indexname'"
    disp "    (Row: Row # in the original data; Pdwise: periodwise)"

    list Row `dmu' `id'  Pdwise  `resvars' , sep(0) 
    di "Note: missing value indicates infeasible problem."

    if `"`saving'"'!=""{
      save `saving'
      gettoken filenames saving:saving, parse(",")
      local filenames `filenames'.dta
      disp _n `"Estimated Results are saved in `filenames'."'
    } 
    

      return local file `filenames'
      restore     



  }
        



end

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

********************************************************************************* 
  
    *local ninv: word count `invars'
    *local ngov: word count `gopvars'
    syntax varlist [if] [in], id(varname) time(varname)  [gx(varlist) ///
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
     



  if "`gx'"!=""{
    local ngx: word count `gx'
    if `ngx'!=`ninp'{
        disp as error "# of input variables != # of variables specified in gx()."
        error 498
    }
    local gmat `gmat' `gx'
  
  }
  else{
  
    if  `"`ort'"'!="out"{
      local invarscopy `invars'
      forv k=1/`ninp'{
        gettoken word invarscopy:invarscopy
        tempvar gx_`k'
        qui gen `gx_`k''=-`word'
        local gmat `gmat' `gx_`k''
      }   
    
    }
    else{
      forv k=1/`ninp'{
        tempvar gx_`k'
        qui gen `gx_`k''=0
        local gmat `gmat' `gx_`k''
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
  
  }
  else{
  
    if  `"`ort'"'=="out"{
      local gopvarscopy `gopvars'
      forv k=1/`ngo'{
        gettoken word gopvarscopy:gopvarscopy
        tempvar gy_`k'
        qui gen `gy_`k''=`word'
        local gmat `gmat' `gy_`k''
      }   
    
    }
    else{
      forv k=1/`ngo'{
        tempvar gy_`k'
        qui gen `gy_`k''=0
        local gmat `gmat' `gy_`k''
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
  
  }
  else{
  
    if `"`ort'"'=="out"{
  
      local bopvarscopy `bopvars'
    forv k=1/`nbo'{
        gettoken word bopvarscopy:bopvarscopy
        tempvar gb_`k'
      qui gen `gb_`k''=-`word'
      local gmat `gmat' `gb_`k''
    } 
  
    }
    else{
      forv k=1/`nbo'{
        tempvar gb_`k'
        qui gen `gb_`k''=0
        local gmat `gmat' `gb_`k''
      }     
    
    
    }

  
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
  
    qui gen `flag'=0
  
  sort `period' `dmu'
  
  if `"`techtype'"'=="contemporaneous"{
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t')
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui drop `temp'
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'==`t'+1) 
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D21'=`temp' if `period'==`t'
            qui drop `temp'
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1)
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui drop `temp'
        }       

    }
  
  
  }
  
  
    if `"`techtype'"'=="sequential"{
  
    
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t')
            _ddf if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui replace `flag'=0
            qui drop `temp'
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'<=`t'+1) 
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D21'=`temp' if `period'==`t'
            qui drop `temp'
        }  

        forv t=2/`tmax'{
            qui replace `flag'= (`period'<=`t'-1) 
            _ddf if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui drop `temp'
        }       

    
  
  
  }
  
  
  
     if `"`techtype'"'=="window"{
    local band=(`window'-1)/2
   
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t'+`band' & `period'>=`t'-`band') 
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui cap drop `temp'
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'= (`period'<=`t'+1+`band' &  `period'>=`t'-`band'+1)
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D21'=`temp' if `period'==`t'
            qui cap drop `temp'
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'<=`t'-1+`band' & `period'>=`t'-1-`band') 
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui cap drop `temp'
        }       

    
  
  
  }



  if `"`techtype'"'=="biennial" {
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t' | `period'==`t'+1)
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui drop `temp'
        }     

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1 | `period'==`t')
            _ddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui drop `temp'
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
      _ddf  if `touse' & `period'==`t', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
      qui replace `DD'=`temp' if `period'==`t'
      qui cap drop `temp'
    }
    
    qui bys `dmu' (`period'): gen TECH=`DD'/`DD'[_n-1]  
    qui bys `dmu' (`period'): gen BPC=TFPCH/TECH      
  
    label var TECH  "Technical efficiency change" 
    label var BPC "Best practice gap change"
    local resvars TFPCH TECH  BPC
    
    
  }

  else if `"`techtype'"'=="biennial"{

    qui bys `dmu' (`period'): gen TFPCH=`D12'/`DD'[_n-1] if _n>1
    label var TFPCH "Total factor productivity change"
    cap drop `temp' 

    sort `period' `dmu'
    forv t=1/`tmax'{
      qui replace `flag'=(`period'==`t')
      _ddf  if `touse' & `period'==`t', rflag(`flag') gen(`temp') gv(`gmat') `vrs' ort(`ort') in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
      qui replace `DD'=`temp' if `period'==`t'
      qui cap drop `temp'
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
    mark `touse2' if `rflag' // `rflag' might be empty
    markout `touse2' `invars' `opvars' `badvars'
    //qui gen `touse2'=`rflag'  
        qui gen `gen'=.

        local data `invars' `opvars' `badvars'
        local num1: word count `invars'   
        local num2: word count `opvars'
    
*******************************************************************************
/////////This section is from  Yong-bae Ji and Choonjoo Lee's DEA.ado//////////   
    // default orientation - Input Oriented
    if ("`ort'" == "") local ort = "IN"
    else {
      local ort = upper("`ort'")
      if ("`ort'" == "I" | "`ort'" == "IN" | "`ort'" == "INPUT") {
        local ort = "IN"
      }
      else if ("`ort'" == "O" | "`ort'" == "OUT" | "`ort'" == "OUTPUT") {
        local ort = "OUT"
      }
      else {
        di as err "option ort allows for case-insensitive " _c
        di as err "(i|in|input|o|out|output) or nothing."
        exit 198
      }
    }
    
*******************************************************************************   
    if "`vrs'"!=""{
      local rts=1
    }
    else{
      local rts=0
    
    }
    
    mata: _ddf("`data'",`num1',`num2',"`touse'", "`touse2'","`gvec'","`gen'",`rts',`maxiter',`tol')

    if "`ort'" =="OUT"{
      
      qui replace `gen'=1/(1+`gen')
      }
    else{
      qui replace `gen'=(1-`gen')
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


  if "`gx'"!=""{
    local ngx: word count `gx'
    if `ngx'!=`ninp'{
        disp as error "# of input variables != # of variables specified in gx()."
        error 498
    }
    local gmat `gmat' `gx'
  
  }
  else{
  
    if  `"`ort'"'!="out"{
      local invarscopy `invars'
      forv k=1/`ninp'{
        gettoken word invarscopy:invarscopy
        tempvar gx_`k'
        qui gen `gx_`k''=-`word'
        local gmat `gmat' `gx_`k''
      }   
    
    }
    else{
      forv k=1/`ninp'{
        tempvar gx_`k'
        qui gen `gx_`k''=0
        local gmat `gmat' `gx_`k''
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
  
  }
  else{
  
    if  `"`ort'"'=="out"{
      local gopvarscopy `gopvars'
      forv k=1/`ngo'{
        gettoken word gopvarscopy:gopvarscopy
        tempvar gy_`k'
        qui gen `gy_`k''=`word'
        local gmat `gmat' `gy_`k''
      }   
    
    }
    else{
      forv k=1/`ngo'{
        tempvar gy_`k'
        qui gen `gy_`k''=0
        local gmat `gmat' `gy_`k''
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
  
  }
  else{
  
    if `"`ort'"'=="out"{
  
      local bopvarscopy `bopvars'
    forv k=1/`nbo'{
        gettoken word bopvarscopy:bopvarscopy
        tempvar gb_`k'
      qui gen `gb_`k''=-`word'
      local gmat `gmat' `gb_`k''
    } 
  
    }
    else{
      forv k=1/`nbo'{
        tempvar gb_`k'
        qui gen `gb_`k''=0
        local gmat `gmat' `gb_`k''
      }     
    
    
    }

  
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
  
    qui gen `flag'=0
  
  sort `period' `dmu'
  
  if `"`techtype'"'=="contemporaneous"{
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t')
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui drop `temp'
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'==`t'+1) 
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D21'=`temp' if `period'==`t'
            qui drop `temp'
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1)
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui drop `temp'
        }       

    }
  
  
  }


  if `"`techtype'"'=="biennial"{
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t' | `period'==`t'+1)
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui drop `temp'
        }    

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1 | `period'==`t')
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui drop `temp'
        }       

    }
  
  
  }

  
  
    if `"`techtype'"'=="sequential"{
  
    
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t')
            _ddf_luen if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui replace `flag'=0
            qui drop `temp'
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'<=`t'+1) 
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D21'=`temp' if `period'==`t'
            qui drop `temp'
        }  

        forv t=2/`tmax'{
            qui replace `flag'= (`period'<=`t'-1) 
            _ddf_luen if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui drop `temp'
        }       

    
  
  
  }
  
  
  
     if `"`techtype'"'=="window"{
    local band=(`window'-1)/2
   
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t'+`band' & `period'>=`t'-`band') 
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui cap drop `temp'
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'= (`period'<=`t'+1+`band' &  `period'>=`t'-`band'+1)
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D21'=`temp' if `period'==`t'
            qui cap drop `temp'
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'<=`t'-1+`band' & `period'>=`t'-1-`band') 
            _ddf_luen  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui cap drop `temp'
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
      _ddf_luen  if `touse' & `period'==`t', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
      qui replace `DD'=`temp' if `period'==`t'
      qui cap drop `temp'
    }
    
    qui bys `dmu' (`period'): gen TECH=`DD'[_n-1]-`DD'
    qui bys `dmu' (`period'): gen BPC=TFPCH-TECH      
  
    label var TECH  "Technical efficiency change" 
    label var BPC "Best practice gap change"
    local resvars TFPCH TECH  BPC
    
    
  }
  else if `"`techtype'"'="biennial"{

     qui bys `dmu' (`period'): gen TFPCH=`DD'[_n-1]-`D12' if _n>1
     label var TFPCH "Total factor productivity change"
     qui cap drop `temp'
     sort `period' `dmu'
    forv t=1/`tmax'{
      qui replace `flag'=(`period'==`t')
      _ddf_luen  if `touse' & `period'==`t', rflag(`flag') gen(`temp') gv(`gmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
      qui replace `DD'=`temp' if `period'==`t'
      qui cap drop `temp'
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
    
    tempvar touse2
    mark `touse2' if `rflag' // `rflag' might be empty
    markout `touse2' `invars' `opvars' `badvars'
    //qui gen `touse2'=`rflag'  
        qui gen `gen'=.

        local data `invars' `opvars' `badvars'
        local num1: word count `invars'   
        local num2: word count `opvars'
    
*******************************************************************************
    
*******************************************************************************   
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

********************************************************************************* 
  
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
     


  if "`gx'"!=""{
    local ngx: word count `gx'
    if `ngx'!=`ninp'{
        disp as error "# of input variables != # of variables specified in gx()."
        error 498
    }
    local gmat `gmat' `gx'
  
  }
  else{
  
    if  `"`ort'"'!="out"{
      local invarscopy `invars'
      forv k=1/`ninp'{
        gettoken word invarscopy:invarscopy
        tempvar gx_`k'
        qui gen `gx_`k''=-`word'
        local gmat `gmat' `gx_`k''
      }   
    
    }
    else{
      forv k=1/`ninp'{
        tempvar gx_`k'
        qui gen `gx_`k''=0
        local gmat `gmat' `gx_`k''
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
  
  }
  else{
  
    if  `"`ort'"'=="out"{
      local gopvarscopy `gopvars'
      forv k=1/`ngo'{
        gettoken word gopvarscopy:gopvarscopy
        tempvar gy_`k'
        qui gen `gy_`k''=`word'
        local gmat `gmat' `gy_`k''
      }   
    
    }
    else{
      forv k=1/`ngo'{
        tempvar gy_`k'
        qui gen `gy_`k''=0
        local gmat `gmat' `gy_`k''
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
  
  }
  else{
  
    if `"`ort'"'=="out"{
  
      local bopvarscopy `bopvars'
    forv k=1/`nbo'{
        gettoken word bopvarscopy:bopvarscopy
        tempvar gb_`k'
      qui gen `gb_`k''=-`word'
      local gmat `gmat' `gb_`k''
    } 
  
    }
    else{
      forv k=1/`nbo'{
        tempvar gb_`k'
        qui gen `gb_`k''=0
        local gmat `gmat' `gb_`k''
      }     
    
    
    }

  
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
  
    qui gen `flag'=0
  
  sort `period' `dmu'


  if `"`techtype'"'=="contemporaneous"{
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t')

            noi _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui drop `temp'
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'==`t'+1) 
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D21'=`temp' if `period'==`t'
            qui drop `temp'
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1)
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui drop `temp'
        }       

    }
  
  
  }


  if `"`techtype'"'=="biennial"{
  
      qui{
        forv t=1/`tmax'{
            qui replace `flag'= (`period'==`t' | `period'==`t'+1)

             _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui drop `temp'
        }    

        forv t=2/`tmax'{
            qui replace `flag'=(`period'==`t'-1 | `period'==`t')
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui drop `temp'
        }       

    }
  
  
  }


  
  
  
    if `"`techtype'"'=="sequential"{
  
    
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t')
            _nddf if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui replace `flag'=0
            qui drop `temp'
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'=(`period'<=`t'+1) 
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D21'=`temp' if `period'==`t'
            qui drop `temp'
        }  

        forv t=2/`tmax'{
            qui replace `flag'= (`period'<=`t'-1) 
            _nddf if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui drop `temp'
        }       

    
  
  
  }
  
  
  
     if `"`techtype'"'=="window"{
    local band=(`window'-1)/2
   
        forv t=1/`tmax'{
            qui replace `flag'=(`period'<=`t'+`band' & `period'>=`t'-`band') 
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui cap drop `temp'
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
            qui replace `flag'= (`period'<=`t'+1+`band' &  `period'>=`t'-`band'+1)
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D21'=`temp' if `period'==`t'
            qui cap drop `temp'
        }  

        forv t=2/`tmax'{
            qui replace `flag'=(`period'<=`t'-1+`band' & `period'>=`t'-1-`band') 
            _nddf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui cap drop `temp'
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
      _nddf  if `touse' & `period'==`t', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
      qui replace `DD'=`temp' if `period'==`t'
      qui cap drop `temp'
    }
    
    qui bys `dmu' (`period'): gen TECH=`DD'[_n-1]-`DD'
    qui bys `dmu' (`period'): gen BPC=TFPCH-TECH      
  
    label var TECH  "Technical efficiency change" 
    label var BPC "Best practice gap change"
    local resvars TFPCH TECH  BPC
    
    
  }
  else  if `"`techtype'"'=="biennial"{

    qui bys `dmu' (`period'): gen TFPCH=`DD'[_n-1]-`D12' 
    label var TFPCH "Total factor productivity change"
    cap drop `temp'   
    
    sort `period' `dmu'
    forv t=1/`tmax'{
      qui replace `flag'=(`period'==`t')
      _nddf  if `touse' & `period'==`t', rflag(`flag') gen(`temp') gv(`gmat') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
      qui replace `DD'=`temp' if `period'==`t'
      qui cap drop `temp'
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

    syntax [if] [in], gen(string) INvars(varlist) OPvars(varlist) BADvars(varlist) wmat(string) GVec(varlist) [rflag(varname)  VRS maxiter(numlist) tol(numlist)]
        

        marksample touse 
    markout `touse' `invars' `opvars' `badvars' `gvec'
    
    tempvar touse2
    mark `touse2' if `rflag' // `rflag' might be empty
    markout `touse2' `invars' `opvars' `badvars'
    //qui gen `touse2'=`rflag'  
        qui gen `gen'=.

        local data `invars' `opvars' `badvars'
        local num1: word count `invars'   
        local num2: word count `opvars'
    


    
*******************************************************************************   
    if "`vrs'"!=""{
      local rts=1
    }
    else{
      local rts=0
    
    }
    
    mata: _nddf("`data'",`num1',`num2',"`touse'", "`touse2'","`gvec'","`wmat'","`gen'",`rts',`maxiter',`tol')
     

end 





********************************************




* 
capture program drop tfp_sbm
program define tfp_sbm, rclass prop(xt)
    version 16
    
    qui mata mata mlib index
  
    _xt, trequired 
    local id=r(ivar)
    local time=r(tvar)

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
  
  
    syntax varlist [if] [in],    [dmu(varname)   gx(varlist) gy(varlist) gb(varlist) ///
                                   BIennial      SEQuential GLOBAL ///
                       WINdow(numlist intege max=1 >=1) SAVing(string)   Wmat(string)     ///
                       maxiter(numlist integer >0 max=1) tol(numlist max=1)]
    
  if `"`wmat'"'!=""  confirm matrix `wmat'
  
  
    preserve
    marksample touse

    local bopvars `varlist'
    local ninp: word count `invars'
    local ngo: word count `gopvars'
    local nbo: word count `bopvars'
    local nvar: word count `invars' `gopvars' `bopvars'
    qui keep `invars' `gopvars' `bopvars' `touse' `id' `time' `dmu'  `gx' `gy' `gb'

    qui gen Row=_n
    label var Row "Row # in the original dataset"
  


  if "`gx'"!=""{
    local ngx: word count `gx'
    if `ngx'!=`ninp'{
        disp as error "# of input variables != # of variables specified in gx()."
        error 498
    }
    local gmat `gmat' `gx'
  
  }
  else{
  
      local invarscopy `invars'
      forv k=1/`ninp'{
        gettoken word invarscopy:invarscopy
        tempvar gx_`k'
        qui gen `gx_`k''=-`word'
        local gx `gx' `gx_`k''
      }
      local gmat `gmat' `gx'   
    
    }
  

  
  }
  
  if "`gy'"!=""{
    local ngy: word count `gy'
    if `ngy'!=`ngo'{
        disp as error "# of desriable output variables != # of variables specified in gy()."
        error 498
    }
    local gmat `gmat' `gy'
  
  }
  else{
  
      local gopvarscopy `gopvars'
      forv k=1/`ngo'{
        gettoken word gopvarscopy:gopvarscopy
        tempvar gy_`k'
        qui gen `gy_`k''=`word'
        local gy `gy' `gy_`k''
      }  
       local gmat `gmat' `gy'
    
  }
    
  if "`gb'"!=""{
    local ngb: word count `gb'
    if `ngb'!=`nbo'{
        disp as error "# of undesriable output variables != # of variables specified in gb()."
        error 498
    }
    local gmat `gmat' `gb'
  
  }
  else{
  
  
    local bopvarscopy `bopvars'
    forv k=1/`nbo'{
        gettoken word bopvarscopy:bopvarscopy
        tempvar gb_`k'
      qui gen `gb_`k''=-`word'
      local gb `gb' `gb_`k''
    } 
    local gmat `gmat' `gb'
  
  
  } 



  
    ltfp_sbmdf `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') ///
                               `global' `sequential' window(`window') ///
                               maxiter(`maxiter') tol(`tol')   wmat(`wmat') ///
                               gx(`gx') gy(`gy') gb(`gb')
  
    local resvars `r(rvars)'
      
      
        di 
        format `resvars' %9.4f
        qui keep if  `touse'
        qui cap bys `id' (`time'): gen Pdwise=`time'[_n-1]+"~"+`time' if _n>1
        qui cap bys `id' (`time'): gen Pdwise=string(`time'[_n-1])+"~"+string(`time') if _n>1
        label var Pdwise "Period wise"       
        order Row `dmu' `id' Pdwise  `resvars' 
        qui keep  Row `dmu' `id' Pdwise  `resvars'    
        qui keep if !missing(Pdwise)
        disp _n(2) " Luenberger Productivity Inducator:"
        disp "    (Row: Row # in the original data; Pdwise: periodwise)"

        list Row `dmu' `id'  Pdwise  `resvars', sep(0) 
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



cap program drop ltfp_sbmdf
program define ltfp_sbmdf,rclass
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
    syntax varlist [if] [in], id(varname) time(varname) gx(varlist) gy(varlist) gb(varlist)  [wmat(string)  ///
                            BIennial GLOBAL SEQuential WINdow(numlist intege max=1 >=1) ///
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
     
     local techtype "global"
  
  } 
  
   

   if "`sequential'"!=""{
 
     if "`window'"!=""{
     
       disp as error "sequential and window() cannot be specified together."
       error 498     
     
     }     
     
     local techtype "sequential"
  
  } 
    
 
     if "`window'"!=""{
     
         local techtype "window"   
     
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

    tempvar flag temp DD D21 D12  DDv D21v D12v
    
    qui gen `DD'=.
    qui gen `D21'=.
    qui gen `D12'=.
  
    qui gen `DDv'=.
    qui gen `D21v'=.
    qui gen `D12v'=.  
  
  
  
    qui gen `flag'=0
  
   sort `period' `dmu'
  
 
  if `"`techtype'"'=="global"{

      qui replace `flag'=1
      sbmdf if `touse', rflag(`flag') gen(`temp')  wmat(`wmat')   in(`invars') ///
                     op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol') ///
                     gx(`gx') gy(`gy') gb(`gb')
        
    qui bys `dmu' (`period'): gen TFPCH=`temp'[_n-1]-`temp' 
    label var TFPCH "Total factor productivity change"
    cap drop `temp'   
    
    sort `period' `dmu'
    forv t=1/`tmax'{
      qui replace `flag'=(`period'==`t')
      sbmdf  if `touse' & `period'==`t', rflag(`flag') gen(`temp') ///
       wmat(`wmat')  in(`invars') op(`gopvars') bad(`bopvars') ///
        maxiter(`maxiter') tol(`tol') gx(`gx') gy(`gy') gb(`gb')
      qui replace `DD'=`temp' if `period'==`t'
      qui cap drop `temp'
    }
    
    qui bys `dmu' (`period'): gen TECH=`DD'[_n-1]-`DD'
    qui bys `dmu' (`period'): gen BPC=TFPCH-TECH      
  
    label var TECH  "Technical efficiency change" 
    label var BPC "Best practice gap change"
    local resvars TFPCH TECH  BPC
    
    
  }
  else{
  
  
  
  
        forv t=1/`tmax'{
      
      
      if `"`techtype'"'=="contemporaneous"  local contt   (`period'<`t')
        
      
      if `"`techtype'"'=="sequential"   local contt   (`period'<=`t')
  
        
      
       if `"`techtype'"'=="window"{
        local band=(`window'-1)/2
        local contt   ( `period'<=`t'+`band' & `period'>=`t'-`band')      
        
       }      
          
          
      
            qui replace `flag'=`contt' 
            sbmdf  if `period'==`t' & `touse', rflag(`flag') gen(`temp')  wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `DD'=`temp' if `period'==`t'
            qui cap drop `temp'
            sbmdf  if `period'==`t' & `touse', rflag(`flag') gen(`temp')  wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol') vrs
            qui replace `DDv'=`temp' if `period'==`t'
            qui cap drop `temp'     
      
        }    
        local tt=`tmax'-1
        forv t=1/`tt'{
      
      if `"`techtype'"'=="contemporaneous" local con21  ( `period'<`t'+1)
      
      
      if `"`techtype'"'=="sequential"   local con21  ( `period'<=`t'+1)   
      
        
      
       if `"`techtype'"'=="window"{
        local band=(`window'-1)/2
        local con21   (`period'<=`t'+1+`band' &  `period'>=`t'-`band'+1)    
        
       }      
                
      
            qui replace `flag'= `con21'
            sbmdf if `period'==`t' & `touse', rflag(`flag') gen(`temp') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D21'=`temp' if `period'==`t'
            qui cap drop `temp'
            sbmdf if `period'==`t' & `touse', rflag(`flag') gen(`temp') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol') vrs
            qui replace `D21v'=`temp' if `period'==`t'
            qui cap drop `temp'     
      
        }  

        forv t=2/`tmax'{
      
      
      if `"`techtype'"'=="contemporaneous"  local con12   (`period'<`t'-1)    
      if `"`techtype'"'=="sequential"   local con12   (`period'<=`t'-1)   
        
      
       if `"`techtype'"'=="window"{
        local band=(`window'-1)/2
        local con12   (`period'<=`t'-1+`band' & `period'>=`t'-1-`band')       
        
       }      
                
      
            qui replace `flag'=`con12'
            sbmdf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol')
            qui replace `D12'=`temp' if `period'==`t'
            qui cap drop `temp'
            sbmdf  if `period'==`t' & `touse', rflag(`flag') gen(`temp') wmat(`wmat') `vrs'  in(`invars') op(`gopvars') bad(`bopvars') maxiter(`maxiter') tol(`tol') vrs
            qui replace `D12v'=`temp' if `period'==`t'
            qui cap drop `temp'     
      
        }  
  
  
      //su `DD' `D12' `D21'
    qui {
      sort `dmu' `period'
      bys `dmu' (`period'): gen LPEC=`DDv'[_n-1]-`DDv'
      bys `dmu' (`period'): gen LTPT=0.5*(`DDv'-`D12v'+`D21v'[_n-1]-`DDv'[_n-1])
    bys `dmu' (`period'): gen LSEC=`DD'[_n-1]-`DD'-(`DDv'[_n-1]-`DDv')
    bys `dmu' (`period'): gen LTPSC=0.5*(`DD'-`D12'+`D21'[_n-1]-`DD'[_n-1])-LTPT
      gen LTFP= 0.5*(`DD'-`D12'+`D21'[_n-1]-`DD'[_n-1])+`DD'[_n-1]-`DD'
      local resvars LTFP  LPEC LTPT LSEC LTPSC
        //label var TFPCH "Total factor productivity change"
        //label var TECH  "Technical efficiency change" 
        //label var TECCH "Techological change"       
    } 
  
  }
  
  return local rvars "`resvars'"

end  



capture program drop sbmdf
program define sbmdf
    version 16

    syntax [if] [in], gen(string) gx(varlist) gy(varlist) gb(varlist) INvars(varlist) OPvars(varlist) BADvars(varlist)   [Slacks(string) rflag(varname) wmat(string) VRS maxiter(numlist) tol(numlist)]
        

    marksample touse 
    markout `touse' `invars' `opvars' `badvars' 
    
  tempvar touse2
  //qui gen byte `touse2'=1
  
  if "`rflag'"!=""{
      
        //qui replace `touse2'= `rflag' // `rflag' might be empty 
    mark `touse2' if  `rflag' & !missing(`rflag')
    
  }
  else{
    mark `touse2'
  }

    markout `touse2' `invars' `opvars' `badvars'
    //qui gen `touse2'=`rflag'  
  confirm new var `gen'
        qui gen `gen'=.

        local data `invars' `opvars' `badvars'
        local num1: word count `invars'   
        local num2: word count `opvars'
    
  if `"`slacks'"'!=""{
    foreach v in `invars' `opvars' `badvars'{
      qui gen `slacks'_`v'=.
      label var `slacks'_`v' `"slack in `v'"'
      local sname `sname' `slacks'_`v'
    }
  }

    
    if "`vrs'"!=""{
      local rts=1
    }
    else{
      local rts=0
    
    }
    
    mata: _sbmdf("`data'",`num1',`num2',"`gx'","`gy'","`gb'","`touse'", "`touse2'","`wmat'","`gen'","`sname'",`rts',-1,-1)
     
   
    if `"`slacks'"'!=""{
    foreach v in `invars' `opvars' `badvars'{
      qui replace `slacks'_`v'=`slacks'_`v'*`v'
    }
  }

end 


cap mata mata drop _sbmdf()
cap mata mata drop _sbmdf2()
cap mata mata drop lpres()

mata:

  struct lpres   { real scalar    fval 
           real matrix    coeff
           real scalar    converged
           real scalar    returncode
            }


void function _sbmdf(   string scalar     d, ///
                      real scalar       nx, ///
                      real scalar       ny, ///
                      string scalar   gx0,  ///
                      string scalar   gy0, ///
                      string scalar   gb0, ///
            string scalar     touse1,  ///
            string scalar     touse2, ///
            string scalar     wname, ///
            string scalar     bname, ///
            string scalar     sname, ///            
            real scalar       rstype, ///         
            real scalar       maxiter, ///
            real scalar       tol)
  {
    
     
     struct lpres scalar sbmres
    

          //wmat=st_matrix(wname)
          data=st_data(.,d,touse1)
          data=data' 
          gx=st_data(.,gx0,touse1)
          gx=gx'
          gy=st_data(.,gy0,touse1)
          gy=gy'
          gb=st_data(.,gb0,touse1)
          gb=gb'         

          M=rows(data)
          dataref=st_data(.,d,touse2)
          dataref=dataref'
          Xref=dataref[1..nx,.]
          Yref=dataref[nx+1..(nx+ny),.]
          Bref=dataref[(nx+ny+1)..M,.]
          X=data[1..nx,.]
          Y=data[nx+1..(nx+ny),.]
          B=data[(nx+ny+1)..M,.]
          //=cols(dataref)  
      //wmat=(0.5:/J(1,nx,nx),0.5:/J(1,M-nx,M-nx))
          theta=J(cols(data),1,.)
      slacks=J(cols(data),M,.)
      
      if(wname!=""){
        wmat=st_matrix(wname)
      }
      else{
        wmat=(0.5:/J(1,nx,nx),0.5:/J(1,M-nx,M-nx))
      }     
      
          for(j=1;j<=cols(data);j++){

         sbmres=_sbmdf2(X[.,j],Y[.,j],B[.,j],gx[.,j],gy[.,j],gb[.,j],Xref,Yref,Bref,wmat,rstype,maxiter,tol)    
           if(sbmres.converged==1){
             theta[j]=sbmres.fval 
             slacks[j,.]=sbmres.coeff[1,1..M]      
           }         
         
         
                 //q.getInequality()
         }

          st_view(gen=.,.,bname,touse1)
          gen[.,.]=theta
    
     if(sname!=""){
        st_view(SLMAT=.,.,sname,touse1)
            SLMAT[.,.]=slacks 
     }
    
      
}



  function _sbmdf2(   real colvector     X, ///
                      real colvector     Y, ///
                      real colvector     B, ///
                      real colvector     gX, ///
                      real colvector     gY, ///
                      real colvector     gB, ///                     
            real matrix        Xref,  ///
            real matrix        Yref,  ///
            real matrix        Bref,  ///
            real rowvector     wmat, ///
            real scalar       rstype, ///         
            real scalar       maxiter, ///
            real scalar       tol)
  {

          
          nx=length(X)
      ny=length(Y)
      nb=length(B)
      M=nb+nx+ny
          //gX=-X
          //gY=Y
          //gB=-B     
          N=cols(Xref)
      class LinearProgram scalar q
      q = LinearProgram()
    c=(wmat,J(1,N,0))
    lowerbd=J(1,length(c),0)
    upperbd=J(1,length(c),.)    
    q.setCoefficients(c)
    q.setBounds(lowerbd, upperbd)

    if(maxiter!=-1){
      q.setMaxiter(maxiter)
    }
    if (tol!=-1){
      q.setTol(tol)
    }     

          
  

         Aie1=diag(-gX),J(nx,ny+nb,0),Xref
         bie1=X
         Aie2=J(ny,nx,0),diag(gY),J(ny,nb,0),-Yref
         bie2=-Y
         Aec=J(nb,nx+ny,0),diag(-gB),Bref
         bec=B

        if(rstype==0){
          q.setEquality((Aie1 \ Aie2 \ Aec),(bie1 \ bie2 \ bec))      
        }
        else{
           lsum=J(1,M,0),J(1,N,1)
           q.setEquality((Aie1 \ Aie2 \ Aec \ lsum), (bie1 \ bie2 \ bec \ 1))     
         }

                 beta=q.optimize()  
         //q.optimize()
         struct lpres scalar retres
         retres.fval=q.value()
         retres.coeff=q.parameters()
         retres.converged=q.converged()
         retres.returncode=q.returncode()
         return(retres)
        

}




end

