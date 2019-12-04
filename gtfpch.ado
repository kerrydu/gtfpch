*! version 1.0
* By Kerry Du, 3 Dec 2019 
**
* 
capture program drop gtfpch
program define gtfpch, rclass prop(xt)
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
	                                       SEQuential GLOBAL FGNZ RD  LUENberger ort(string) ///
										   SAVing(string)   Wmat(string) NONRadial       ///
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
                               `global' `sequential' window(`window') ///
                               maxiter(`maxiter') tol(`tol')
            local indexname "Malmquist-Luenberger Productivity Index"

        }
        else{

            _luen_ddf `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') gx(`gx') ///
                               gy(`gy') gb(`gb')  ort(`ort')   ///
                               `global' `sequential' window(`window') ///
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
                               `global' `sequential' window(`window') ///
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
                               `global' `sequential' window(`window') ///
                               maxiter(`maxiter') tol(`tol')
            //local indexname "Malmquist-Luenberger Productivity Index"

        }
        else{

            _luen_ddf `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') gx(`gx') ///
                               gy(`gy') gb(`gb')  ort(`ort')  vrs  ///
                               `global' `sequential' window(`window') ///
                               maxiter(`maxiter') tol(`tol')
            //local indexname "Luenberger Productivity Index (base on DDF)"
        }
    }
    else{

             _luen_nddf `invars'=`gopvars':`bopvars' if `touse', id(`id') time(`time') gx(`gx') ///
                               gy(`gy') gb(`gb')  ort(`ort')  wmat(`weightvec') ///
                               `global' `sequential' window(`window') vrs ///
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
                               gy(varlist) gb(varlist) VRS ort(string) ///
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
                               gy(varlist) gb(varlist) VRS ort(string) ///
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
                               gy(varlist) gb(varlist) VRS ort(string) ///
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




