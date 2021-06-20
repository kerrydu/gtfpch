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
										   maxiter(numlist integer >0 max=1) tol(numlist max=1) NODOTS SAVing(string) noPRINT]
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
	
	
///////////////////////////////////////////////////////	


/*
version 16
cap mata mata drop _ddf()
cap mata mata drop _nddf()
cap mata mata drop _ddf2()
cap mata mata drop _nddf2()
cap mata mata drop _tenddf()
cap mata mata drop _tenddf2()
cap mata mata drop ddfcontemporaneous()
cap mata mata drop ddfcontemporaneousbt()
cap mata mata drop ddfbiennial()
cap mata mata drop ddfbiennialbt()
cap mata mata drop ddfsequential()
cap mata mata drop ddfsequentialbt()
cap mata mata drop ddfwindow()
cap mata mata drop ddfwindowbt()
cap mata mata drop ddfglobal()
cap mata mata drop ddfglobalbt()
cap mata mata drop nddfcontemporaneous()
cap mata mata drop nddfcontemporaneousbt()
cap mata mata drop nddfbiennial()
cap mata mata drop nddfbiennialbt()
cap mata mata drop nddfsequential()
cap mata mata drop nddfsequentialbt()
cap mata mata drop nddfwindow()
cap mata mata drop nddfwindowbt()
cap mata mata drop nddfglobal()
cap mata mata drop nddfglobalbt()
cap mata mata drop biasandciddf()
cap mata mata drop displaydots()
cap mata mata drop calddf()
cap mata mata drop calnddf()
cap mata mata drop bootdea()

mata:
mata clear 

void function displaydots(real scalar k)
{
	if(mod(k,50)==0){
		printf("..........%2.0f\n",k)
	}
	//else{
	//	printf(".")
	//}
}



// define structure
	struct bootdea { 
					   real colvector tebc
					   real colvector LL
					   real colvector UU
					   real colvector vari	
				       real colvector BovV			   
				       real colvector bias					   				   
					
					  }
					  

real matrix function  calddf(real matrix data, 
                             real matrix dataref, 
                             real scalar  M,
							 real scalar  N,
							 real matrix gmat,
					         real scalar rstype, 					
					         real scalar maxiter, 
					         real scalar  tol    )
{
    ///nb=rows(data)-M-N
	X=data[1::M,.]
    Y=data[(M+1)::(M+N),.]	
	B=data[(M+N+1)::rows(data),.]
	Xref=dataref[1::M,.]
    Yref=dataref[(M+1)::(M+N),.]	
	Bref=dataref[(M+N+1)::rows(data),.]
	gX=gmat[1::M,.]
    gY=gmat[(M+1)::(M+N),.]	
	gB=gmat[(M+N+1)::rows(data),.]	
	
	class LinearProgram scalar q
	q = LinearProgram()
	
	c=(1,J(1,cols(dataref),0))
		lowerbd=.,J(1,length(c)-1,0)
		upperbd=J(1,length(c),.)		
		q.setCoefficients(c)
		q.setBounds(lowerbd, upperbd)
		if(maxiter!=-1){
		  q.setMaxiter(maxiter)
		}
		if (tol!=-1){
		  q.setTol(tol)
		}
	beta=J(cols(data),1,.)
	
	for(i=1;i<=cols(data);i++){
		Aie1=-gX[.,i],Xref
		bie1=X[.,i]
		Aie2=gY[.,i],-Yref
		bie2=-Y[.,i]
		Aec=-gB[.,i],Bref
		bec=B[.,i]
		if(rstype==0){
			q.setEquality(Aec,bec)		 	
		}
		else{
			lsum=0,J(1,N,1)
			q.setEquality(Aec \ lsum, bec \ 1)		 	
		}
		q.setInequality((Aie1 \ Aie2 ), (bie1 \ bie2))
        beta[i,1]=q.optimize()	    
	}
	
	return(beta)
	
	
}


real matrix function calnddf(real matrix data, 
                             real matrix dataref, 
                             real scalar  M,
							 real scalar  N,
							 real rowvector wmat,
							 real matrix gmat,
					         real scalar rstype, 					
					         real scalar maxiter, 
					         real scalar  tol    )
{
    ///nb=rows(data)-M-N
	X=data[1::M,.]
    Y=data[(M+1)::(M+N),.]	
	B=data[(M+N+1)::rows(data),.]
	Xref=dataref[1::M,.]
    Yref=dataref[(M+1)::(M+N),.]	
	Bref=dataref[(M+N+1)::rows(data),.]
	gX=gmat[1::M,.]
    gY=gmat[(M+1)::(M+N),.]	
	gB=gmat[(M+N+1)::rows(data),.]	
	
	class LinearProgram scalar q
	q = LinearProgram()	

	beta=J(cols(data),rows(data)+1,.)

    c=(wmat,J(1,cols(dataref),0))
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

          
    for(i=1;i<=cols(data);i++){

		Aie1=diag(-gX[.,i]),J(M,rows(data)-M,0),Xref
		bie1=X[.,i]
		Aie2=J(N,M,0),diag(gY[.,i]),J(N,rows(data)-M-N,0),-Yref
		bie2=-Y[.,i]
		Aec=J(rows(data)-M-N,M+N,0),diag(-gB[.,i]),Bref
		bec=B[.,i]
		

		q.setInequality((Aie1 \ Aie2 ), (bie1 \ bie2))

		if(rstype==0){
		q.setEquality(Aec,bec)		 	
		}
		else{
			lsum=J(1,length(wmat),0),J(1,cols(dataref),1)
			q.setEquality(Aec \ lsum, bec \ 1)		 	
		}
		
		b0=q.optimize()	

		if (q.converged()==1){
			bslk=q.parameters()
			beta[i,]=b0,bslk[1..rows(data)]
		}
		else{
			beta[i,]=J(1,1+rows(data),.)
		}				
        

	}

     return(beta)

	
	
}




									
real matrix function nddfcontemporaneous(   real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real rowvector    wmat, ///
											real matrix       gmat, ///
											real scalar       band, ///	
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol)									
	{


		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
		  
          
	      t0uniq=uniqrows(t0) 
		  D0=J(rows(data0),Mc-1,.)  
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],t1:==t0uniq[j])
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,.]= calnddf(XYB',XYBref',k1,k2,wmat,gmat2',rstype,maxiter,tol)
		  }
		  
		  return(D0)


}



       void function nddfcontemporaneousbt( real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real rowvector    wmat, ///
											real matrix       gmat, ///
											real scalar       band, ///	
										    real scalar       gamma, ///
										    real scalar       level, ///											
											real scalar       brep,
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol,  ///
											real scalar       nodots, ///
											string scalar     Dv)									
	{
         // sort by id t
		  struct bootdea scalar cires
          msub=ceil(length(uniqrows(dataref[.,1]))^gamma)
		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
	      t0uniq=uniqrows(t0) 
		  //Kr=length(uniqrows(dataref[.,1]))
		  //D0=J(rows(data0),1,.)
         st_view(D0=.,.,Dv)
         //g1[.,.]=D0
         st_view(Dbc=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         st_view(LL=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         st_view(UU=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         st_view(vari=.,.,"bootvar")
         //g5[.,.]=vari	            

		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],t1:==t0uniq[j])
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,.]= calnddf(XYB',XYBref',k1,k2,wmat,gmat2',rstype,maxiter,tol)
		  }
		  
		  id = dataref[.,1]
		  mm_panels(id,idinfo=.)
		  
		  Dboot=J(rows(data0),brep,.)
		  
		  for(i=1;i<=brep;i++){
			  if(nodots!=1){
				  displaydots(i)
			  }
		      idx=mm_sample(msub,.,idinfo,.,1)  // require sorting data by id time
			  datastar=dataref[idx,.]
			  Dbooti=nddfcontemporaneous(data0,datastar,k1,k2,wmat,gmat,band,rstype,maxiter,tol)
			  Dboot[.,i]=Dbooti[.,1]
		  }
		  
		  
		 for(j=1;j<=length(t0uniq);j++){
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  Kr=length(select(dataref[.,1],t1:==t0uniq[j]))
			  cires=biasandciddf(D0[pos,1],Dboot[pos,.], k1, cols(data0)-k1-2, length(pos), Kr, level,0,msub)
			  Dbc[pos,1]=cires.tebc
			  LL[pos,1] =cires.LL
			  UU[pos,1] =cires.UU
			  //bias[pos,1]=cires.bias
			  vari[pos,1]=cires.vari
		  }

         //st_view(g1=.,.,"Dv")
         //g1[.,.]=D0
         //st_view(g2=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         //st_view(g3=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         //st_view(g4=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         //st_view(g5=.,.,"bootvar")
         //g5[.,.]=vari	                 



}


									
									
real matrix function nddfbiennial(          real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real rowvector    wmat, ///
											real matrix       gmat, ///
											real scalar       band, ///	
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol)									
	{


		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
		  
          
	      t0uniq=uniqrows(t0) 
		  D0=J(rows(data0),Mc-1,.)  
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],(t1:>=t0uniq[j]):*(t1:<=t0uniq[j]+1))
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,.]= calnddf(XYB',XYBref',k1,k2,wmat,gmat2',rstype,maxiter,tol)
		  }
		  
		  return(D0)


}



real matrix function nddfbiennialbt(        real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real rowvector    wmat, ///
											real matrix       gmat, ///
											real scalar       band, ///	
										    real scalar       gamma, ///
										    real scalar       level, ///											
											real scalar       brep,
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol,  ///
											real scalar       nodots, ///
											string scalar     Dv)									
	{
         // sort by id t
		 
		  struct bootdea scalar cires
          msub=ceil(length(uniqrows(dataref[.,1]))^gamma)
		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
	      t0uniq=uniqrows(t0) 
		  //Kr=length(uniqrows(dataref[.,1]))
		  //D0=J(rows(data0),1,.)
         st_view(D0=.,.,Dv)
         //g1[.,.]=D0
         st_view(Dbc=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         st_view(LL=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         st_view(UU=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         st_view(vari=.,.,"bootvar")
         //g5[.,.]=vari	    
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],(t1:>=t0uniq[j]):* (t1:<=t0uniq[j]+1))
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,.]= calnddf(XYB',XYBref',k1,k2,wmat,gmat2',rstype,maxiter,tol)
		  }
		  
		  id = dataref[.,1]
		  mm_panels(id,idinfo=.)
		  
		  Dboot=J(rows(data0),brep,.)
		  
		  for(i=1;i<=brep;i++){
			  if(nodots!=1){
				  displaydots(i)
			  }			  
		      idx=mm_sample(msub,.,idinfo,.,1)
			  datastar=dataref[idx,.]
			  Dbooti=nddfbiennial(data0,datastar,k1,k2,wmat,gmat,band,rstype,maxiter,tol)
			  Dboot[.,i]=Dbooti[.,1]
		  }
		  
		  
		 for(j=1;j<=length(t0uniq);j++){
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  Kr=length(select(dataref[.,1],(t1:>=t0uniq[j]):* (t1:<=t0uniq[j]+1)))
			  cires=biasandciddf(D0[pos,1],Dboot[pos,.], k1, cols(data0)-k1-2, length(pos), Kr, level,0,msub)
			  Dbc[pos,1]=cires.tebc
			  LL[pos,1] =cires.LL
			  UU[pos,1] =cires.UU
			  //bias[pos,1]=cires.bias
			  vari[pos,1]=cires.vari
		  }
		  


}






real matrix function nddfsequential(        real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real rowvector    wmat, ///
											real matrix       gmat, ///
											real scalar       band, ///	
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol)									
	{


		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
          
	      t0uniq=uniqrows(t0) 
		  D0=J(rows(data0),Mc-1,.)  
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],t1:<=t0uniq[j])
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,.]= calnddf(XYB',XYBref',k1,k2,wmat,gmat2',rstype,maxiter,tol)
		  }
		  
		  return(D0)


}



real matrix function nddfsequentialbt(      real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real rowvector    wmat, ///
											real matrix       gmat, ///
											real scalar       band, ///	
										    real scalar       gamma, ///
										    real scalar       level, ///											
											real scalar       brep,  ///
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol, ///
											real scalar       nodots, ///
											string scalar     Dv)									
	{
         // sort by id t
		 
		  struct bootdea scalar cires
          msub=ceil(length(uniqrows(dataref[.,1]))^gamma)
		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
	      t0uniq=uniqrows(t0) 
		  //Kr=length(uniqrows(dataref[.,1]))
		  //D0=J(rows(data0),1,.)
         st_view(D0=.,.,Dv)
         //g1[.,.]=D0
         st_view(Dbc=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         st_view(LL=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         st_view(UU=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         st_view(vari=.,.,"bootvar")
         //g5[.,.]=vari	    
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],t1:<=t0uniq[j])
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,.]= calnddf(XYB',XYBref',k1,k2,wmat,gmat2',rstype,maxiter,tol)
		  }
		  
		  id = dataref[.,1]
		  mm_panels(id,idinfo=.)
		  
		  Dboot=J(rows(data0),brep,.)
		  
		  for(i=1;i<=brep;i++){
			  if(nodots!=1){
				  displaydots(i)
			  }			  
		      idx=mm_sample(msub,.,idinfo,.,1)
			  datastar=dataref[idx,.]
			  mata:Dbooti=nddfsequential(data0,datastar,k1,k2,wmat,gmat,band,rstype,maxiter,tol)
			  Dboot[.,i]=Dbooti[.,1]
		  }
		  
		  
		 for(j=1;j<=length(t0uniq);j++){
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  Kr=length(select(dataref[.,1],t1:<=t0uniq[j]))
			  cires=biasandciddf(D0[pos,1],Dboot[pos,.], k1, cols(data0)-k1-2, length(pos), Kr, level,0,msub*sum(t1:<=t0uniq[j]))
			  Dbc[pos,1]=cires.tebc
			  LL[pos,1] =cires.LL
			  UU[pos,1] =cires.UU
			  //bias[pos,1]=cires.bias
			  vari[pos,1]=cires.vari
		  }
		  


}




real matrix function nddfglobal(            real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real rowvector    wmat, ///
											real matrix       gmat, ///
											real scalar       band, ///	
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol)									
{


		XYB    = data0[.,3..cols(data0)]
	    XYBref = dataref[.,3..cols(dataref)]
		D0 = calnddf(XYB',XYBref',k1,k2,wmat,gmat',rstype,maxiter,tol)
		//D0
		return(D0)
		


}



real matrix function nddfglobalbt(          real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real rowvector    wmat, ///
											real matrix       gmat, ///
											real scalar       band, ///	
										    real scalar       gamma, ///
										    real scalar       level, ///											
											real scalar       brep,
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol,   ///
											real scalar       nodots, ///
											string scalar     Dv)									
	{
         // sort by id t
		 
		  struct bootdea scalar cires
          msub=ceil(length(uniqrows(dataref[.,1]))^gamma)
		  t0  = data0[.,2]
		  //t1  = dataref[.,2]
          //Mc=cols(data0) 
	      t0uniq=uniqrows(t0) 
		  Kr=rows(dataref)
		
		  //D0=J(rows(data0),1,.)
         st_view(D0=.,.,Dv)
         //g1[.,.]=D0
         st_view(Dbc=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         st_view(LL=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         st_view(UU=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         st_view(vari=.,.,"bootvar")
         //g5[.,.]=vari	    		  
		XYB    = data0[.,3..cols(data0)]
		XYBref = dataref[.,3..cols(dataref)]
		D0[.,.] = calnddf(XYB',XYBref',k1,k2,wmat,gmat',rstype,maxiter,tol)
		  
		  
		  id = dataref[.,1]
		  mm_panels(id,idinfo=.)
		  
		  Dboot=J(rows(data0),brep,.)
		  msubmat=J(1,brep,.)
		  for(i=1;i<=brep;i++){
			  if(nodots!=1){
				  displaydots(i)
			  }			      
			  
		      idx=mm_sample(msub,.,idinfo,.,1)
			  datastar=dataref[idx,.]
			  //rows(datastar)
			  msubmat[1,i]=rows(datastar)
			  Dbooti=nddfglobal(data0,datastar,k1,k2,wmat,gmat,band,rstype,maxiter,tol)
			  Dboot[.,i]=Dbooti[.,1]
		  }
		  
		  
		 for(j=1;j<=length(t0uniq);j++){
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  cires=biasandciddf(D0[pos,1],Dboot[pos,.], k1, cols(data0)-k1-2, length(pos), Kr, level,0,msubmat[1,j])
			  Dbc[pos,1]=cires.tebc
			  LL[pos,1] =cires.LL
			  UU[pos,1] =cires.UU
			 // bias[pos,1]=cires.bias
			  vari[pos,1]=cires.vari
		  }
		  


}







real matrix function nddfwindow(            real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real rowvector    wmat, ///
											real matrix       gmat, ///
											real scalar       band, ///	
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol)									
	{


		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
		  
          
	      t0uniq=uniqrows(t0) 
		  D0=J(rows(data0),Mc-1,.)  
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],t1:<=(t0uniq[j]+band):*(t1:>=t0uniq[j]-band))
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,.]= calnddf(XYB',XYBref',k1,k2,wmat,gmat2',rstype,maxiter,tol)
		  }
		  
		  return(D0)


}



real matrix function nddfwindowbt(          real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real scalar       wmat, ///
											real matrix       gmat, ///
											real scalar       band, ///	
										    real scalar       gamma, ///
										    real scalar       level, ///											
											real scalar       brep,
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol, ///
											real scalar       nodots, ///
											string scalar     Dv)									
	{
         // sort by id t
		 
		  struct bootdea scalar cires
          msub=ceil(length(uniqrows(dataref[.,1]))^gamma)
		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
	      t0uniq=uniqrows(t0) 
		  //Kr=length(uniqrows(dataref[.,1]))
		  //D0=J(rows(data0),1,.)
         st_view(D0=.,.,Dv)
         //g1[.,.]=D0
         st_view(Dbc=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         st_view(LL=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         st_view(UU=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         st_view(vari=.,.,"bootvar")
         //g5[.,.]=vari	     
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],(t1:<=t0uniq[j]+band):*(t1:>=t0uniq[j]-band))
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,.]= calnddf(XYB',XYBref',k1,k2,wmat,gmat2',rstype,maxiter,tol)
		  }
		  
		  id = dataref[.,1]
		  mm_panels(id,idinfo=.)
		  
		  Dboot=J(rows(data0),brep,.)
		  
		  for(i=1;i<=brep;i++){
			  if(nodots!=1){
				  displaydots(i)
			  }			  
		      idx=mm_sample(msub,.,idinfo,.,1)
			  datastar=dataref[idx,.]
			  Dbooti=nddfwindow(data0,datastar,k1,k2,wmat,gmat,band,rstype,maxiter,tol)
			  Dboot[.,i]=Dbooti[.,1]
		  }
		  
		  
		 for(j=1;j<=length(t0uniq);j++){
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  Kr=length(select(dataref[.,1],(t1:<=t0uniq[j]+band):*(t1:>=t0uniq[j]-band)))
			  cires=biasandciddf(D0[pos,1],Dboot[pos,.], k1, cols(data0)-k1-2, length(pos), Kr, level,0,msub*sum((t1:<=t0uniq[j]+band):* (t1:>=t0uniq[j]-band)))
			  Dbc[pos,1]=cires.tebc
			  LL[pos,1] =cires.LL
			  UU[pos,1] =cires.UU
			  //bias[pos,1]=cires.bias
			  vari[pos,1]=cires.vari
		  }
		  


}

//////////////////////////////////////


void function _ddf( string scalar     d, ///
                    real scalar       k1, ///
                    real scalar       k2, ///
					string scalar     touse1,  ///
					string scalar     touse2, ///
					string scalar     gname, ///
					string scalar     bname, ///
					real scalar       rstype, ///					
					real scalar       maxiter, ///
					real scalar       tol)
	{


          data=st_data(.,d,touse1)
          data=data'
          M=rows(data)          
          g=st_data(.,gname,touse1)
          g=g'
        
          dataref=st_data(.,d,touse2)
          dataref=dataref'

          Xref=dataref[1..k1,.]
          Yref=dataref[k1+1..(k1+k2),.]
          Bref=dataref[(k1+k2+1)..M,.]
          X=data[1..k1,.]
          Y=data[k1+1..(k1+k2),.]
          B=data[(k1+k2+1)..M,.]
          theta=J(cols(data),1,.)         
		 
		 for(j=1;j<=cols(X);j++){
			theta[j]=_ddf2(X[.,j],Y[.,j],B[.,j],Xref,Yref,Bref,g[.,j],rstype,maxiter,tol)
		 }

          st_view(gen=.,.,bname,touse1)
          gen[.,.]=theta
		

}




									
									
real matrix function ddfcontemporaneous(    real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real matrix       gmat, ///
											real scalar       band, ///	
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol)									
	{


		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
		  
          
	      t0uniq=uniqrows(t0) 
		  D0=J(rows(data0),1,.)  
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],t1:==t0uniq[j])
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,1]= calddf(XYB',XYBref',k1,k2,gmat2',rstype,maxiter,tol)
		  }
		  
		  return(D0)


}



       void function ddfcontemporaneousbt(  real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real matrix       gmat, ///
											real scalar       band, ///	
										    real scalar       gamma, ///
										    real scalar       level, ///											
											real scalar       brep,
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol,  ///
											real scalar       nodots)									
	{
         // sort by id t
		  struct bootdea scalar cires
          msub=ceil(length(uniqrows(dataref[.,1]))^gamma)
		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
	      t0uniq=uniqrows(t0) 
		  //Kr=length(uniqrows(dataref[.,1]))
		  //D0=J(rows(data0),1,.)
         st_view(D0=.,.,"Dv")
         //g1[.,.]=D0
         st_view(Dbc=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         st_view(LL=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         st_view(UU=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         st_view(vari=.,.,"bootvar")
         //g5[.,.]=vari	            

		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],t1:==t0uniq[j])
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,1]= calddf(XYB',XYBref',k1,k2,gmat2',rstype,maxiter,tol)
		  }
		  
		  id = dataref[.,1]
		  mm_panels(id,idinfo=.)
		  
		  Dboot=J(rows(data0),brep,.)
		  
		  for(i=1;i<=brep;i++){
			  if(nodots!=1){
				  displaydots(i)
			  }
		      idx=mm_sample(msub,.,idinfo,.,1)
			  datastar=dataref[idx,.]
			  Dboot[.,i]=ddfcontemporaneous(data0,datastar,k1,k2,gmat,band,rstype,maxiter,tol)
		  }
		  
		  
		 for(j=1;j<=length(t0uniq);j++){
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  Kr=length(select(dataref[.,1],t1:==t0uniq[j]))
			  cires=biasandciddf(D0[pos,1],Dboot[pos,.], k1, cols(data0)-k1-2, length(pos), Kr, level,0,msub)
			  Dbc[pos,1]=cires.tebc
			  LL[pos,1] =cires.LL
			  UU[pos,1] =cires.UU
			  //bias[pos,1]=cires.bias
			  vari[pos,1]=cires.vari
		  }

         //st_view(g1=.,.,"Dv")
         //g1[.,.]=D0
         //st_view(g2=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         //st_view(g3=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         //st_view(g4=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         //st_view(g5=.,.,"bootvar")
         //g5[.,.]=vari	                 



}


									
									
real matrix function ddfbiennial(           real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real matrix       gmat, ///
											real scalar       band, ///	
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol)									
	{


		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
		  
          
	      t0uniq=uniqrows(t0) 
		  D0=J(rows(data0),1,.)  
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],(t1:>=t0uniq[j]):*(t1:<=t0uniq[j]+1))
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,1]= calddf(XYB',XYBref',k1,k2,gmat2',rstype,maxiter,tol)
		  }
		  
		  return(D0)


}



real matrix function ddfbiennialbt(         real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real matrix       gmat, ///
											real scalar       band, ///	
										    real scalar       gamma, ///
										    real scalar       level, ///											
											real scalar       brep,
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol, ///
											real scalar       nodots)									
	{
         // sort by id t
		 
		  struct bootdea scalar cires
          msub=ceil(length(uniqrows(dataref[.,1]))^gamma)
		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
	      t0uniq=uniqrows(t0) 
		  //Kr=length(uniqrows(dataref[.,1]))
		  //D0=J(rows(data0),1,.)
         st_view(D0=.,.,"Dv")
         //g1[.,.]=D0
         st_view(Dbc=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         st_view(LL=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         st_view(UU=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         st_view(vari=.,.,"bootvar")
         //g5[.,.]=vari	    
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],(t1:>=t0uniq[j]):* (t1:<=t0uniq[j]+1))
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,1]= calddf(XYB',XYBref',k1,k2,gmat2',rstype,maxiter,tol)
		  }
		  
		  id = dataref[.,1]
		  mm_panels(id,idinfo=.)
		  
		  Dboot=J(rows(data0),brep,.)
		  
		  for(i=1;i<=brep;i++){
			  if(nodots!=1){
				  displaydots(i)
			  }			  
		      idx=mm_sample(msub,.,idinfo,.,1)
			  datastar=dataref[idx,.]
			  Dboot[.,i]=ddfbiennial(data0,datastar,k1,k2,gmat,band,rstype,maxiter,tol)
		  }
		  
		  
		 for(j=1;j<=length(t0uniq);j++){
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  Kr=length(select(dataref[.,1],(t1:>=t0uniq[j]):* (t1:<=t0uniq[j]+1)))
			  cires=biasandciddf(D0[pos,1],Dboot[pos,.], k1, cols(data0)-k1-2, length(pos), Kr, level,0,msub)
			  Dbc[pos,1]=cires.tebc
			  LL[pos,1] =cires.LL
			  UU[pos,1] =cires.UU
			  //bias[pos,1]=cires.bias
			  vari[pos,1]=cires.vari
		  }
		  


}






real matrix function ddfsequential(         real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real matrix       gmat, ///
											real scalar       band, ///	
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol)									
	{


		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
          
	      t0uniq=uniqrows(t0) 
		  D0=J(rows(data0),1,.)  
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],t1:<=t0uniq[j])
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,1]= calddf(XYB',XYBref',k1,k2,gmat2',rstype,maxiter,tol)
		  }
		  
		  return(D0)


}



real matrix function ddfsequentialbt(       real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real matrix       gmat, ///
											real scalar       band, ///	
										    real scalar       gamma, ///
										    real scalar       level, ///											
											real scalar       brep,
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol, ///
											real scalar       nodots)									
	{
         // sort by id t
		 
		  struct bootdea scalar cires
          msub=ceil(length(uniqrows(dataref[.,1]))^gamma)
		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
	      t0uniq=uniqrows(t0) 
		  //Kr=length(uniqrows(dataref[.,1]))
		  //D0=J(rows(data0),1,.)
         st_view(D0=.,.,"Dv")
         //g1[.,.]=D0
         st_view(Dbc=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         st_view(LL=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         st_view(UU=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         st_view(vari=.,.,"bootvar")
         //g5[.,.]=vari	    
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],t1:<=t0uniq[j])
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,1]= calddf(XYB',XYBref',k1,k2,gmat2',rstype,maxiter,tol)
		  }
		  
		  id = dataref[.,1]
		  mm_panels(id,idinfo=.)
		  
		  Dboot=J(rows(data0),brep,.)
		  
		  for(i=1;i<=brep;i++){

			  if(nodots!=1){
				  displaydots(i)
			  }			  
		      idx=mm_sample(msub,.,idinfo,.,1)
			  datastar=dataref[idx,.]
			  Dboot[.,i]=ddfsequential(data0,datastar,k1,k2,gmat,band,rstype,maxiter,tol)
		  }
		  
		  
		 for(j=1;j<=length(t0uniq);j++){
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  Kr=length(select(dataref[.,1],t1:<=t0uniq[j]))
			  cires=biasandciddf(D0[pos,1],Dboot[pos,.], k1, cols(data0)-k1-2, length(pos), Kr, level,0,msub*sum(t1:<=t0uniq[j]))
			  Dbc[pos,1]=cires.tebc
			  LL[pos,1] =cires.LL
			  UU[pos,1] =cires.UU
			  //bias[pos,1]=cires.bias
			  vari[pos,1]=cires.vari
		  }
		  


}




real matrix function ddfglobal(             real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real matrix       gmat, ///
											real scalar       band, ///	
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol)									
{


		XYB    = data0[.,3..cols(data0)]
	    XYBref = dataref[.,3..cols(dataref)]
		D0 = calddf(XYB',XYBref',k1,k2,gmat',rstype,maxiter,tol)
		//D0
		return(D0)
		


}



real matrix function ddfglobalbt(           real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real matrix       gmat, ///
											real scalar       band, ///	
										    real scalar       gamma, ///
										    real scalar       level, ///											
											real scalar       brep,  ///
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol, ///
											real scalar       nodots)									
	{
         // sort by id t
		 
		  struct bootdea scalar cires
          msub=ceil(length(uniqrows(dataref[.,1]))^gamma)
		  t0  = data0[.,2]
		  //t1  = dataref[.,2]
          //Mc=cols(data0) 
	      t0uniq=uniqrows(t0) 
		  Kr=rows(dataref)
		
		  //D0=J(rows(data0),1,.)
         st_view(D0=.,.,"Dv")
         //g1[.,.]=D0
         st_view(Dbc=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         st_view(LL=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         st_view(UU=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         st_view(vari=.,.,"bootvar")
         //g5[.,.]=vari	    		  
		XYB    = data0[.,3..cols(data0)]
		XYBref = dataref[.,3..cols(dataref)]
		D0[.,1] = calddf(XYB',XYBref',k1,k2,gmat',rstype,maxiter,tol)
		  
		  
		  id = dataref[.,1]
		  mm_panels(id,idinfo=.)
		  
		  Dboot=J(rows(data0),brep,.)
		  msubmat=J(1,brep,.)
		  for(i=1;i<=brep;i++){
			  if(nodots!=1){
				  displaydots(i)
			  }			  
		      idx=mm_sample(msub,.,idinfo,.,1)
			  datastar=dataref[idx,.]
			  //rows(datastar)
			  msubmat[1,i]=rows(datastar)
			  Dboot[.,i]=ddfglobal(data0,datastar,k1,k2,gmat,band,rstype,maxiter,tol)
		  }
		  
		  
		 for(j=1;j<=length(t0uniq);j++){
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  cires=biasandciddf(D0[pos,1],Dboot[pos,.], k1, cols(data0)-k1-2, length(pos), Kr, level,0,msubmat[1,j])
			  Dbc[pos,1]=cires.tebc
			  LL[pos,1] =cires.LL
			  UU[pos,1] =cires.UU
			 // bias[pos,1]=cires.bias
			  vari[pos,1]=cires.vari
		  }
		  


}







real matrix function ddfwindow(             real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real matrix       gmat, ///
											real scalar       band, ///	
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol)									
	{


		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
		  
          
	      t0uniq=uniqrows(t0) 
		  D0=J(rows(data0),1,.)  
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],t1:<=(t0uniq[j]+band):*(t1:>=t0uniq[j]-band))
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,1]= calddf(XYB',XYBref',k1,k2,gmat2',rstype,maxiter,tol)
		  }
		  
		  return(D0)


}



real matrix function ddfwindowbt(           real matrix       data0, ///
											real matrix       dataref, ///
											real scalar       k1, ///
											real scalar       k2, ///
											real matrix       gmat, ///
											real scalar       band, ///	
										    real scalar       gamma, ///
										    real scalar       level, ///											
											real scalar       brep,
											real scalar       rstype, ///					
											real scalar       maxiter, ///
											real scalar       tol, ///
											real scalar       nodots)									
	{
         // sort by id t
		 
		  struct bootdea scalar cires
          msub=ceil(length(uniqrows(dataref[.,1]))^gamma)
		  t0  = data0[.,2]
		  t1  = dataref[.,2]
          Mc=cols(data0) 
	      t0uniq=uniqrows(t0) 
		  //Kr=length(uniqrows(dataref[.,1]))
		  //D0=J(rows(data0),1,.)
         st_view(D0=.,.,"Dv")
         //g1[.,.]=D0
         st_view(Dbc=.,.,"Dv_bc")
         //g2[.,.]=Dbc
         st_view(LL=.,.,"Dv_lower")
         //g3[.,.]=LL		 
         st_view(UU=.,.,"Dv_upper")
         //g4[.,.]=UU	  
         st_view(vari=.,.,"bootvar")
         //g5[.,.]=vari	     
		  for(j=1;j<=length(t0uniq);j++){
		      XYB    = select(data0[.,3..Mc],t0:==t0uniq[j])
			  XYBref = select(dataref[.,3..Mc],(t1:<=t0uniq[j]+band) + (t1:>=t0uniq[j]-band))
			  gmat2  = select(gmat,t0:==t0uniq[j])
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  D0[pos,1]= calddf(XYB',XYBref',k1,k2,gmat2',rstype,maxiter,tol)
		  }
		  
		  id = dataref[.,1]
		  mm_panels(id,idinfo=.)
		  
		  Dboot=J(rows(data0),brep,.)
		  
		  for(i=1;i<=brep;i++){
			  if(nodots!=1){
				  displaydots(i)
			  }			  
		      idx=mm_sample(msub,.,idinfo,.,1)
			  datastar=dataref[idx,.]
			  Dboot[.,i]=ddfwindow(data0,datastar,k1,k2,gmat,band,rstype,maxiter,tol)
		  }
		  
		  
		 for(j=1;j<=length(t0uniq);j++){
			  pos=select(1::rows(data0),t0:==t0uniq[j])
			  Kr=length(select(dataref[.,1],t1:<=t0uniq[j]))
			  cires=biasandciddf(D0[pos,1],Dboot[pos,.], k1, cols(data0)-k1-2, length(pos), Kr, level,0,msub*sum((t1:<=t0uniq[j]+band):* (t1:>=t0uniq[j]-band)))
			  Dbc[pos,1]=cires.tebc
			  LL[pos,1] =cires.LL
			  UU[pos,1] =cires.UU
			  //bias[pos,1]=cires.bias
			  vari[pos,1]=cires.vari
		  }
		  


}





function biasandciddf(  real colvector te,real matrix teboot,  ///
                        real scalar M, real scalar N, real scalar K, real scalar Kr, ///
                        real scalar level,real scalar smoothed, real scalar msub)
{
    if(smoothed){
	    con1=1
		con2=1
	}
	else{
	    con1=msub^(2/(M+N+1))
		con2=Kr^(-2/(M+N+1))
	}
	con3=con1*con2
    bias=J(K,1,.)
	vari=J(K,1,.)
	BovV=J(K,1,.)
	tebc=J(K,1,.)
	UU=J(K,1,.)
	LL=J(K,1,.)
	
	for(i=1;i<=K;i++){
	    TEi=teboot[i,.]
		TEi=select(TEi,TEi:!=.)
		if(length(TEi)<100){
		    bias[i]=.
			vari[i]=.
			BovV[i]=.
			tebc[i]=.
			LL[i]=.
			UU[i]=.
		}
		else{
			
		    bias[i]=con3*(mean(TEi')-te[i])
			vari[i]=variance(TEi')
			BovV[i]=(bias[i])^2 / vari[i] * 3
            tebc[i] = te[i] - bias[i]
            teOverTEhat = (TEi':- te[i] ) * con1	
		    quans =mm_quantile(teOverTEhat,1, ((0.5-level/200)\(0.5+level/200) ))
			//quans* con2
		    LL[i] = te[i] - quans[2] * con2
		    UU[i] = te[i] - quans[1] * con2			
			//quans =mm_quantile(TEi',1, ((0.5-level/200)\(0.5+level/200) ))
			//LL[i] =  quans[1] 
		    //UU[i] =  quans[2] 	
		}
	}
		struct bootdea scalar res
				   res.tebc      = tebc
				   res.vari      = vari
				   res.BovV    = BovV			   
				   res.bias    = bias
				   res.LL 	   = LL
				   res.UU	   = UU
 
        return(res)
}









real scalar function _ddf2( real colvector    X, ///
                            real colvector    Y, ///
						    real colvector    B, ///
						    real matrix       Xref, ///
						    real matrix       Yref, ///
						    real matrix       Bref, ///
						    real colvector    g, ///
					        real scalar       rstype, ///					
					        real scalar       maxiter, ///
					        real scalar       tol)
	{


          k1=length(X)
		  k2=length(Y)
		  M=length(B)+k1+k2
          gX=g[1..k1,1]
          gY=g[k1+1..(k1+k2),1]
          gB=g[(k1+k2+1)..M,1]         
          N=cols(Xref)

	    class LinearProgram scalar q
	    q = LinearProgram()
		c=(1,J(1,N,0))
		lowerbd=.,J(1,length(c)-1,0)
		upperbd=J(1,length(c),.)		
		q.setCoefficients(c)
		q.setBounds(lowerbd, upperbd)

		if(maxiter!=-1){
		  q.setMaxiter(maxiter)
		}
		if (tol!=-1){
		  q.setTol(tol)
		}			

         
  
		Aie1=-gX,Xref
		bie1=X
		Aie2=gY,-Yref
		bie2=-Y
		Aec=-gB,Bref
		bec=B
		if(rstype==0){
			q.setEquality(Aec,bec)		 	
		}
		else{
			lsum=0,J(1,N,1)
			q.setEquality(Aec \ lsum, bec \ 1)		 	
		}
		q.setInequality((Aie1 \ Aie2 ), (bie1 \ bie2))
        beta=q.optimize()		  
         
        return(beta)
		

}




void function _nddf( string scalar     d, ///
                    real scalar       nx, ///
                    real scalar       ny, ///
					string scalar     touse1,  ///
					string scalar     touse2, ///
					string scalar     gname, ///
					string scalar     wname, ///
					string scalar     bname, ///
					real scalar       rstype, ///					
					real scalar       maxiter, ///
					real scalar       tol)
	{

          wmat=st_matrix(wname)
          data=st_data(.,d,touse1)
          data=data'
		  M=rows(data)
          g=st_data(.,gname,touse1)
          g=g'       
          dataref=st_data(.,d,touse2)
          dataref=dataref'
          Xref=dataref[1..nx,.]
          Yref=dataref[nx+1..(nx+ny),.]
          Bref=dataref[(nx+ny+1)..M,.]
          X=data[1..nx,.]
          Y=data[nx+1..(nx+ny),.]
          B=data[(nx+ny+1)..M,.]
 
          theta=J(cols(data),1,.)
  
          for(j=1;j<=cols(data);j++){

                 theta[j]=_nddf2(X[.,j],Y[.,j],B[.,j],Xref,Yref,Bref,g[.,j],wmat,rstype,maxiter,tol)	  
               
         }

          st_view(gen=.,.,bname,touse1)
          gen[.,.]=theta
		

}



real function _nddf2( real colvector     X, ///
                      real colvector     Y, ///
                      real colvector     B, ///
					  real matrix        Xref,  ///
					  real matrix        Yref,  ///
					  real matrix        Bref,  ///
					  real colvector      g, ////
					  real rowvector     wmat, ///
					  real scalar       rstype, ///					
					  real scalar       maxiter, ///
					  real scalar       tol)
	{

          
          nx=length(X)
		  ny=length(Y)
		  nb=length(B)
		  M=nb+nx+ny
          gX=g[1..nx,1]
          gY=g[nx+1..(nx+ny),1]
          gB=g[(nx+ny+1)..M,1]        
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
				 lsum=J(1,length(wmat),0),J(1,N,1)

				if(rstype==0){
				  q.setEquality(Aec,bec)		 	
				}
				else{
					 lsum=J(1,M,0),J(1,N,1)
					 q.setEquality(Aec \ lsum, bec \ 1)		 	
				 }

				 q.setInequality((Aie1 \ Aie2 ), (bie1 \ bie2))

                 beta=q.optimize()		  
                 //q.getInequality()
                 return(beta)
		

}


void function _tenddf( string scalar     d, ///
                    real scalar       nx, ///
                    real scalar       ny, ///
					string scalar     touse1,  ///
					string scalar     touse2, ///
					string scalar     gname, ///
					string scalar     wname, ///
					string scalar     bname, ///
					real scalar       rstype, ///					
					real scalar       maxiter, ///
					real scalar       tol)
	{

          wmat=st_matrix(wname)
          data=st_data(.,d,touse1)
          data=data'
		  M=rows(data)
          g=st_data(.,gname,touse1)
          g=g'       
          dataref=st_data(.,d,touse2)
          dataref=dataref'
          Xref=dataref[1..nx,.]
          Yref=dataref[nx+1..(nx+ny),.]
          Bref=dataref[(nx+ny+1)..M,.]
          X=data[1..nx,.]
          Y=data[nx+1..(nx+ny),.]
          B=data[(nx+ny+1)..M,.]
 
          theta=J(cols(data),rows(data)+1,.)
  
          for(j=1;j<=cols(data);j++){

                 theta[j,.]=_tenddf2(X[.,j],Y[.,j],B[.,j],Xref,Yref,Bref,g[.,j],wmat,rstype,maxiter,tol)	  
                 //q.getInequality()
         }

          st_view(gen=.,.,bname,touse1)
          gen[.,.]=theta
		

}



real function _tenddf2( real colvector     X, ///
                      real colvector     Y, ///
                      real colvector     B, ///
					  real matrix        Xref,  ///
					  real matrix        Yref,  ///
					  real matrix        Bref,  ///
					  real colvector      g, ////
					  real rowvector     wmat, ///
					  real scalar       rstype, ///					
					  real scalar       maxiter, ///
					  real scalar       tol)
	{

          
          nx=length(X)
		  ny=length(Y)
		  nb=length(B)
		  M=nb+nx+ny
          gX=g[1..nx,1]
          gY=g[nx+1..(nx+ny),1]
          gB=g[(nx+ny+1)..M,1]        
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
	lsum=J(1,length(wmat),0),J(1,N,1)

	if(rstype==0){
		 q.setEquality(Aec,bec)		 	
	}
	else{
		 lsum=J(1,M,0),J(1,N,1)
		 q.setEquality(Aec \ lsum, bec \ 1)		 	
	 }

	 q.setInequality((Aie1 \ Aie2 ), (bie1 \ bie2))

	 dval=q.optimize()

	 if (q.converged()==1){
	 	bslk=q.parameters()
	 	beta=dval,bslk[1..M]
	 }
	 else{
	 	beta=J(1,1+M,.)
	 }
              	  
     //q.getInequality()
     return(beta)
		

}



end