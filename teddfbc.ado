*! version 2.1, 7 Aug 2020
*! version 2.0
* By Kerry Du & Daoping Wang, 30 Jul 2020 
**
* 
capture program drop teddfbc
program define teddfbc, rclass
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
	
	
    syntax varlist [if] [in], [rf(varname)  gx(varlist) gy(varlist) gb(varlist)  ///
	                                         VRS  NONRadial  ///
										   Wmat(string) SAVing(string)  BREP(integer 499) Level(real 95)   ///
										   maxiter(numlist integer >0 max=1) tol(numlist max=1)]
 
									   
	marksample touse 
	markout `touse' `invars' `gopvars' `gx' `gy' `gb'
	local bopvars `varlist'
	local invars: list uniq invars
	local gopvars: list uniq gopvars
	local bopvars: list uniq bopvars
	
	local ninp: word count `invars'
    local ngo: word count `gopvars'
    local nbo: word count `bopvars'
	local nvar=`ninp'+`ngo'+`nbo'
	
	confirm numeric var `invars' `gopvars' `bopvars' `rf'

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
			disp as red "wmat() should be only used when nonradial is specified."
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
            //mat list `weightvec',noblank  noheader		
	}
		
	
	if "`maxiter'"==""{
		local maxiter=-1
	}
	if "`tol'"==""{
		local tol=-1
	}	
	  
	
	preserve
	


	if "`gx'"!=""{
		local ngx: word count `gx'
		if `ngx'!=`ninp'{
		    disp as error "# of input variables != # of variables specified in gx()."
		    error 498
		}
		*gx should not be large than 0
		foreach v of local invars{
			cap assert `v'<=0
			if _rc {
				dis as red "Error: Some of gx() > 0."
				exit
			}
		}
		*		
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
		foreach v of local gopvars{
			cap assert `v'>=0
			if _rc {
				dis as red "Error: some of gy() < 0."
				exit
			}

		local gmat `gmat' `gy'
		local gmatname `gmatname' `gy'
	  }
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
		foreach v of local bopvars{
			cap assert `v'<=0
			if _rc {
				dis as red "Error: some of gb() > 0."
				exit
			}
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
	

	qui keep   `invars' `gopvars' `bopvars' `dmu' `time' `gmat' `touse' `rf'
	qui gen Row=_n
	
	label var Row "Row #"
	qui gen double Dval=.
	*label var Dval "Value of DDFs"
	
	tempvar touse2  
    qui gen byte `touse2'=1 
    if "`rf'"!=""{
    	qui replace `touse2'=(`rf'!=0) if !missing(`rf')
    }
    markout `touse2' `invars' `opvars' `badvars' `rf'
    *count if `touse2'		
	
	qui mata mata mlib index	
	local data `invars' `gopvars' `bopvars'

		
		mata: te=te_ddf("`data'",`ninp',`ngo',"`touse'", "`touse2'","`gmat'","Dval",`rts',`maxiter',`tol')
		mata: data0   = st_data(.,"`data'","`touse'")
		mata: dataref = st_data(.,"`data'","`touse2'")
		mata: gmat    = st_data(.,"`gmat'","`touse'")
		qui keep if `touse'
		mata: bootddf(te,data0,dataref,gmat,`ninp',`ngo',`level',`brep',`rts',`maxiter',`tol')
		
    
	//list, sep(0) 
	
	restore 
	
	end
	
	
///////////////////////////////////////////////////////	
cap mata mata drop te_ddf()
cap mata mata drop bootddf()
cap mata mata drop calddf()
cap mata mata drop smplHomTEddf()
cap mata mata drop biasandciddf()

mata:

real matrix function te_ddf(string scalar     d, ///
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
		  return(theta)
		
}



 void function bootddf(real matrix te,
                       real matrix data, 
					   real matrix dataref, 
                       real matrix gmat,
					   real scalar M, 
					   real scalar N, 
					   real scalar level,
					   real scalar brep,
					   real scalar rstype, 					   
					   real scalar maxiter, 
					   real scalar  tol    )
 {
     K=rows(data)
	 //Kr=rows(dataref)
	 teboot=J(K,brep,.)
	 teRef    = calddf(dataref',dataref',M,N,gmat',rstype,maxiter,tol)
	 //teRef
	 flag=select(1::rows(dataref),teRef:!=.)
	 teRef=teRef[flag]
	 dataref=dataref[flag,.]
	 Kr=rows(dataref)
	 terfl    = teRef \ 2:- teRef
     mybw     = kdens_bw( terfl, 1, "sjpi")
     scVarHom = ( 1 + mybw^2 / variance(terfl) )^(-1/2)
	 
	 for(j=1;j<=brep;j++){
		 tStar    = smplHomTEddf(terfl, Kr, mybw, scVarHom)
		 datastar = dataref+J(1,cols(data),teRef-tStar):*gmat	  
		 teboot[.,j]   = calddf(data',datastar',M,N,gmat',rstype,maxiter,tol)
		// teboot[.,j]
   }
	
   biasandciddf(te,teboot, M, cols(data)-M,  K,  Kr, level,1,0)	
   
}


real matrix function calddf(real matrix data, 
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
	//lowerbd=.,J(1,length(c)-1,0)
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
		//q.query()
        beta[i,1]=q.optimize()	    
	}
	
	return(beta)
	
	
}


real matrix function smplHomTEddf(real colvector terfl, 
                                  real scalar    Kr, 
                                  real scalar    mybw, 
                                  real scalar scVarHom )
{
 newSample = mm_sample(Kr,2*Kr)
 bStar     =  terfl[newSample]
 epsStar   = rnormal(Kr,1,0,1)
 tStar     = bStar + mybw * epsStar
 tStar     =  (tStar:>= 1):*(2:-tStar)+(tStar:<1):*tStar
 // correct the sequence
 tStar     = scVarHom*(tStar:-mean(bStar)):+mean(bStar)
 return(tStar)
}

 
void function biasandciddf(real colvector te,real matrix teboot,  ///
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
			//teOverTEhat
		    quans =mm_quantile(teOverTEhat,1, ((0.5-level/200)\(0.5+level/200) ))
			//quans* con2
		    //LL[i] = te[i] - quans[2] * con2
		    //UU[i] = te[i] - quans[1] * con2			
			quans =mm_quantile(TEi',1, ((0.5-level/200)\(0.5+level/200) ))
			LL[i] =  quans[1] 
		    UU[i] =  quans[2] 	
		}
	}
	te
	LL
	UU
	/*
	st_addvar("double","Dvalbc")
	st_store(.,"Dvalbc",tebc)
	st_addvar("double","biassqvarc")
	st_store(.,"biassqvar",BovV)	
	st_addvar("double","LL")
	st_store(.,"LL",LL)	
	st_addvar("double","UU")
	st_store(.,"UU",UU)	
	st_addvar("double","varboot")
	st_store(.,"varboot",vari)	
	*/
}






end







///////////////////////////////////////////////////////
/*

	else{

		label var Dval "Value of DDF"
		
		mata: _ddf("`data'",`ninp',`ngo',"`touse'", "`touse2'","`gmat'","Dval",`rts',`maxiter',`tol')
		mata: _ddf("`data'",`ninp',`ngo',"`touse2'", "`touse2'","`gmat'","Dvalref",`rts',`maxiter',`tol')
		
		if "`subsample'"==""{
			local jj=0
			foreach v in `invars'{
				*tempvar `v'star
				*local gg: word `jj' of `gmat'
				*cap assert `jj'==0 
				*if _rc!=0{
				*	qui gen double ``v'star'=`v'+Dvalref*`jj' if `touse2'
				*    local datastar `datastar' ``v'star'
				*}
				*else{
				*	local datastar `datastar' `v'
				*}
				tempvar `v'beta
				local j++
				local gg: word `jj' of `gmat'
				qui gen double ``v'beta'=(`v'+Dvalref*`jj')/`v' if `touse2'
			}
			foreach v in `outvars'{
				*tempvar `v'star
				*local gg: word `jj' of `gmat'
				*cap assert `jj'==0 
				*if _rc!=0{
				*	qui gen double ``v'star'=`v'+Dvalref*`jj' if `touse2'
				*    local datastar `datastar' ``v'star'
				*}
				*else{
				*	local datastar `datastar' `v'
				*}
				tempvar `v'beta
				local j++
				local gg: word `jj' of `gmat'
				qui gen double ``v'beta'=`v'/(`v'+Dvalref*`jj') if `touse2'
			}
			foreach v in `badvars'{
				*tempvar `v'star
				*local gg: word `jj' of `gmat'
				*cap assert `jj'==0 
				*if _rc!=0{
				*	qui gen double ``v'star'=`v'+Dvalref*`jj' if `touse2'
				*    local datastar `datastar' ``v'star'
				*}
				*else{
				*	local datastar `datastar' `v'
				*}
				tempvar `v'beta
				local j++
				local gg: word `jj' of `gmat'
				qui gen double ``v'beta'=(`v'+Dvalref*`jj')/`v' if `touse2'
			}				
								

		}
		

		markout `touse2' `datastar'

		mata: _subsmpbootddf("`data'",`ninp',`ngo',"`touse'","`datastar'","`touse2'","`gmat'",`rts',`maxiter',`tol',`brep')
	}

	
	mata:

	  dvalmat=J(N,brep,.)
	  for(i=1;i<=brep;i++){
	  	sampleid=mm_sample(ceil(M*kappa),M)
	  	dvalmat[.,i]=ddf(xyb,dataref[sampleid,.],ninp,ngo,,gmat,rts,...)
	  }


   end



	  dvalmat=J(N,brep,.)
	  for(i=1;i<=brep;i++){
	  	sampleid=mm_sample(ceil(M*kappa),M)
	  	dvalmat[.,i]=_ddf2(X0,Y0,B0,Xref[sampleid,.],Yref[sampleid,.],Bref[sampleid,.],gmat,rts,mxiter,tol)
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
