version 16

// Versioning
_get_version gtfpch
assert("`package_version'" != "")
mata: string scalar gtfpch_version() return("`package_version'")
mata: string scalar gtfpch_stata_version() return("`c(stata_version)'")
mata: string scalar gtfpch_joint_version() return("`package_version'|`c(stata_version)'")

////////////////////////////////
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
//mata clear 

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
