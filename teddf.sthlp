{smcl}
{* *! version 0.3  29 Oct 2019}{...}
{cmd:help teddf}
{hline}

{title:Title}

{phang}
{bf:teddf} {hline 2} Directional Distance Function with undesirable outputs for Efficiency Measurement in Stata

{title:Syntax}

{p 8 21 2}
{cmd:teddf} {it:{help varlist:inputvars}} {cmd:=} {it:{help varlist:desirable_outputvars}} {cmd::} {it:{help varlist:undesirable_outputvars}} {ifin} 
{cmd:,} {cmdab:d:mu(}{varname}{cmd:)} [{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{cmdab:d:mu:(varname)}}specifies names of DMUs. It is required. 

{synopt:{cmdab:t:ime:(varname)}}specifies time period for contemporaneous production technology. If {opt time:(varname)} is not specified, global production technology is assumed. 
{p_end}

{synopt:{opt gx(varlist)}}specifies direction components for input adjustment. The order of variables specified in gx() should as the same in {it:{help varlist:inputvars}}. The default is gx=(0,..,0).
{p_end}

{synopt:{opt gy(varlist)}}specifies direction components for desirable output adjustment. The order of variables specified in gy() should as the same in {it:{help varlist:desirable_outputvars}}. The default is gy=Y.
{p_end}

{synopt:{opt gb(varlist)}}specifies direction components for undesirable output adjustment. The order of variables specified in gb() should as the same in {it:{help varlist:undesirable_outputvars}}. The default is gb=-B. 
{p_end}

{synopt:{cmdab:nonr:adial}}specifies using non-radial directional distance function.
{p_end}

{synopt:*{opt wmat(name)}}specifies a weight matrix for adjustment of input and output variables. The default is  W=(1,...,1).
{p_end}

{synopt:{opt vrs}}specifies production technology with variable returns to scale. By default, production technology with constant returns to scale is assumed.
{p_end}

{synopt:{opt rf(varname)}}specifies the indicator variable that defines which data points of outputs and inputs form the technology reference set.
{p_end}

{synopt:{opt tone}}specifies Tone's (2004) SBM model.
{p_end}

{synopt:{cmdab:win:dow(#)}}specifies window production technology with the #-period bandwidth.
{p_end}

{synopt:{cmdab:bi:ennial}}specifies biennial production technology.
{p_end}

{synopt:{cmdab:seq:uential}}specifies sequential production technology.
{p_end}

{synopt:{cmdab:glo:bal}}specifies global production technology.
{p_end}

{synopt:{opt sav:ing(filename[,replace])}}specifies that the results be saved in {it:filename}.dta. 
{p_end}

{synopt:{opt maxiter(#)}}specifies the maximum number of iterations, which must be an integer greater than 0. The default value of maxiter is 16000.
{p_end}

{synopt:{opt tol(real)}}specifies the convergence-criterion tolerance, which must be greater than 0.  The default value of tol is 1e-8.
{p_end}

{synoptline}
{p2colreset}{...}
{p 4 6 2}* wmat(name) can only be used when nonradial is specified.{p_end}




{title:Description}

{pstd}
{cmd:teddf} selects the input and output variables in the opened data set and solves directional 
distance function models by options specified. 

{phang}
The teddf program uses the buit-in mata function linearprogram(). Stata 16 or later is required.

{phang}
The teddf program requires initial data set that contains the input and output variables for observed units. 

{phang}
Variable names must be identified by inputvars for input variable, by desirable_outputvars for desirable output variable,  and by undesirable_outputvars for undesirable output variable
 to allow that {cmd:teddf} program can identify and handle the multiple input-output data set.



{title:Examples}

{phang}{"use ...\example_ddf.dta"}

{phang}{cmd:. teddf labor capital energy= gdp: co2, dmu(id)}

{phang}{cmd:. teddf labor capital energy= gdp: co2, dmu(id) time(t) nonr sav(ddf_result)}

{phang}{cmd:. teddf labor capital energy= gdp: co2, dmu(id) nonr vrs sav(ddf_result,replace)}

{phang}{cmd:. teddf labor capital energy= gdp: co2, dmu(id) time(t) seq sav(ddf_result,replace)}

{title:Saved Results}

{psee}
Macro:

{psee}
{cmd: r(file)} the filename stores results of {cmd:teddf}.
{p_end}


{marker references}{...}
{title:References}

{phang}
Chung, Y.H., Färe, R., Grosskopf, S. Productivity and Undesirable Outputs: A Directional Distance Function Approach.
 Journal of Environmental Management, 1997, 51:229-240.
 
{phang}
Färe, R., Grosskopf, S. Directional distance functions and slacks-based measures of efficiency. European Journal of Operational Research, 2010, 200:320-322.

{phang}
Tone, K. (2004). Dealing with Undesirable Outputs in DEA: A Slacks-Based Measure (SBM) Approach. North American Productivity Workshop 2004, Toronto, 23-25 June 2004, 44-45.

{phang}
Oh, D.-h. A global Malmquist-Luenberger productivity index. Journal of Productivity Analysis, 2010, 34:183-197.

{phang}
Oh, D.-h., Heshmati A.  A sequential Malmquist–Luenberger productivity index: Environmentally sensitive productivity growth 
considering the progressive nature of technology. Energy Economics, 2010,3 2:1345-1355.

{phang}
Zhou, P., Ang, B.W., Wang, H. Energy and CO2 emission performance in electricity generation: a non-radial directional distance function approach. Eur. J. Oper. Res., 2012, 221:625-635.

{title:Author}

{psee}
Kerry Du

{psee}
Xiamen University

{psee}
Xiamen, China

{psee}
E-mail: kerrydu@xmu.edu.cn
{p_end}


{psee}
Daoping Wang

{psee}
Shanghai University of Finance and Economics

{psee}
Shanghai, China

{psee}
E-mail: daopingwang@live.sufe.edu.cn
{p_end}


{psee}
Ning Zhang

{psee}
Jinan University

{psee}
Guangzhou, China

{psee}
E-mail: zn928@naver.com
{p_end}
