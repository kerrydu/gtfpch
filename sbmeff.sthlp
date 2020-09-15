{smcl}
{* *! version 0.3  29 Oct 2019}{...}
{cmd:help sbmeff}
{hline}

{title:Title}

{phang}
{bf:sbmeff} {hline 2} Slacks-based Measure of Efficiency in Stata

{title:Syntax}

{p 8 21 2}
{cmd:sbmeff} {it:{help varlist:inputvars}} {cmd:=} {it:{help varlist:desirable_outputvars}} {ifin} 
{cmd:,} {cmdab:d:mu(}{varname}{cmd:)} [{it:options}]

{p 8 21 2}
{cmd:sbmeff} {it:{help varlist:inputvars}} {cmd:=} {it:{help varlist:desirable_outputvars}} {cmd::} {it:{help varlist:undesirable_outputvars}} {ifin} 
{cmd:,} {cmdab:d:mu(}{varname}{cmd:)} [{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{cmdab:d:mu:(varname)}}specifies names of DMUs. It is required. 

{synopt:{cmdab:t:ime:(varname)}}specifies time period for contemporaneous production technology. If {opt time:(varname)} is not specified, global production technology is assumed. 
{p_end}

{synopt:{cmdab:win:dow(#)}}specifies window production technology with the #-period bandwidth.
{p_end}

{synopt:{cmdab:bi:ennial}}specifies biennial production technology.
{p_end}

{synopt:{cmdab:seq:uential}}specifies sequential production technology.
{p_end}

{synopt:{cmdab:glo:bal}}specifies global production technology.
{p_end}

{synopt:{opt vrs}}specifies production technology with variable returns to scale. By default, production technology with constant returns to scale is assumed.
{p_end}

{synopt:{opt rf(varname)}}specifies the indicator variable that defines which data points of outputs and inputs form the technology reference set.
{p_end}

{synopt:{opt sav:ing(filename[,replace])}}specifies that the results be saved in {it:filename}.dta. 
{p_end}

{synopt:{opt maxiter(#)}}specifies the maximum number of iterations, which must be an integer greater than 0. The default value of maxiter is 16000.
{p_end}

{synopt:{opt tol(real)}}specifies the convergence-criterion tolerance, which must be greater than 0.  The default value of tol is 1e-8.
{p_end}

{synoptline}
{p2colreset}{...}

{title:Description}

{pstd}
{cmd:sbmeff} selects the input and output variables from the user designated data file or in the opened data set and solves Slacks-based Measure of Efficiency models by options specified. 

{phang}
The sbmeff program uses the buit-in mata function linearprogram(). Stata 16 or later is required.

{phang}
The sbmeff program requires initial data set that contains the input and output variables for observed units. 

{phang}
Variable names must be identified by inputvars for input variable, by desirable_outputvars for desirable output variable,  and by undesirable_outputvars for undesirable output variable
 to allow that {cmd:sbmeff} program can identify and handle the multiple input-output data set.



{title:Examples}

{phang}{"use ...\example_sbm.dta"}

{phang}{cmd:. sbmeff labor capital energy=gdp, dmu(id) time(t) vrs}

{phang}{cmd:. sbmeff labor capital energy=gdp, dmu(id) sav(sbm_result)}

{phang}{cmd:. sbmeff labor capital energy=gdp:co2, dmu(id) time(t) vrs}

{phang}{cmd:. sbmeff labor capital energy=gdp:co2, dmu(id) sav(sbm_result,replace)}


{title:Saved Results}

{psee}
Macro:

{psee}
{cmd: r(file)} the stored results of {cmd:sbmeff} that have observation rows of DMUs and variable columns with input data, output data, efficiency scores, and slacks.
{p_end}


{marker references}{...}
{title:References}

{phang}
Tone, K. (2001). A slacks-based measure of efficiency in data envelopment analysis. European Journal of Operational Research, 130(3):498-509.

{phang}
Tone, K. (2004). Dealing with Undesirable Outputs in DEA: A Slacks-Based Measure (SBM) Approach. North American Productivity Workshop 2004, Toronto, 23-25 June 2004, 44-45.


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
