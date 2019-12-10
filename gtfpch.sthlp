{smcl}
{* *! version 0.3  29 Oct 2019}{...}
{cmd:help gtfpch}
{hline}

{title:Title}

{phang}
{bf:gtfpch} {hline 2} Total Factor Productivity with Undesirable Outputs in Stata 

{title:Syntax}

{p 8 21 2}
{cmd:gtfpch} {it:{help varlist:inputvars}} {cmd:=} {it:{help varlist:desirable_outputvars}} {cmd::} {it:{help varlist:undesirable_outputvars}} {ifin} 
{cmd:,}  [{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{cmdab:d:mu:(varname)}}specifies names of DMUs. 

{synopt:{opt gx(varlist)}}specifies direction components for input adjustment. The order of variables specified in gx() should as the same in {it:{help varlist:inputvars}}. 
{p_end}

{synopt:{opt gy(varlist)}}specifies direction components for desirable output adjustment. The order of variables specified in gy() should as the same in {it:{help varlist:desirable_outputvars}}. 
{p_end}

{synopt:{opt gb(varlist)}}specifies direction components for undesirable output adjustment. The order of variables specified in gb() should as the same in {it:{help varlist:undesirable_outputvars}}.  
{p_end}

{synopt:{cmdab:seq:uential}}specifies sequential production technology.
{p_end}

{synopt:{opt global}}specifies global production technology.
{p_end}

{synopt:{cmdab:nonr:adial}}specifies using non-radial directional distance function.
{p_end}

{synopt:*{opt wmat(name)}}specifies a weight matrix for adjustment of input and output variables. 
{p_end}

{synopt:*{cmdab:luen:berger}}specifies estimating Luenberger productivity index. The default is Malmquist–Luenberger productivity index.
{p_end}

{synopt:{opt ort(string)}}specifies the oriention. The default is ort(out), 
meaning the output oriented productivity index. ort(i) means the input oriented productivity index.
{p_end}

{synopt:{opt fgnz}}specifies decomposing TFP change following the spirit of Färe, Grosskopf, Norris, and Zhang's (1994).
{p_end}

{synopt:{opt rd}}specifies decomposing TFP change following the spirit of Ray and Desli's (1997) method.
{p_end}

{synopt:{opt sav:ing(filename[,replace])}}specifies that the results be saved in {it:filename}.dta. 
{p_end}

{synopt:{opt maxiter(#)}}specifies the maximum number of iterations, which must be an integer greater than 0. The default value of maxiter is 16000.
{p_end}

{synopt:{opt tol(real)}}specifies the convergence-criterion tolerance, which must be greater than 0.  The default value of tol is 1e-8.
{p_end}

{synoptline}
{p2colreset}{...}
{p 4 6 2}  A panel variable and a time variable must be specified; use xtset. {p_end}
{p 4 6 2}* wmat(name) can only be used when nonradial is specified. {p_end}
{p 4 6 2}* Luenberger productivity index is estimated when nonradial is specified.{p_end}




{title:Description}

{pstd}
{cmd:gtfpch} selects the input and output variables in the opened data set and estimate the TFP index by options specified. 

{phang}
The gtfpch program requires initial data set that contains the input and output variables for observed units. 

{phang}
Variable names must be identified by inputvars for input variable, by desirable_outputvars for desirable output variable,  and by undesirable_outputvars for undesirable output variable
 to allow that {cmd:gtfpch} program can identify and handle the multiple input-output data set.



{title:Examples}

{phang}{cmd:. use "https://raw.githubusercontent.com/kerrydu/gtfpch/master/example_ddf.dta"}

{phang}{cmd:. xtset id t}

{phang}{cmd:. gtfpch labor capital energy= gdp: co2,  nonr sav(ddf_result)}

{phang}{cmd:. gtfpch labor capital energy= gdp: co2,  ort(i)}

{phang}{cmd:. gtfpch labor capital energy= gdp: co2, seq luen sav(ddf_result,replace)}

{title:Saved Results}

{psee}
Macro:

{psee}
{cmd: r(file)} the filename stores results of {cmd:gtfpch}.
{p_end}


{marker references}{...}
{title:References}

{phang}
Chung, Y.H., Färe, R., Grosskopf, S. Productivity and Undesirable Outputs: A Directional Distance Function Approach.
 Journal of Environmental Management, 1997, 51:229-240.
 
{phang}
Färe, R., Grosskopf, S. Directional distance functions and slacks-based measures of efficiency. European Journal of Operational Research, 2010, 200:320-322.

{phang}
Oh, D.-h. A global Malmquist-Luenberger productivity index. Journal of Productivity Analysis, 2010, 34:183-197.

{phang}
Oh, D.-h., Heshmati A.  A sequential Malmquist–Luenberger productivity index: Environmentally sensitive productivity growth 
considering the progressive nature of technology. Energy Economics, 2010,3 2:1345-1355.

{phang}
Zhou, P., Ang, B.W., Wang, H. Energy and CO2 emission performance in electricity generation: a non-radial directional distance function approach. Eur. J. Oper. Res., 2012, 221:625-635.

{phang}
Mahlberg, B., Sahoo, B.K. Radial and non-radial decompositions of Luenberger productivity indicator with an illustrative application,
International Journal of Production Economics, 2011, 131:721-726.

{phang}
Färe, R., Grosskopf, S., Norris, M., & Zhang, Z. (1994). Productivity Growth, Technical Progress, and Efficiency Change in Industrialized Countries. The American Economic Review, 84(1), 66-83.

{phang}
Ray, S., & Desli, E. (1997). Productivity Growth, Technical Progress, and Efficiency Change in Industrialized Countries: Comment. The American Economic Review, 87(5), 1033-1039.


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
