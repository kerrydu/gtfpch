@echo off
rem DOS batch file to LaTeX a Stata Journal insert
rem version 2.0.1  21dec2005

set PS=ps
set MAIN=main
:beginloop
if "%1" EQU "alt" set MAIN=altmain
if "%1" EQU "altmain" set MAIN=altmain
if "%1" EQU "altmain.tex" set MAIN=altmain
if "%1" EQU "nops" set PS=nops
if "%1" EQU "?" goto help
if "%1" EQU "h" goto help
if "%1" EQU "-h" goto help
if "%1" EQU "help" goto help
if "%1" EQU "--help" goto help
shift
if "%1" NEQ "" goto beginloop
:endloop

goto skiphelp

:help
echo usage: %0 [option(s)]
echo.
echo Run latex and bibtex for a Stata Journal insert.
echo.
echo %0 will run bibtex on main.aux.  A log of the results are placed in
echo bibtex.log
echo.
echo Options:
echo         help	Display this information and quit
echo         altmain	latex altmain.tex instead of main.tex
echo         nops	do not generate a PostScript file

goto exit

:skiphelp

latex %MAIN%
if errorlevel 1 goto latexerror
echo.
bibtex %MAIN%
echo.
latex %MAIN%
if errorlevel 1 goto latexerror
echo.
latex %MAIN%
if errorlevel 1 goto latexerror
echo.

if "%PS%" EQU "nops" goto dvimsg
echo Generating PostScript file:
dvips %MAIN%.dvi -o %MAIN%.ps
if errorlevel 1 goto dvipserror
echo.
echo finished...
echo you may new view %MAIN%.ps

goto exit

:dvimsg
echo finished...
echo you may now view %MAIN%.dvi

goto exit

:latexerror
echo error occured whith (latex %MAIN%)
echo exiting now
goto exit

:dvipserror
echo error when running dvips
goto exit

:exit
