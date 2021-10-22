cap program drop gtfpchcheckupdate
program define gtfpchcheckupdate
version 16
global gtfpchcheck done
local url1 https://raw.githubusercontent.com/kerrydu/gtfpch/master/
local url2 https://gitee.com/kerrydu/gtfpch/raw/master/
cap mata: vfile = cat(`"`url1'/`0'.ado"')
if _rc{
	cap mata: vfile = cat(`"`url2'/`0'.ado"')
}
if _rc exit

mata: vfile = select(vfile,vfile:!="")
mata: vfile = usubinstr(vfile,char(9)," ",.)
mata: vfile = select(vfile,!ustrregexm(vfile,"^( )+$"))
mata: st_local("versiongit",vfile[1])
local versiongit = ustrregexrf("`versiongit'","^[\D]+","")
gettoken vers versiongit:versiongit, p(", ")
local versiongit `vers'

qui findfile `0'.ado
mata: vfile = cat("`r(fn)'")
mata: vfile = select(vfile,vfile:!="")
mata: vfile = usubinstr(vfile,char(9)," ",.)
mata: vfile = select(vfile,!ustrregexm(vfile,"^( )+$"))
mata: st_local("versionuse",vfile[1])
local versionuse = ustrregexrf("`versionuse'","^[\D]+","")
gettoken vers versionuse:versionuse, p(", ")
local versionuse `vers'	

if(`versionuse'<`versiongit'){
	di "New version available, `versionuse' =>`versiongit'"
	di "It can be updated by:"
	di "  net install gtfpch,from(https://raw.githubusercontent.com/kerrydu/gtfpch/master/) replace"
	di "or,"
	di "  net install gtfpch,from(https://gitee.com/kerrydu/gtfpch/raw/master/) replace"
}

end
