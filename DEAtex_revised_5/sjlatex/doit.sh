#!/bin/sh
# doit.sh
# version 1.3.0  21dec2005#
# Written by:   Jeff Pitblado
# Description:  shell script to run -latex- on a stata journal insert

self=`basename $0`

# Subroutines ###############################################################

USAGE () {
	cat << EOF 1>&2
Usage: $self [filename]

"filename" is the name of a LaTeX file (the '.tex' extension is optional),
assumed to be 'main.tex' if one is not specified.

EOF
	exit 1
}

RUNTEX () {
	if ! latex --interaction batchmode $MAIN
	then
		echo
		echo "$self: error occured with (latex $MAIN)"
		echo "$self: see $MAIN.log"
		exit 1
	fi
}

# Parsing script args #######################################################

if [ ${#} -gt 1 ]
then
	USAGE
fi

case "$1" in
"")
	MAIN=main
	;;
-h|--help)
	USAGE
	;;
*)
	MAIN=`basename $1 .tex`
	;;
esac

if [ ! -f "$MAIN.tex" ]
then
	echo "$self: $MAIN.tex not found"
	echo
	USAGE
fi

# Document production #######################################################

rm -f *.aux *.lof *.lot *.toc *.mtc* *.bbl *.blg

RUNTEX
echo "Running: bibtex ..."
bibtex $MAIN > bibtex.log 2>&1
RUNTEX
RUNTEX

echo "$self:  ...finished"
echo "$self:  see bibtex.log for possible errors in references"
echo "$self:  ready for viewing: $MAIN.dvi"

# end
