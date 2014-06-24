#!/bin/sh

DO_RUN=0

OUTNAME=$1
ORIGSUF=$2
ORIGNAME=${OUTNAME}${ORIGSUF}
OUTDIR=/tmp
ORIGDIR=Examples

UNSCREENED=".unscreened"
SCREENED=".screened"
SCREENEDDENS=".screened.density"
TRIMMED=".screened.trimmed"
TRIMMEDFILLED=".screened.trimmed.filled"

function do_copy() {
	suf=$1
	ext=$2
	namefrom="${OUTDIR}/${OUTNAME}$suf.$ext"
	nameto="${ORIGDIR}/${ORIGNAME}$suf.$ext"
	if [[ $DO_RUN -ne 0 ]]; then
		cp $namefrom $nameto
	else
		echo "Will copy $namefrom --> $nameto"
	fi
}

do_copy "${UNSCREENED}" ply
do_copy "${SCREENED}" ply
do_copy "${SCREENEDDENS}" ply
do_copy "${TRIMMED}" ply
do_copy "${TRIMMEDFILLED}" ply
