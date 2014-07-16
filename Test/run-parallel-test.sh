#!/bin/sh

DO_CHECK=1

INNAME=$1
OUTNAME=$2
OUTDIR=/tmp

PARALLELSINGLE=".parallel.single"
PARALLELMULTI=".parallel.multi"
NOCOMMENTS=".nocomments"

function run_poisson() {
	greeting=$1
	shift
   suf=$1
   shift
	echo "Running ${greeting} version"
   Bin/PoissonRecon --in "$INNAME" --out "${OUTDIR}/${OUTNAME}${suf}.ply" $@
   sed '/^comment.*$/d' "${OUTDIR}/${OUTNAME}${suf}.ply" > \
		"${OUTDIR}/${OUTNAME}${suf}${NOCOMMENTS}.ply"
}

function run_prog() {
	echo "Running for ${INNAME}"
	run_poisson "Single-threaded" $PARALLELSINGLE --threads 1 $@
	run_poisson "Multi-threaded" $PARALLELMULTI --threads 4 $@
	cmp "${OUTDIR}/${OUTNAME}${PARALLELSINGLE}${NOCOMMENTS}.ply" \
		"${OUTDIR}/${OUTNAME}${PARALLELMULTI}${NOCOMMENTS}.ply"
	res=$?
   if [[ $DO_CHECK -ne 0 ]] && [[ $res -ne 0 ]]; then
      exit $res
   fi
}

run_prog --depth 10
