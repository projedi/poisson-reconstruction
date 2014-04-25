#!/bin/sh

INNAME=$1
OUTNAME=$2
ORIGSUF=$3
ORIGNAME=${OUTNAME}${ORIGSUF}
OUTDIR=/tmp
ORIGDIR=Examples

UNSCREENED=".unscreened.ply"
SCREENED=".screened.ply"
SCREENEDDENS=".screened.density.ply"
TRIMMED=".screened.trimmed.ply"
TRIMMEDFILLED=".screened.trimmed.filled.ply"

echo "Running for $INNAME"

function run() {
   name=$1
   shift
   inp=$1
   shift
   suf=$1
   shift
   "$name" --in "$inp" --out "${OUTDIR}/${OUTNAME}$suf" $@
   sed '/^comment.*$/d' "${ORIGDIR}/${ORIGNAME}$suf" > "${OUTDIR}/orig"
   sed '/^comment.*$/d' "${OUTDIR}/${OUTNAME}$suf" > "${OUTDIR}/new"
   cmp "${OUTDIR}/orig" "${OUTDIR}/new"
}

function run_poisson() {
   suf=$1
   shift
   # Parallelism is undeterministic here. Therefore it must be run on a single thread.
   run "Bin/PoissonRecon" "$INNAME" "$suf" --threads 1 $@
}

function run_trimmer() {
   suf=$1
   shift
   run "Bin/SurfaceTrimmer" "${OUTDIR}/${OUTNAME}${SCREENEDDENS}" "$suf" $@
}

echo "Unscreened"
run_poisson "${UNSCREENED}" --depth 10 --pointWeight 0

echo "Screened"
run_poisson "${SCREENED}" --depth 10

echo "Screened with density"
run_poisson "${SCREENEDDENS}" --depth 10 --density

echo "Trimming surface"
run_trimmer "${TRIMMED}" --trim 7 --aRatio 0

echo "Trimming surface and closing holes"
run_trimmer "${TRIMMEDFILLED}" --trim 7
