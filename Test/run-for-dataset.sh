#!/bin/sh

DO_CHECK=1

INNAME=$1
OUTNAME=$2
ORIGSUF=$3
ORIGNAME=${OUTNAME}${ORIGSUF}
OUTDIR=/tmp
ORIGDIR=Examples

UNSCREENED=".unscreened"
SCREENED=".screened"
SCREENEDDENS=".screened.density"
TRIMMED=".screened.trimmed"
TRIMMEDFILLED=".screened.trimmed.filled"

echo "Running for $INNAME"

function run() {
   name=$1
   shift
   inp=$1
   shift
   suf=$1
   shift
   "$name" --in "$inp" --out "${OUTDIR}/${OUTNAME}$suf.ply" $@
   sed '/^comment.*$/d' "${ORIGDIR}/${ORIGNAME}$suf.ply" > "${OUTDIR}/orig"
   sed '/^comment.*$/d' "${OUTDIR}/${OUTNAME}$suf.ply" > "${OUTDIR}/new"
   cmp "${OUTDIR}/orig" "${OUTDIR}/new"
   if [[ $? -ne 0 ]]; then
      return 1
   fi
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

function run_prog() {
   prog=$1
   shift
   greeting=$1
   shift
   echo "$greeting"
   case "$prog" in
      "poisson")
         run_poisson $@
         res=$?
         ;;
      "trimmer")
         run_trimmer $@
         res=$?
         ;;
      *)
         echo "Unrecognised prog: $prog"
         exit 255
         ;;
   esac
   if [[ $DO_CHECK -ne 0 ]] && [[ $res -ne 0 ]]; then
      exit $res
   fi
}

run_prog poisson "Unscreened" "${UNSCREENED}" --depth 10 --pointWeight 0

run_prog poisson "Screened" "${SCREENED}" --depth 10

run_prog poisson "Screened with density" "${SCREENEDDENS}" --depth 10 --density

run_prog trimmer "Trimming surface" "${TRIMMED}" --trim 7 --aRatio 0

run_prog trimmer "Trimming surface and closing holes" "${TRIMMEDFILLED}" --trim 7
