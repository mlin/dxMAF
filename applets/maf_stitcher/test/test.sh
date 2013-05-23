#!/bin/bash -e

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ] ; do SOURCE="$(readlink "$SOURCE")"; done
HERE="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

function test_in_project {
  dx-build-applet $HERE/.. -d $TEST_PROJECT:/
  dx="dx --project-context-id=$TEST_PROJECT"
  job2=$($dx run :/maf_stitcher -i species_set=dm3_12flies \
          -i regions="MAF whole-genome alignments:/zzz data for testing/tal-AA" \
          -i exon_type=CDS -y --brief)
  job=$($dx run :/maf_stitcher -i species_set=hg19_33mammals \
          -i regions="MAF whole-genome alignments:/zzz data for testing/ccds.hg19.chr21" \
          -i transcripts_per_job=50 -i exon_type=CDS -y --brief)
  curl -sL ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs103/CCDS_nucleotide.20121025.fna.gz > $HERE/ccds.fna.gz &
  curlpid=$!
  $dx wait $job2
  $dx find jobs --origin $job2
  dna=$($dx head $job2:alignments -w 1024 | awk 'BEGIN { FS="│"; } NR==4 {print $5;}' | tr -d '-')
  talAA="ATGCTGGATCCCACTGGAACATACCGGCGACCACGCGACACGCAGGACTCCCGCCAAAAGAGGCGACAGGACTGCCTGGATCCAACCGGGCAGTACTAG"
  if [ "$dna" != "$talAA" ]; then
    echo 'tal-AA mismatch!'
    echo $dna
    echo $talAA
    exit 1
  fi
  dx wait $job
  $dx find jobs --origin $job
  wait $curlpid
  $HERE/verify.py $job $HERE/ccds.fna.gz
}

if [ "$1" = "inner" ]; then
  test_in_project
  exit
fi

TEST_PROJECT=$(dx new project --brief "maf_stitcher test ($(date))")
export TEST_PROJECT

function cleanup {
  dx rmproject -y $TEST_PROJECT > /dev/null
}
trap "{ cleanup; exit }" int

$SOURCE inner &
inner=$!

exit_code=0
wait $inner || exit_code=$?

if [ "$exit_code" -ne 0 ]; then
  echo "Leaving behind project for failed test, named \"$(dx describe $TEST_PROJECT --name)\""
  echo "To delete: dx rmproject -y $TEST_PROJECT"
else
  cleanup
fi

exit $exit_code
