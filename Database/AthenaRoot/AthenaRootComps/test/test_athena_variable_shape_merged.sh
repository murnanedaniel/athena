#!/bin/bash

set -e

reffile=$1

testname=test_athena_variable_shape
mkdir -p $testname
cd $testname
rm -f *.txt *.xml

echo "::: generate merged.root..."
athena.py \
    -c 'EVTMAX=30; BRANCHES=["RunNumber", "EventNumber", "el_n", "el_eta","el_phi"]; OUTBRANCHES=["el_n",]; FNAMES=["shape1/f1.root","shape2/f2.root","shape3/f3.root"]' \
    AthenaRootComps/test_athena_variable_shape_ntuple.py \
    >| log.004.txt 2>| log.004.stderr.txt \
    || exit 1
/bin/mv d3pd.root f4.root || exit 1
/bin/mv data.var.txt data.4.txt || exit 1
(acmd.py dump-root f4.root 2> f4.stderr.txt) \
    | grep "^egamma" \
    | egrep "RunNumber|EventNumber|el_n" \
    | tee f4.ascii \
    || exit 1
cat f4.ascii| cut -d. -f3 >| f4.ascii.todiff || exit 1

echo "::: compare py-alg outputs..."
cat shape[123]/data.*.txt > data.merged.txt || exit 1
diff -urN data.merged.txt data.4.txt || exit 1
echo "::: compare py-alg outputs... [ok]"

echo "::: compare py-alg output to reference..."
diff -urN data.merged.txt $reffile || exit 1
echo "::: compare py-alg output to reference... [ok]"

echo "::: compare dump-root outputs..."
cat shape[123]/f*.ascii.todiff > merged.ascii.todiff || exit 1
diff -urN merged.ascii.todiff f4.ascii.todiff || exit 1
echo "::: compare dump-root outputs... [ok]"

exit 0
