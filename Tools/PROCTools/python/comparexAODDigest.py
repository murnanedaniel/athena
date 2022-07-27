#!/usr/bin/env python

# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
import argparse
import sys

def extractData(filename):

    result=dict()

    rein=open(filename)
    for line in rein:
        l=line.split()
        if l[0]=='run':
            header=l
        else:
            run=int(l[0])
            evt=int(l[1])
            result[(run,evt)]=[float(s) for s in l[2:]]
    rein.close()
    return result,header


def compare2Files(file1, file2, ignoredColumns=None):

    res1,h1=extractData(file1)
    res2,h2=extractData(file2)


    if (h1 != h2):
        print("ERROR, headers don't match")
        print(h1)
        print(h2)
        return 1

    runEvts = res1.keys()
    if len(res2) < len(res1):
        runEvts = res2.keys()

    diffCounter=dict()
    for h in h1[2:]:
        diffCounter[h]=0

    # Loop over events:
    for runEvt in runEvts:
        values1 = res1[runEvt]
        values2 = res2[runEvt]
        for i, name in enumerate(h1[2:]):
            #print (name,i,len(values1),len(values2))
            if values1[i] != values2[i]:
                ignored = False
                if ignoredColumns and name in ignoredColumns:
                    ignored = True
                suffix = ' (ignored)' if ignored else ''
                print ("Diff: Run {} Evt {} {} {} -> {}{}".format(runEvt[0],runEvt[1],name,values1[i],values2[i],suffix))
                diffCounter[name]+=1
                pass
            pass
        pass

    print("Summary of differences:")
    noChanges=""
    nEvt=len(runEvts)
    changed = False
    for (name,count) in diffCounter.items():
        if (count>0):
            #print (name,":",count,"(of ",nEvt,")")
            print ("{}: {} events (out of {})".format(name,count,nEvt))
            if not ignoredColumns or name not in ignoredColumns:
                changed = True
        else:
            noChanges+=" "+name
    if not changed:
        print("No changes")
    else:
        print("No changes for:",noChanges)
    return changed


def compareDigest(filelist):
    if len(filelist)<2:
        print("Got only %i files. Can't compare")
        return None

    runevtset=set()

    summary=dict() #key is the datestamp

    for f in filelist:
        datestamp=f.split('/')[9]
        print("Fond file for %s" % datestamp)
        header=None
        if datestamp in summary:
            print("ERROR, duplicate date-stamp %s" % datestamp)
            continue

        res,hdr=extractData(f)
        if header is None:
            header=hdr
        elif (header!=hdr):
            print("ERROR, headers of file %s doesn't match!" % f)
            continue

        summary[datestamp]=res
        runevtset |= set(res.keys())
        pass


    #collected all data, now invert the matrix, sort per run/event
    nValues=len(header)-2

    perEvt=dict()
    for runevt in runevtset:
        perEvt[runevt]=[]
        for i in range(nValues):
            perEvt[runevt].append(set())

    for day,data in summary.items():
        for runevt,v in data.items():
            for i in range(nValues):
                perEvt[runevt][i].add(v[i])


    row_format ="{:>12}" * len(header)
    #row_format+=os.linesep
    print (row_format.format(*header))
    for runevt,v in perEvt.items():
        updates=[runevt[0],runevt[1]]
        updates+=[len(x)-1 for x in v]
        print (row_format.format(*updates))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A script to compare 2 xAODDigest files")
    parser.add_argument("digestFile", nargs=2, type=str,
                        help="digest filename", action="store")
    parser.add_argument("--ignoreMuons", help="Ignore muons",
                        action="store_true", default=False)
    parser.add_argument("--ignoreMET", help="Ignore MET",
                        action="store_true", default=False)

    args = parser.parse_args()

    ignored = []
    if args.ignoreMuons:
        ignored += ['nMuons', 'muon1pt', 'muon1eta', 'muon1phi']
    if args.ignoreMET:
        ignored += ['nmet', 'metx', 'mety', 'sumet']

    retval = compare2Files(args.digestFile[0], args.digestFile[1], ignoredColumns=ignored)

    sys.exit(retval)
