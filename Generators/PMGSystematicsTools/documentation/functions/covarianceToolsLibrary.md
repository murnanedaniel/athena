Help on module covarianceToolsLibrary

# NAME
covarianceToolsLibrary

# FILE
local/bin/covarianceToolsLibrary.py

# DESCRIPTION
The `covarianceToolsLibrary.py` can be used directly to manipulate covariance info and evaluate goodness of fit, inmported into your own python script.
Examples of how to do this are provided in the `examples/covarianceTool-examples/` directory

Author: Louie D. Corpe (UCL)  
Email: l.corpe@cern.ch

# FUNCTIONS
## addSystematicToCorrInfo(h, systName, binErrorsToAdd)

## calculateChi2(data, mc, cov=None, verbosity=0)

## chi2ContribsByRow(chi2contribs)

## chi2ContributionMatrix(diff, prec)

## correlation(...)
correlation(sample1, sample2)
(float, list[float]) -> int
Return the unweighted correlation of the two provided sample lists.

## correlationMatrix(covMatrix)

## covariance(...)
covariance(sample1, sample2)
(list[float], list[float]) -> float
Return the unweighted covariance of the two provided sample lists.

## divide(...)
divide(ao1, ao2)
(AnalysisObject, AnalysisObject) -> Scatter{1,2,3}D
Divide one AnalysisObject by another, producing a Scatter of appropriate dimension by using the logic of the bound divideBy methods.

## drawMatrix(matrix, outfile, xLabel=None)

## getCorrInfo(h, names=[])

## getDiagonalCovMatricesFromHisto(scat)

## getDiagonalCovMatricesFromScatter(scat)

## hasAsymmetricErrors(ao)

## histoToMatrix(h, bw=False)

## index_between(...)
index_between(x, binedges)
(float, list[float]) -> int
Return the index of the bin which would contain x, or -1 if there is no enclosing
bin in the given set of n+1 bin edges.

## linspace(...)
linspace(nbins, start, end)
(int, float, float) -> list[float]
Make a list of n+1 bin edges linearly spaced between start and end, with the first and
last edges on those boundaries.

## logspace(...)
logspace(nbins, start, end)
(int, float, float) -> list[float]
Make a list of n+1 bin edges linearly spaced on the interval log(start..end), with
the first and last edges on those boundaries.

## makeCovarianceMatrix(ao, ignore_corrs=False)

## makeCovarianceMatrixFromToys(ao, ntoys=10000, ignore_corrs=False)

## makeSuperAO(aoList, verbosity=0)

## manualChi2(data, mc, cov, verbosity=0)

## mean(...)
mean(sample)
(list[float]) -> float
Return the unweighted mean of the entries in the provided sample list.

## mkScatter(...)
mkScatter(ao, usefocus=False, usestddev=False)
AnalysisObject -> Scatter{1,2,3}D
Convert an AnalysisObject to a Scatter, using the logic of the bound mkScatter methods.

@todo This falls back on use of optional args until we find one that works: is there a nicer way?

## plotVariations(ao, outdir, label, threshold=0.05, xLabel=None, title=None)

## printMatrix(matrix, outfile)

## read(...)
read(filename, asdict=True, patterns=None, unpatterns=None)

Read data objects from the provided filename, auto-determining the format
from the file extension.

The loaded data objects can be filtered on their path strings, using the
optional patterns and unpatterns arguments. These can be strings, compiled
regex objects with a 'match' method, or any iterable of those types. If
given, only analyses with paths which match at least one pattern, and do not
match any unpatterns, will be returned.

Returns a dict or list of analysis objects depending on the asdict argument.

## readAIDA(...)
readAIDA(file_or_filename, asdict=True, patterns=None, unpatterns=None)

Read data objects from the provided AIDA-format file.

The loaded data objects can be filtered on their path strings, using the
optional patterns and unpatterns arguments. These can be strings, compiled
regex objects with a 'match' method, or any iterable of those types. If
given, only analyses with paths which match at least one pattern, and do not
match any unpatterns, will be returned.

Returns a dict or list of analysis objects depending on the asdict argument.

DEPRECATED: AIDA is a dead format. At some point we will stop supporting it.

## readFLAT(...)
readFLAT(file_or_filename, asdict=True, patterns=None, unpatterns=None)

Read data objects from the provided FLAT-format file.

The loaded data objects can be filtered on their path strings, using the
optional patterns and unpatterns arguments. These can be strings, compiled
regex objects with a 'match' method, or any iterable of those types. If
given, only analyses with paths which match at least one pattern, and do not
match any unpatterns, will be returned.

Returns a dict or list of analysis objects depending on the asdict argument.

## readYODA(...)
readYODA(file_or_filename, asdict=True, patterns=None, unpatterns=None)

Read data objects from the provided YODA-format file.

The loaded data objects can be filtered on their path strings, using the
optional patterns and unpatterns arguments. These can be strings, compiled
regex objects with a 'match' method, or any iterable of those types. If
given, only analyses with paths which match at least one pattern, and do not
match any unpatterns, will be returned.

Returns a dict or list of analysis objects depending on the asdict argument.

## scatterToMatrix(h, bw=False)

## setCorrInfo(h, binErrors)
`h` AO (the Scatter you want to set the correlation info for) 
`binErrors`

## updateAOinYODA(ao, infile)

## version(...)
version()
Return YODA library version as a string

## write(...)
write(ana_objs, filename)

Write data objects to the provided filename,
auto-determining the format from the file extension.

## writeAIDA(...)
writeAIDA(ana_objs, file_or_filename)

Write data objects to the provided file in AIDA format.

## writeFLAT(...)
writeFLAT(ana_objs, file_or_filename)

Write data objects to the provided file in FLAT format.

## writeYODA(...)
writeYODA(ana_objs, file_or_filename)

Write data objects to the provided file in YODA format.

# DATA
## HAS_ROOT_SUPPORT = False


