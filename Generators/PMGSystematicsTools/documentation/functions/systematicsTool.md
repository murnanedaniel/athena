Help on module systematicsTool

# NAME
systematicsTool

# FILE
local/bin/systematicsTool.py

# DESCRIPTION
This module contains helper functions which allow users to manipluate multiple 
YODA files corresponding to different subsamples of a given process/generator,
different sources of theory uncertainty, etc. and combine them, plotting the 
theory error bands along the way. 

Can be used either as a standalone exectuable or the functions can be imported
intot a custom python script.

This module uses a loosely-defined datatype which shall be referred to as 
AnalysisObject, which is effectively a dict of np.arrays and strings.
This encodes the information from a yoda.Scatter2D(), yoda.Histo1D(),
root.TGraphAsymmErrors() or root.TH1D() in the same format, and which
can be easily manipulated since it consists of np.arrays or strings
An AO has the following keys:

`ao['name']`   the name of the object a la YODA   
`ao['path']`   the path of the object a la YODA  
`ao['annotations']`  a dict of annotations a la YODA  
`ao['bw']`     np.array of widths of each bin  
`ao['x']`      np.array of mid-points of bins  
`ao['xup']`    np.array of high bin edges  
`ao['xdn']`    np.array of low bin edges  
`ao['y']`      np.array of y values  
`ao['yup']`    np.array of y+err_up values  
`ao['ydn']`    np.array of y-err_dn values  
`ao['nbins']`  number of bins  

This should probably be replaced by a Class in the next version !

For a full list of functions and their usage, type:
pydoc systematicsTool

TODO (perspectives for future development):
- Also handle ROOT as input/output instead of YODA
- Replace AnalysisObject by an actual Class rather than loosely-defined dicts
- Implement combineVariationsReplicas()
- Implement includeErrorsInEnvelope for combineVariationsEnvelope() 
- Support other types in readFromYODA


Author: Louie D. Corpe (CERN)  
Email: l.corpe@cern.ch

# FUNCTIONS
## arrayDictToTGraph(ao, isData=False, setYErrorsToZero=False, nominalAOForRatio=None)
`ao` AnalysisObject (the one to convert to a TGraph)  
`isData` Bool [optional] (modify style of TGraph for data points)  
`setYErrorsToZero` Bool [optional](Might want to do this if you don't care about the
error bands for some uncertainty component)  
`nominalAOForRatio` [optional]AnalysisObject or None (if not None, make a ratio plot
with respect to this AnalysisObject)  
`return` ROOT.TGraphAsymmErrors()  

Fill a TGraphAsymmErrors from an AO, for plotting purposes. 
Can format as MC or Data. If nominalAOForRatio is specified (ie not None), 
the TGraphAsymmErrors is divided by the 'y' value of the nominalAOForRatio 
in each bin.

## combineAllVariations(weightList, indir, outdir, regexFilter=None, regexVeto=None, combinationRecipe=None, returnOnlyVariationsInComination=True, schema='!INDIR/!WEIGHTNAME.yoda:!AONAME', orderedWeights=None, inFile=None)
`weightList` dict in output format of the readDatabase.getWeight() tool.
See https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/PmgSystTool
for more info.  
`indir` String (input dir, where to find samples for each weight (merged
across jobs/DSIDs if needed)  
`outdir` String (output dir, where to write the output files where the
individual weights are combined into systematic variations)  
`regexFilter` String [optional] (AOs whose names match regex are processed)  
`regexVeto` String [optional](AOs whose names match regex are NOT processed) 
`combinationRecipe` String [optional] (if you want to use a specific 
combination recipe, specify it here. If None, then the code will try to
auto-determine the correct recipe from the $SYSTTOOLSPATH/data/Syst_Database.yaml.
supports String:String where the first part is the yaml file containing the recipe)
`returnOnlyVariationsInComination` Bool [optional] (by default the code will return 
only a list of files for variations which take part in the combination. Set to False
to return all possible variations)
`schema` String [optional] (tells the code how the naming convention for histograms is  
structured within the input file.)
`schema` List[weight names] [optional] (if the `!WEIGHTINDEX!` keyword is used in `schema`  
this option is needed so the code can convert from weight name to index)
`inFile` String [optional] (needed if the `!INFILE!` keyword is used in `schema`)
`return` dict{String,String} (format is as follows {N!label:filename},
where N is the order in which to plot the variations, labels is what to put on
the legend, and filename is the name of the YODA/ROOT file to use as input)  

Using the information provided by the readDatabase tool (see 
https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/PmgSystTool),
combine individual YODA/ROOT files for each Matrix Element weights into
individual YODA/ROOT files for each systematic uncertainty.

## combineVariation(wName, wInfo, fOut, regexFilter=None, regexVeto=None)
`wName` String (Name of the systemtic weight to combine)  
`wInfo` dict in output format of the readDatabase.getWeight() tool.
See https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/PmgSystTool
for more info.  
`regexFilter` String [optional] (AOs whose names match regex are processed)  
`regexVeto` String [optional] (AOs whose names match regex are NOT processed)  
`return` None  

Produce aand write a YODA file for the systematioc uncertainty wName by 
combining the weights listed in wInfo['weights'], according to the procedure
specified in wInfo['combination'] for each AO.

The following values of wInfo['combination'] are suppored: stat,lhapdf,
replicas, alphas, envelope, hessian, and customFunction ( in this case, the
formula given in wInfo['function'] is evaluated).

## combineVariationsAlphaS(nom, variations)
`nom` AnalysisObject (of the nominal variation)  
`variations` dict(String, AnalysisObject)  
`return` AnalysisObject  

Get the alphaS variation. This is just a special case of Envelope
where the error is symmetrized!

## combineVariationsEnvelope(nom, variations, asym=True, includeErrorsInEnvelope=False)
`nom` AnalysisObject (of the nominal variation)  
`variations` dict(String, AnalysisObject)  
`asym` Bool [optional] (return asymmetric errors? Or symmetrize?)  
`includeErrorsInEnvelope` Bool [optional] (TODO to be implemented)  
`return` AnalysisObject  

Take the min/max envelope of the arguments as the up/down errors of the
resulting AO. The central value is taken from 'nominal'.

## combineVariationsFromFormula(nom, variations, formula)
`nom` AnalysisObject (of the nominal variation)  
`variations` dict(String, AnalysisObject)  
`formula` String (custom formula to combine uncertainties)  
`return` AnalysisObject  

Combine the specified variations according to custom formula given by formula

## combineVariationsHessian(nom, variations, asym=True)
`nom` AnalysisObject (of the nominal variation)  
`variations` dict(String, AnalysisObject)  
`asym` Bool [optional] (return asymmetric errors? Or symmetrize?)  
`return` AnalysisObject  

Combine the specified variations according to the Hession prescription
Central value given by nom.

## combineVariationsLHAPDF(nom, variations, pset, asym=False)
`nom` AnalysisObject (of the nominal variation)  
`variations` dict(String, AnalysisObject)  
`pset` LHAPDF PDFSet (obtained from lhapdf.getPDFSet(pdfSetName))  
`asym` Bool [optional] (return asymmetric errors? Or symmetrize?)  
`return` AnalysisObject  

Combines PDF variations according to the LHAPDF PDFSet prescription.

## combineVariationsReplicas(nom, variations, asym=True)
`nom` AnalysisObject (of the nominal variation)  
`variations` dict(String, AnalysisObject)  
`asym` Bool [optional] (return asymmetric errors? Or symmetrize?)  
`return` AnalysisObject  

Takes the standard deviation of the variations... in principle
This is a placeholder for now... needs to be implemented TODO

## combineVariationsStat(nom)
`nom` AnalysisObject (of the nominal variation)  
`return` AnalysisObject  

Dummy function which currently just grabs the dn/up errs as the stat errs.

## extractTarballsFromDirectory(fulldir, force=False, verbose=False, rootOrYoda='yoda')
`fulldir` String (directory containing tarballs to unpack)  
`force` Bool [optional] (re-do untarring even if already done)  
`verbose` Bool [optional] (Print more output if set to True)  
`rootOrYoda` String [optional] (specify file type of tarred objects)   

`return` list of [sample,newfn,filepath]  

This function goes through a directpry and unpacks all the tarballs it finds.
By default if a matching unpacked file has been found,
a tarball is not needlessly re-unpacked unless option `force` is used.

## getAverageUncertaintySizePerBin(fIn, regexFilter=None, regexVeto=None)
`fIn` String (input file to evaluate)  
`regexFilter` String [optional] (AOs whose names match regex are processed)  
`regexVeto` String [optional] (AOs whose names match regex are NOT processed)  
`return` Float  

Get a rough estimate of the average symmetric relative error per bin, used to
determine what order to plot the uncertaities in.

## getCombinationRecipe(systWeights, combinationRecipeFile=None, combinationRecipeName=None)
`systWeights` list(String) (list of weight types which are avaliable
for a given dsid, to identify the correct combination recipe)  
`combinationRecipeName` String [optional] (specify if you want to use a 
specific comboination uncertainty, otherwise, it will try to auto-
determine frim the systWeights)
`combinationRecipeFile` String [optional] (specify if the target file is in
a different location to $SYSTTOOLSPATH/data/Syst_Database.yaml
`return` String,dict{components,formula}  

Each sample type should have a unique set of weight types/names which allow us 
to identify in the Syst_Database.yaml what the correct combination recipe is.

## getCrossSectionCorrection(xsList, systFiles, nomFiles, rivetNormalised=False, applyKFactorCorrection=False, varWeightName='', nominalWeightName='')
`xsList` list[Float] (list of nominal XS values for the subsamples)  
`systFiles` list[String] (list of file names for the variated subsamples 
(same order as xsList)  
`nomFiles` list[String] (list of file names for the nominal subsamples 
`rivetNormalised` Bool [optional] (whether or not to apply the total yield  
correction for variation AOs, since by default rivet normalises by sum-of-weights  
for that variation rather than sum-of-weight of the nominal)
`applyKFactorCorrection` Bool [optional] (If the sample to merge has a k-factor, there  
is an additional correction needed to avoid double counting the normalisation uncertainty)  
`varWeightName` String [optional] (name of the variation to get weight correction for)  
`nominalWeightName` String [optional] (name of the nominal variation to get weight correction with respect to)  
(same order as xsList)  

This function takes a list of nominal and variated YODA/ROOT files, and a list of
nominal cross-sections in the same order, to calculate the on-the-fly cross-
section corrections needed for each systFile. For yoda files produced in the Rivet 2.x  
series, a correction is needed to undo the fact that Rivet normalises systematics  
by the varied sum of weights instead of the nominal sum of weights.  
If using a k-factor, an additional correction is needed to avoid double-counting  
the normalisation uncertainty. For details, see  
https://twiki.cern.ch/twiki/bin/view/AtlasProtected/PmgSystematicUncertaintyRecipes#On_the_fly_systematic_variations

## getFileKeys(d, basepath='/')
`d` TDirectory or  TFile  (name of file or directory to be read)  
`basepath` String [optional] (absolute path leading to `d` in the file)  
`return` list[String,AnalysisObject]  

Recursively gets the list of Histogram/TGraph names inside a ROOT file directory (supports nested directories)

## getFormulaComponents(formula)
`formula` String  
`return` list(String)  

Take this formula written in this string and recusively optain a list of 
basic components which is needs eg:  
Funcion1(foo, Function2( bar, foo), Function3(foo, Function4 (bar, foo)))  
--> [foo, bar, foo, foo, bar, foo]

## getPlotInfo(aoName, pathInRivetEnv)
`aoName` String (name to access plot info for)  
`pathInRivetEnv` String (.plot file where to get the plot info from)  
`return` dict{String,String} (the list of plotting specifications from the 
.plot file, as a dict)  

Rivet uses a separate .plot file to format plots. We want to use the same
information to format our plots, and this function is a helper to retrieve
that info in a more usable format from a given .plot file.

## getSumOfWeights(path, nominalWeight='')
THIS IS A PLACEHOLDER, @TODO need to use the schema to get this properly working.  
user should also have a command line argument to specify name of sum-of-weights histo

## getXS(dsid, campaign=15)
`dsid` Int (dataset ID of the ATLAS dataset you want to get the
cross-section for  
`campaign` String [optional] (the MC campaign number 15 for mc15_13TeV  
datasets, 16 for mc16_13TeV. If campaign > 16 look instead for a local  
file `data/PMGxsecDB_manual.txt` where the user can input the XS info)  
`return` Float  

Fetch the cross-section (XS) of the dataset whose ID you specified.
This function checks for the DSID in the central PMG XS database on cvmfs
but by setting campaign=999 you can force it to check in a local file:
PMGxsecDB_manual.txt where the user can input the XS info manually.

## lookupPDFSetName(lhapdfid)
`lhapdfid` Int  
`return` String  

Takes a LHAPDF ID code and figures out what is the name of the PDF Set
which it belons to.

## main(argv)
This module can also be run as a standalone executable. 
For info about the options try:
systematicsTool.py -h

This is useful for unpacking large GRID job outputs, as it will do all the 
book-keeping, untarring of tarballs, merging across jobs, across subsamples,
across samples weighted by DSID, and then finally processing the systematics
and plotting them. Multithreading available.

## makeDummyHisto(tg, isLog=False, ratioZoom=None, isRatio=False)
`tg` TGraphAsymmErrors (use this to build the dummy TH1D)  
`isLog` Bool [optional] (Max/Min will need to be adjusted differently
in case of a log plot)  
`ratioZoom` [Float,Float] or None (if not None, use this factor as the low/high limits of the
ratio plot)  
`isRatio` Bool [optional] (Specify whether or not this is a ratio plot)  
`return` TH2D (with appropriate max/min on each axis)  

In order to control better what axis ranges to use, construct a dummy TH1D
with the desired max/min on each axis to be able to show the TGraphs nicely

## makeSystematicsPlotsWithRIVET(mergedSystDict)
`mergedSystDict` dict{String,String} (format is as follows {N!label:filename},
where N is the order in which to plot the variations, labels is what to put on
the legend, and filename is the name of the YODA/ROOT file to use as input)  
`return` None   

Make some ugly plots using modified rivet make-plots. I think they are kinda 
ugly and would suggest using makeSystematicsPlotsWithROOT instead.

## makeSystematicsPlotsWithROOT(mergedSystDict, outdir, nominalName='Nominal', ratioZoom=None, regexFilter=None, regexVeto=None, label='', plotInfo=None)
`mergedSystDict` dict{String,String} (format is as follows {N!label:filename},
where N is the order in which to plot the variations, labels is what to put on
the legend, and filename is the name of the YODA/ROOT file to use as input)  
`outdir` String (the output dir for the plots)  
`nominalName` String (which entry of mergedSystDict is the nominal? One of the
entries in mergedSystDict should have this as a label)  
`ratioZoom` [Float,Float] or None (if not None, use this factor as the low/high limits of the
ratio plot)  
`regexFilter` String [optional] (AOs whose names match regex are processed)  
`regexVeto` String [optional](AOs whose names match regex are NOT processed)  
`label` String [optional](additional label to add to the plots name)  
`plotInfo` String [optional](path to plot info file (will try to find it dynamically if not specified)  
`return` None  

This is a generic macro which will make the output plots of the AOs contained
in the files listed in mergedSystDict.

## mergeInChunks(outfn, paths, progressDict=None, nFilesPerChunk=100, force=False, rootOrYoda='yoda')
`outfn` String (The name of the output Yoda file you wish to produce, with full path)  
`paths` list[String] (List of Yoda files to merge, with full paths)  
`progressText` String [optional] (An extra string to print at the end of the progress message)  
`nFilesPerChunk` int [optional] (How many files to do per chunk)  
`force` Bool [optional] by default, if there is already a matching merged file, this function does nothing
but `force` forces the merge again  
`rootOrYoda` String [optional] (specify file type of objects to merge)   
`return` None

This function safely merges multiple yoda files in chunks of nFilesToProcess at a time (if too many, yodamerge can fail!)
Recommended is nFilesPerChunk=100, but definitely less than 300.

## printProgress(progress)
`progress` dict (a dictionary of the job labels:status)  
`return` Void  

This is a helper function print the progress update in multi-
threader processes

## readFromFile(filename, regexFilter=None, regexVeto=None)
`filename` String (path to file which is to be read)  
`regexFilter` String [optional] (AOs whose names match regex are processed)  
`regexVeto` String [optional](AOs whose names match regex are NOT processed)  
`return` dict{String,AnalysisObject}  

decides whether to process a file as ROOT or YODA depending on the file extension

## readFromROOT(filename, regexFilter=None, regexVeto=None)
`filename` String (path to file which is to be read)  
`regexFilter` String [optional] (AOs whose names match regex are processed)  
`regexVeto` String [optional](AOs whose names match regex are NOT processed)  
`return` dict{String,AnalysisObject}  

Open a ROOT file, read the contents and return them as AnalysisObjects.
Control which AOs to select/reject using the optional regexFilter/regexVeto
arguments.
Only supports TH1D and TGraphAsymmErrors types for now. TODO Support other types.

## readFromYODA(filename, regexFilter=None, regexVeto=None)
`filename` String (path to file which is to be read)  
`regexFilter` String [optional] (AOs whose names match regex are processed)  
`regexVeto` String [optional](AOs whose names match regex are NOT processed)  
`return` dict{String,AnalysisObject}  

Open a YODA file, read the contents and return them as AnalysisObjects.
Control which AOs to select/reject using the optional regexFilter/regexVeto
arguments.
Only supports Histo1D and Scatter2D types for now. TODO Support other types.

## renameFilesWithoutPrefix(directory)
`directory` String  
`return` None  

Renames the files in a given directory, such that the longest common prefix
which occurs in all filenames is ommitted. Useful if your GRID jobs were
submitted with the format <prefix>.<weights_names>.yoda and you want to 
lose the <prefix>. bit !

## resolveFormula(nominal, formula, componentsMap, level=0, verbose=0)
`nominal` AnalysisObject (nominal object, which gives us the centeal values.
The typ is a dict of np.arrays() encoding information equivalent
to a Scatter2D, Histo1D, TGraphAsymmErrors or TH1D. See module description 
for more information about this format. TODO: dedicated class for this?)  
`formula` String (The formula with which to combine the components)  
`componentsMap` dict(filenames,AnalysisObjects) (this is a map between the
file which the component corresponds to and the AnalysisObject which it 
corresponds to.  
`level` Int [optional] (keeps track of how deep in the recursiom we have gone)  
`verbose` Int [optional] (0 or 1, whether or not to print a lot of debug messages)  
`return` AnalysisObject  

Resolves the formula iteratively, using the AnalysisObjects listed in the
componentsMap. resolveFormula supports the following functions:  

`Central(arg1)`: Use the central value of arg1, and set errors to 0  

`DownUpNominal(arg0, arg1, arg2)`: Produce a new AO where ydn is taken from 
the central value ('y') of arg0, yup is taken from the central value of arg1,
and y is taken from the central value of arg2.  

`Envelope(arg1,...,argN)`: Take the min/max envelope of the arguments as the 
up/down errors of the resulting AO. The central value is taken from 'nominal'.  

`Inverse(arg1)`: Take the inverse of the AO, so y--> 1/y, and the errors are 
propagated as yerr --> |yerr/y**2|  

`Product(arg1,...,argN)`: Take the product of the arguments of the function. 
The 'y' value is the product of the 'y' values of the arguments, the errors 
are propagated according the relative errors in quadrature.  

`QuadSum(arg1,..,argN)`: Sum in quadrature of the errors of the arguments,
and central value taken from nominal AO.  

`Value(arg1)`: Used to scan an AO and its errors by the float specified in argi1.

## safeDiv(numerator, denominator)
`numerator` np.array()  
`denominator` np.array()  
`return` np.array  

Avoid NaNs when doing division. Replaces NaN with 0 in the array.

## safeFileName(name, removeExtension=True)
`name` String  
`removeExtension` Bool [optional] (remove the file extension, eg .root, or .yoda)  
`return` String  

Converts an input string (eg a weight name) which may contain non-standard
characters into a string which is safe to use as a filename.
In particular, the following substitutions are made:  
'.' --> 'p'  
' ' --> '_'  
':' --> '_'

## safeRootLatex(unsafeLatex)
`unsafeLatex` String (unsafe Latex string to be converted)  
`return` String (safe TLatex string which can be used on ROOT plots)  

TLatex is not quite the same as regular latex, and won't compiled properly
out of the box unless a few changes are made. This function does that
hopefully in the majority of cases! No promises though... *sigh*

## splitOperand(operand, bracket='()')
`operand` String (The operand which you want to decompose into chunks)  
`bracket` String [optional] (The type of bracket you want to avoid splitting, 
eg (), {}, []...)  
`return` list(String)  

Splits the operand of a formula into comma-separated chunks without splitting
operands of nested functions.
eg: Funcion1(foo, Function2( bar, foo), Function3(foo, Function4 (bar, foo)))  
--> [foo,  Function2( bar, foo), Function3(foo, Function4 (bar, foo))]

## updateProgress(progress, key, message)
`progress` dict (a dictionary of the job labels:status)  
`key` String (the job label to update)  
`message` String (the status to update)  
`return` Void  

This is a helper function to update and print a progress update
when multi-threading

## weightCorrection(var, nom, sampleDir='', varWeightName='', nominalWeightName='')
`var` String (name of YODA file for variation to get the correction for)  
`nom` String (name of YODA file for nominal to get correct away from)  
`sampleDir` String [optional] (directory path for a given subsample)  
`varWeightName` String [optional] (name of the variation to get weight correction for)  
`nominalWeightName` String [optional] (name of the nominal variation to get weight correction with respect to)  
`return` Float  

Computes the on-the-fly weight correction of the systematic variation with
respect to the nominal. This is needed for YODA files because by default,
Rivet normalises by the sum of weights for a given instance (ie variation)
rather than by the nominal. This weight is used to fix the issue.

## writeToFile(histDict, fOut)
`histDict` dict{String,AnalysisObject} (AOs to write out)  
`fOut` String (output file to write to)  
`return` None  

Write AOs out to a fiel, auto-determines root or yoda from file extension of fOut.

## writeToROOT(histDict, fOut)
`histDict` dict{String,AnalysisObject} (AOs to write out)  
`fOut` String (output file to write to)  
`return` None  

Write AOs out to a ROOT file.

## writeToYODA(histDict, fOut)
`histDict` dict{String,AnalysisObject} (AOs to write out)  
`fOut` String (output file to write to)  
`return` None  

Write AOs out to a YODA file.


