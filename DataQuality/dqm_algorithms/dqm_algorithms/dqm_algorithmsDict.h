/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef DQM_ALGORITHMS_DQM_ALGORITHMSDICT_H
#define DQM_ALGORITHMS_DQM_ALGORITHMSDICT_H

#include "dqm_algorithms/AddReference.h"
#include "dqm_algorithms/AddReference_All_Bins_Filled.h"
#include "dqm_algorithms/AddReference_BinContentComp.h"
#include "dqm_algorithms/AddReference_Bins_Diff_FromAvg.h"
#include "dqm_algorithms/AddReference_Bins_Equal_Threshold.h"
#include "dqm_algorithms/AddReference_Bins_GreaterThanEqual_Threshold.h"
#include "dqm_algorithms/AddReference_Bins_GreaterThan_Threshold.h"
#include "dqm_algorithms/AddReference_Bins_LessThanEqual_Threshold.h"
#include "dqm_algorithms/AddReference_Bins_LessThan_Threshold.h"
#include "dqm_algorithms/AddReference_Bins_NotEqual_Threshold.h"
#include "dqm_algorithms/All_Bins_Filled.h"
#include "dqm_algorithms/AveragePrint.h"
#include "dqm_algorithms/BasicGraphCheck.h"
#include "dqm_algorithms/BasicHistoCheck.h"
#include "dqm_algorithms/BasicHistoCheckModuleStatus.h"
#include "dqm_algorithms/BasicStatCheck.h"
#include "dqm_algorithms/BinContentComp.h"
#include "dqm_algorithms/BinContentDump.h"
#include "dqm_algorithms/OccupancyHoleFinder.h"
#include "dqm_algorithms/MDT_OccupancyHoleFinder.h"
#include "dqm_algorithms/RPC_OccupancyHoleFinder.h"
#include "dqm_algorithms/BinDump.h"
#include "dqm_algorithms/BinPrint.h"
#include "dqm_algorithms/BinThresh.h"
#include "dqm_algorithms/BinThreshold.h"
#include "dqm_algorithms/BinHeightThreshold.h"
#include "dqm_algorithms/BinHeight_GreaterThan_Threshold.h"
#include "dqm_algorithms/BinHeight_GreaterThanEqual_Threshold.h"
#include "dqm_algorithms/BinHeight_LessThan_Threshold.h"
#include "dqm_algorithms/BinHeight_LessThanEqual_Threshold.h"
#include "dqm_algorithms/BinHeight_redEqual_yellowGreaterThan_Threshold.h"
#include "dqm_algorithms/BinHeight_redEqual_yellowLessThan_Threshold.h"
#include "dqm_algorithms/BinHeight_Equal_Threshold.h"
#include "dqm_algorithms/BinsDiffByStrips.h"
#include "dqm_algorithms/BinsDiffFromStripMedian.h"
#include "dqm_algorithms/BinsDiffFromPreviousLBs.h"
#include "dqm_algorithms/BinsDiffFromStripMedianOnline.h"
#include "dqm_algorithms/BinsFilledOutRange.h"
#include "dqm_algorithms/BinsOutOfRange.h"
#include "dqm_algorithms/BinsSymmetric.h"
#include "dqm_algorithms/Bins_Diff_FromAvg.h"
#include "dqm_algorithms/Bins_Equal_Threshold.h"
#include "dqm_algorithms/Bins_GreaterThanAbs_Threshold.h"
#include "dqm_algorithms/Bins_GreaterThanEqual_Threshold.h"
#include "dqm_algorithms/Bins_GreaterThanNonZeroMedian_Threshold.h"
#include "dqm_algorithms/Bins_GreaterThan_Threshold.h"
#include "dqm_algorithms/Bins_LessThanAbs_Threshold.h"
#include "dqm_algorithms/Bins_LessThanEqual_Threshold.h"
#include "dqm_algorithms/Bins_LessThanNonZeroMedian_Threshold.h"
#include "dqm_algorithms/Bins_LessThan_Threshold.h"
#include "dqm_algorithms/Bins_NotEqual_Threshold.h"
#include "dqm_algorithms/BlackBin.h"
#include "dqm_algorithms/BlackBin1D.h"
#include "dqm_algorithms/CheckHisto_Mean.h"
#include "dqm_algorithms/CheckHisto_RMS.h"
#include "dqm_algorithms/CheckMean.h"
#include "dqm_algorithms/Chi2Test.h"
#include "dqm_algorithms/Chi2Test_2D.h"
#include "dqm_algorithms/Chi2Test_Chi2.h"
#include "dqm_algorithms/Chi2Test_Chi2_per_NDF.h"
#include "dqm_algorithms/Chi2Test_Prob.h"
#include "dqm_algorithms/Chi2Test_ProbUW.h"
#include "dqm_algorithms/Chi2Test_ProbWW.h"
#include "dqm_algorithms/Chi2Test_Scatterplot.h"
#include "dqm_algorithms/ChiComp.h"
#include "dqm_algorithms/CorrelationYX.h"
#include "dqm_algorithms/CountsBinsGreaterThan.h"
#include "dqm_algorithms/CSCNoisyDead.h"
#include "dqm_algorithms/DivideBin.h"
#include "dqm_algorithms/DivideReference.h"
#include "dqm_algorithms/DivideReference_All_Bins_Filled.h"
#include "dqm_algorithms/DivideReference_BinContentComp.h"
#include "dqm_algorithms/DivideReference_Bins_Diff_FromAvg.h"
#include "dqm_algorithms/DivideReference_Bins_Equal_Threshold.h"
#include "dqm_algorithms/DivideReference_Bins_GreaterThanEqual_Threshold.h"
#include "dqm_algorithms/DivideReference_Bins_GreaterThan_Threshold.h"
#include "dqm_algorithms/DivideReference_Bins_LessThanEqual_Threshold.h"
#include "dqm_algorithms/DivideReference_Bins_LessThan_Threshold.h"
#include "dqm_algorithms/DivideReference_Bins_NotEqual_Threshold.h"
#include "dqm_algorithms/GatherData.h"
#include "dqm_algorithms/GraphPrint.h"
#include "dqm_algorithms/GraphTest.h"
#include "dqm_algorithms/GrubbsOutlierTest.h"
#include "dqm_algorithms/HLTMETComponents.h"
#include "dqm_algorithms/HLTMETStatus.h"
#include "dqm_algorithms/Histogram_Effective_Empty.h"
#include "dqm_algorithms/Histogram_Empty.h"
#include "dqm_algorithms/Histogram_Not_Empty.h"
#include "dqm_algorithms/IterativeGaussianFit.h"
#include "dqm_algorithms/JarqueBeraTest.h"
#include "dqm_algorithms/JarqueBeraTest_JB.h"
#include "dqm_algorithms/JarqueBeraTest_Prob.h"
#include "dqm_algorithms/KillBinsByStrip.h"
#include "dqm_algorithms/KolmogorovTest.h"
#include "dqm_algorithms/KolmogorovTest_MaxDist.h"
#include "dqm_algorithms/KolmogorovTest_MaxDistPlusNorm.h"
#include "dqm_algorithms/KolmogorovTest_Norm.h"
#include "dqm_algorithms/KolmogorovTest_Prob.h"
#include "dqm_algorithms/KurtosisTest.h"
#include "dqm_algorithms/KurtosisTest_GreaterThan.h"
#include "dqm_algorithms/KurtosisTest_GreaterThanAbs.h"
#include "dqm_algorithms/KurtosisTest_LessThan.h"
#include "dqm_algorithms/KurtosisTest_LessThanAbs.h"
#include "dqm_algorithms/LastBinThreshold.h"
#include "dqm_algorithms/L1Calo_OutlierAndFlatnessTest.h"
#include "dqm_algorithms/MDTADCSpectrum.h"
#include "dqm_algorithms/MDTChi2.h"
#include "dqm_algorithms/MDTCluster.h"
#include "dqm_algorithms/MDTMLOverview.h"
#include "dqm_algorithms/MDTMultiplicity.h"
#include "dqm_algorithms/MDTOverview.h"
#include "dqm_algorithms/MDTOverview_Global.h"
#include "dqm_algorithms/MDTOverview_Station.h"
#include "dqm_algorithms/MDTPercentUnderThresh.h"
#include "dqm_algorithms/MDTTDCOfflineSpectrum.h"
#include "dqm_algorithms/MDTTDCSpectrum.h"
#include "dqm_algorithms/MDTTubeCheck.h"
#include "dqm_algorithms/MDTTubeCheckError.h"
#include "dqm_algorithms/MaskedBinRow.h"
#include "dqm_algorithms/MaximumBin.h"
#include "dqm_algorithms/ModuleStatus_All_Bins_Filled.h"
#include "dqm_algorithms/No_OverFlows.h"
#include "dqm_algorithms/No_UnderFlows.h"
#include "dqm_algorithms/OutlierAndFlatnessTest.h"
#include "dqm_algorithms/PassInput.h"
#include "dqm_algorithms/ReferenceMasking.h"
#include "dqm_algorithms/ReferenceMasking_Bins_Diff_FromAvg.h"
#include "dqm_algorithms/ReferenceMasking_Bins_GreaterThan_Threshold.h"
#include "dqm_algorithms/RepeatAlgorithm.h"
#include "dqm_algorithms/RootFit.h"
#include "dqm_algorithms/RootFitGraph.h"
#include "dqm_algorithms/SCTTrackTiming.h"
#include "dqm_algorithms/SideBand.h"
#include "dqm_algorithms/SideBand_Absolute.h"
#include "dqm_algorithms/SideBand_Relative.h"
#include "dqm_algorithms/Simple_doublegaus_Fit.h"
#include "dqm_algorithms/Simple_erf_Fit_Graph.h"
#include "dqm_algorithms/Simple_fermi_Fit.h"
#include "dqm_algorithms/Simple_fermi_Fit_Graph.h"
#include "dqm_algorithms/Simple_gaus_Fit.h"
#include "dqm_algorithms/Simple_gausplusexpo_Fit.h"
#include "dqm_algorithms/Simple_gauspluspol1_Fit.h"
#include "dqm_algorithms/Simple_landau_Fit.h"
#include "dqm_algorithms/Simple_pol1_Fit.h"
#include "dqm_algorithms/Simple_sinusoid_Fit.h"
#include "dqm_algorithms/SkewnessTest.h"
#include "dqm_algorithms/SkewnessTest_GreaterThan.h"
#include "dqm_algorithms/SkewnessTest_GreaterThanAbs.h"
#include "dqm_algorithms/SkewnessTest_LessThan.h"
#include "dqm_algorithms/SkewnessTest_LessThanAbs.h"
#include "dqm_algorithms/TRTCheckPeakSimple.h"
#include "dqm_algorithms/TRTHistogramHasNonZeroEntries.h"
#include "dqm_algorithms/TripleGaussCollFit.h"
#include "dqm_algorithms/LastBinThresholdAction.h"
#endif // DQM_ALGORITHMS_DQM_ALGORITHMSDICT_H
