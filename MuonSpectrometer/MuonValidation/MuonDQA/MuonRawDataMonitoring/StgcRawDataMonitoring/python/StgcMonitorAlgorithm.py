#
#Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
#

from AthenaConfiguration.ComponentFactory import CompFactory
#from StgcRawDataMonitoring.StgcMonUtils import *

def StgcMonitoringConfig(inputFlags):
    '''Function to configures some algorithms in the monitoring system.'''
    ### STEP 1 ###
    # Define one top-level monitoring algorithm. The new configuration 
    # framework uses a component accumulator.
    from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
    result = ComponentAccumulator()
    # Make sure muon geometry is configured
    from MuonConfig.MuonGeometryConfig import MuonGeoModelCfg
    result.merge(MuonGeoModelCfg(inputFlags))

    # The following class will make a sequence, configure algorithms, and link
    # them to GenericMonitoringTools

    from AthenaCommon.AppMgr import ServiceMgr
    ServiceMgr.Dump = False

    from AthenaMonitoring import AthMonitorCfgHelper
    helper = AthMonitorCfgHelper(inputFlags,'StgcAthMonitorCfg')

    # Adding an algorithm to the helper.

    sTGCMonAlg = helper.addAlgorithm(CompFactory.StgcRawDataMonAlg,'sTGCMonAlg')
    sTGCMonAlg.DoSTGCESD = True  

    # Add a generic monitoring tool (a "group" in old language). The returned      
    # object here is the standard GenericMonitoringTool. 

    sTGCGroup = helper.addGroup(sTGCMonAlg,'sTGCMonitor','Muon/MuonRawDataMonitoring/sTGC/')

    # Configure histograms
    # Overview histograms
    #sTGCGroup.defineHistogram('residual;Residuals',  type='TH1F',  title='Residuals;res[mm];Number of Entries', path='Overview',   xbins=200, xmin=-10, xmax=10.) #yes

    #sTGCGroup.defineHistogram('residual,eta_trk;Res_vs_eta', type='TH2F', title="Res vs Eta;Res;Eta;", path='Overview',xbins=100, xmin=-10, xmax=10., ybins=100, ymin=1.24, ymax=1.4) #yes

    #sTGCGroup.defineHistogram('residual,phi_trk;Res_vs_phi', type='TH2F', title="Res vs Phi;Res;Phi;", path='Overview',xbins=100, xmin=-10, xmax=10., ybins=16, ymin=-3.14, ymax=3.14) #yes
    
    #sTGCGroup.defineHistogram('residual,stPhi_mon;Res_vs_stPhi', type='TH2F', title="Res vs stPhi;Res;stPhi;", path='Overview',xbins=100, xmin=-10, xmax=10., ybins=16, ymin=0, ymax=16) #yes
    
    
    sTGCGroup.defineHistogram('strip_times;Strip_Time', type = 'TH1F', title = 'Strip Time; Strip Time [ns];Number of Entries', path = 'Overview', xbins = 20, xmin = 0., xmax = 100.) #new
    sTGCGroup.defineHistogram('strip_charges;Strip_Charge', type = 'TH1F', title = 'Strip Charge; Strip Charge [fC];Number of Entries', path = 'Overview', xbins = 200, xmin = 0., xmax = 1000.) #new
    sTGCGroup.defineHistogram('strip_number;Strip_Number', type = 'TH1F', title = 'Strip Number; Strip number;Number of Entries', path = 'Overview', xbins = 20, xmin = 0., xmax = 400.) #new
    sTGCGroup.defineHistogram('charge_all;Charge', type = 'TH1F', title = 'Charge;Charge[fC];Number of Entries', path = 'Overview', xbins = 100, xmin = 0., xmax = 1000.) #yes
    
    sTGCGroup.defineHistogram('x_mon,y_mon;Posx_vs_Posy', type = 'TH2F', title="Posx vs Posy;sTGC-GlobalX [mm];sTGC-GlobalY [mm];", path = 'Overview', xbins = 500, xmin = -5000, xmax = 5000., ybins = 500, ymin = -5000., ymax = 5000.) #yes
    
    sTGCGroup.defineHistogram('R_mon,z_mon;R_vs_Posz', type = 'TH2F', title = "R vs Posz; sTGC-GlobalR [mm]; sTGC-GlobalZ [mm];", path = 'Overview', xbins = 500, xmin = 2500., xmax = 5000., ybins = 1000, ymin = 6800 ,ymax = 8000) #yes
    
    sTGCGroup.defineHistogram('numberofstrips_percluster;Number_of_strips_percluster', type = 'TH1F', title = 'Number of strips per cluster;Number of strips;Number of Entries', path = 'Overview', xbins = 12, xmin = 0., xmax = 12.) #yes
#    mmGroup.defineHistogram('mu_TPC_angle;uTPC_angle',  type='TH1F',
#                            title='#mu TPC angle;#mu TPC angle [degrees];Number of Entries',
#                            path='Overview',   xbins=2000, xmin=-100, xmax=100)

#    mmGroup.defineHistogram('mu_TPC_chi2;uTPC_chi2',  type='TH1F',
#                        title='#mu TPC #chi2; #mu TPC #chi2;Number of Entries',
#                        path='Overview',   xbins=100, xmin=0., xmax=1.)

    sTGCGroup.defineHistogram('time_all;Time', type = 'TH1F', title = 'Time;Time[ns];Number of Entries', path = 'Overview', xbins = 5, xmin = 0., xmax = 5.) #yes

    side = ["CSide", "ASide"]

    for iside in side:
        sTGC_SideGroup = "sTGC_sideGroup{0}".format(iside)
        stgcSideGroup    = helper.addGroup(sTGCMonAlg, sTGC_SideGroup, "Muon/MuonRawDataMonitoring/sTGC/" + iside)
        for multip in range(1, 3):
            for gasgap in range(1, 5):
                title_chargePad_phi_vs_eta = f'Charge (pad): {iside} Multiplet {multip} Gas gap {gasgap}; stationPhi; stationEta; Total charge [fC]'
                var_chargePad_phi_vs_eta = f'sector_{iside}_phi_multiplet_{multip}_gasgap_{gasgap}, sector_{iside}_eta_multiplet_{multip}_gasgap_{gasgap};ChargePad_vs_phi_vs_eta_{iside}_multiplet_{multip}_gasgap_{gasgap}'
                if (f'{iside}' == 'ASide'):
                    stgcSideGroup.defineHistogram(var_chargePad_phi_vs_eta, type = 'TH2F', title = title_chargePad_phi_vs_eta, path = 'Summary', xbins = 8, xmin = 1., xmax = 9., ybins = 3, ymin = 1., ymax = 4., opt = 'kAlwaysCreate', weight = f'charge_pad_{iside}_multiplet_{multip}_gasgap_{gasgap}')
                else:
                    stgcSideGroup.defineHistogram(var_chargePad_phi_vs_eta, type = 'TH2F', title = title_chargePad_phi_vs_eta, path = 'Summary', xbins = 8, xmin = 1., xmax = 9., ybins = 3, ymin = -4., ymax = -1., opt = 'kAlwaysCreate', weight = f'charge_pad_{iside}_multiplet_{multip}_gasgap_{gasgap}')
                
                title_chargeStrip_phi_vs_eta = f'Charge (strip): {iside} Multiplet {multip} Gas gap {gasgap}; stationPhi; stationEta; Total charge [fC]'
                var_chargeStrip_phi_vs_eta = f'sector_{iside}_phi_multiplet_{multip}_gasgap_{gasgap}, sector_{iside}_eta_multiplet_{multip}_gasgap_{gasgap};ChargeStrip_vs_phi_vs_eta_{iside}_multiplet_{multip}_gasgap_{gasgap}'
               
                if (f'{iside}' == 'ASide'):
                    stgcSideGroup.defineHistogram(var_chargeStrip_phi_vs_eta, type = 'TH2F', title = title_chargeStrip_phi_vs_eta, path = 'Summary', xbins = 8, xmin = 1., xmax = 9., ybins = 3, ymin = 1., ymax = 4., opt = 'kAlwaysCreate', weight = f'charge_strip_{iside}_multiplet_{multip}_gasgap_{gasgap}')
                else:
                    stgcSideGroup.defineHistogram(var_chargeStrip_phi_vs_eta, type = 'TH2F', title = title_chargeStrip_phi_vs_eta, path = 'Summary', xbins = 8, xmin = 1., xmax = 9., ybins = 3, ymin = -4., ymax = -1., opt = 'kAlwaysCreate', weight = f'charge_strip_{iside}_multiplet_{multip}_gasgap_{gasgap}')

                title_chargeWire_phi_vs_eta = f'Charge (wire): {iside} Multiplet {multip} Gas gap {gasgap}; stationPhi; stationEta; Total charge [fC]'
                var_chargeWire_phi_vs_eta = f'sector_{iside}_phi_multiplet_{multip}_gasgap_{gasgap}, sector_{iside}_eta_multiplet_{multip}_gasgap_{gasgap};ChargeWire_vs_phi_vs_eta_{iside}_multiplet_{multip}_gasgap_{gasgap}'
 
                if (f'{iside}' == 'ASide'):
                    stgcSideGroup.defineHistogram(var_chargeWire_phi_vs_eta, type = 'TH2F', title = title_chargeWire_phi_vs_eta, path = 'Summary', xbins = 8, xmin = 1., xmax = 9., ybins = 3, ymin = 1., ymax = 4., opt = 'kAlwaysCreate', weight = f'charge_wire_{iside}_multiplet_{multip}_gasgap_{gasgap}')
                else:
                    stgcSideGroup.defineHistogram(var_chargeWire_phi_vs_eta, type = 'TH2F', title = title_chargeWire_phi_vs_eta, path = 'Summary', xbins = 8, xmin = 1., xmax = 9., ybins = 3, ymin = -4., ymax = -1., opt = 'kAlwaysCreate', weight = f'charge_wire_{iside}_multiplet_{multip}_gasgap_{gasgap}')

                for phiStation in range(1, 9):
                    title_stationEta_vs_stripNumber_vs_chargeStrip_eachPhi = f'Station eta vs strip number vs charge (strip): {iside} Multiplet {multip} Gas gap {gasgap} stationPhi {phiStation}; stationEta; Strip number; Total charge [fC]'
                    var_stationEta_vs_stripNumber_vs_chargeStrip_eachPhi = f'sector_{iside}_eta_multiplet_{multip}_gasgap_{gasgap}_stationPhi_{phiStation}, stripNumber_strip_{iside}_multiplet_{multip}_gasgap_{gasgap}_stationPhi_{phiStation};StationEta_vs_stripNumber_vs_chargePad_{iside}_multiplet_{multip}_gasgap_{gasgap}_stationPhi_{phiStation}'
                    
                    if (f'{iside}' == 'ASide'):
                        stgcSideGroup.defineHistogram(var_stationEta_vs_stripNumber_vs_chargeStrip_eachPhi, type = 'TH2F', title = title_stationEta_vs_stripNumber_vs_chargeStrip_eachPhi, path = 'Summary', xbins = 3, xmin = 1., xmax = 4., ybins = 400, ymin = 0., ymax = 400., opt = 'kAlwaysCreate', weight = f'charge_strip_{iside}_multiplet_{multip}_gasgap_{gasgap}_stationPhi_{phiStation}')
                    else:
                        stgcSideGroup.defineHistogram(var_stationEta_vs_stripNumber_vs_chargeStrip_eachPhi, type = 'TH2F', title = title_stationEta_vs_stripNumber_vs_chargeStrip_eachPhi, path = 'Summary', xbins = 3, xmin = -4., xmax = -1., ybins = 400, ymin = 0., ymax = 400., opt = 'kAlwaysCreate', weight = f'charge_strip_{iside}_multiplet_{multip}_gasgap_{gasgap}_stationPhi_{phiStation}')
    """
    
    side = ["CSide","ASide"]
    
    for iside in side:
        if iside == "ASide":
            thisLabelx11, thisLabely11 = getsTGCLabel("x_lab_occ_ASide", "y_lab_occ_ASide")
        if iside == "CSide":
            thisLabelx11, thisLabely11 = getsTGCLabel("x_lab_occ_CSide", "y_lab_occ_CSide")

        sTGC_SideGroup = "sTGC_sideGroup{0}".format(iside)
        stgcSideGroup    = helper.addGroup(sTGCMonAlg, sTGC_SideGroup, "Muon/MuonRawDataMonitoring/sTGC/" + iside)
        
        phimax = 16

        for phi in range(1, phimax + 1):
            thisLabely=getsTGCLabelY("y_lab_occ_lb")
            stgcSideGroup.defineHistogram(f'lb_mon,sector_lb_{iside}_phi{phi};Occupancy_lb_{iside}_phi{phi}', type='TH2F', title=f'Occupancy wrt lb sector {phi}; LB; PCB', path='Occupancy',  xbins=20, xmin=-0.5, xmax=99.5, opt='kAddBinsDynamically,kAlwaysCreate', ybins=20, ymin=0., ymax=64, ylabels=thisLabely)
            
        for gas1 in range(1, 5):
            for multi1 in range(1, 3):
                title_ontrack=f'Posy vs Posx {iside} multiplet{multi1} gap{gas1} ontrack; sTGC-GlobalX [mm]; sTGC-GlobalY [mm]'
                var_ontrack=f'x_{iside}_multiplet{multi1}_gas_gap_{gas1}_ontrack,y_{iside}_multiplet{multi1}_gas_gap_{gas1}_ontrack;Posy_vs_Posx_{iside}_multiplet{multi1}_gas_gap_{gas1}_ontrack'
                stgcSideGroup.defineHistogram(var_ontrack, type='TH2F', title=title_ontrack, path='PosY_vs_Posx_perLayer_ontrack',xbins=500, xmin=-5000, xmax=5000., ybins=500, ymin=-5000.,ymax=5000., opt='kAlwaysCreate')
    """
    
    #sTGCGroup.defineHistogram('statEta_strip,strip_number;Strip_Numbers_vs_StationEta', type='TH2F',title='Strip Numbers vs Station Eta;; Strip Numbers;',path='Overview', xbins=5, xmin=-2, xmax=3., xlabels=['#eta-2','#eta-1','','#eta1','#eta2'], ybins=5120, ymin=0., ymax=5120.)

#    thisLabelx,thisLabely=getsTGCLabel("x_lab_occ_etaminus1","y_lab_occ_etaminus1")

#    mmGroup.defineHistogram('sector_CSide_eta1,stationPhi_CSide_eta1;Occupancy_CSide_eta1_PCB', type='TH2F', title='Occupancy CSide eta1 PCB; ; ;', path='Occupancy', xbins=40, xmin=0, xmax=40., ybins=16, ymin=1, ymax=17,xlabels=thisLabelx, ylabels=thisLabely)

#    mmGroup.defineHistogram('sector_CSide_eta1_ontrack,stationPhi_CSide_eta1_ontrack;Occupancy_CSide_eta1_PCB_ontrack', type='TH2F', title='Occupancy CSide eta1 PCB ontrack; ; ;', path='Occupancy_ontrack', xbins=40, xmin=0, xmax=40., ybins=16, ymin=1, ymax=17,xlabels=thisLabelx, ylabels=thisLabely)

#    thisLabelx1,thisLabely1=getsTGCLabel("x_lab_occ_etaminus2","y_lab_occ_etaminus2")

#    mmGroup.defineHistogram('sector_CSide_eta2,stationPhi_CSide_eta2;Occupancy_CSide_eta2_PCB', type='TH2F', title='Occupancy CSide eta2 PCB; ; ;', path='Occupancy', xbins=24, xmin=0, xmax=24., ybins=16, ymin=1, ymax=17,xlabels=thisLabelx1, ylabels=thisLabely1)

#    mmGroup.defineHistogram('sector_CSide_eta2_ontrack,stationPhi_CSide_eta2_ontrack;Occupancy_CSide_eta2_PCB_ontrack', type='TH2F', title='Occupancy CSide eta2 PCB ontrack; ; ;', path='Occupancy_ontrack', xbins=24, xmin=0, xmax=24., ybins=16, ymin=1, ymax=17,xlabels=thisLabelx1, ylabels=thisLabely1)

#    thisLabelx2,thisLabely2=getsTGCLabel("x_lab_occ_eta1","y_lab_occ_eta1")

#    mmGroup.defineHistogram('sector_ASide_eta1,stationPhi_ASide_eta1;Occupancy_ASide_eta1_PCB', type='TH2F', title='Occupancy ASide eta1 PCB; ; ;', path='Occupancy', xbins=40, xmin=0, xmax=40., ybins=16, ymin=1, ymax=17,xlabels=thisLabelx2, ylabels=thisLabely2)

#    mmGroup.defineHistogram('sector_ASide_eta1_ontrack,stationPhi_ASide_eta1_ontrack;Occupancy_ASide_eta1_PCB_ontrack', type='TH2F', title='Occupancy ASide eta1 PCB ontrack; ; ;', path='Occupancy_ontrack', xbins=40, xmin=0, xmax=40., ybins=16, ymin=1, ymax=17,xlabels=thisLabelx2, ylabels=thisLabely2)

#    thisLabelx3,thisLabely3=getsTGCLabel("x_lab_occ_eta2","y_lab_occ_eta2")

#    mmGroup.defineHistogram('sector_ASide_eta2,stationPhi_ASide_eta2;Occupancy_ASide_eta2_PCB', type='TH2F', title='Occupancy ASide eta2 PCB; ; ;', path='Occupancy', xbins=24, xmin=0, xmax=24., ybins=16, ymin=1, ymax=17,xlabels=thisLabelx3, ylabels=thisLabely3)

#    mmGroup.defineHistogram('sector_ASide_eta2_ontrack,stationPhi_ASide_eta2_ontrack;Occupancy_ASide_eta2_PCB_ontrack', type='TH2F', title='Occupancy ASide eta2 PCB ontrack; ; ;', path='Occupancy_ontrack', xbins=24, xmin=0, xmax=24., ybins=16, ymin=1, ymax=17,xlabels=thisLabelx3, ylabels=thisLabely3)

#    thisLabely4=getsTGCLabelY("y_lab_lb_CSide_eta2")

#    mmGroup.defineHistogram('lb_mon,sector_lb_CSide_eta2;Occupancy_lb_CSide_eta2_PCB', type='TH2F', title="Occupancy lb CSide eta2 PCB; ; ;", path='Occupancy', xbins=100,xmin=-0.5,xmax=99.5,ymin=1, ymax=49,ylabels=thisLabely4,opt='kAddBinsDynamically')

#    mmGroup.defineHistogram('lb_ontrack,sector_lb_CSide_eta2_ontrack;Occupancy_lb_CSide_eta2_PCB_ontrack', type='TH2F', title="Occupancy lb CSide eta2 PCB ontrack; ; ;", path='Occupancy_ontrack', xbins=100,xmin=-0.5,xmax=99.5,ymin=1, ymax=49, ylabels=thisLabely4, opt='kAddBinsDynamically')

#    thisLabely5=getsTGCLabelY("y_lab_lb_CSide_eta1")

#    mmGroup.defineHistogram('lb_mon,sector_lb_CSide_eta1;Occupancy_lb_CSide_eta1_PCB', type='TH2F', title="Occupancy lb CSide eta1 PCB; ; ;", path='Occupancy', xbins=100,xmin=-0.5,xmax=99.5,ymin=1, ymax=81,ylabels=thisLabely5,opt='kAddBinsDynamically')

#    mmGroup.defineHistogram('lb_ontrack,sector_lb_CSide_eta1_ontrack;Occupancy_lb_CSide_eta1_PCB_ontrack', type='TH2F', title="Occupancy lb CSide eta1 PCB ontrack; ; ;", path='Occupancy_ontrack', xbins=100,xmin=-0.5,xmax=99.5,ymin=1, ymax=81,ylabels=thisLabely5,opt='kAddBinsDynamically')

#    thisLabely6=getsTGCLabelY("y_lab_lb_ASide_eta1")

#    mmGroup.defineHistogram('lb_mon,sector_lb_ASide_eta1;Occupancy_lb_ASide_eta1_PCB', type='TH2F', title="Occupancy lb ASide eta1 PCB; ; ;", path='Occupancy', xbins=100,xmin=-0.5,xmax=99.5,ymin=1, ymax=81,ylabels=thisLabely6,opt='kAddBinsDynamically')

#    mmGroup.defineHistogram('lb_ontrack,sector_lb_ASide_eta1_ontrack;Occupancy_lb_ASide_eta1_PCB_ontrack', type='TH2F', title="Occupancy lb ASide eta1 PCB ontrack; ; ;", path='Occupancy_ontrack', xbins=100,xmin=-0.5,xmax=99.5,ymin=1, ymax=81,ylabels=thisLabely6,opt='kAddBinsDynamically')

#    thisLabely7=getsTGCLabelY("y_lab_lb_ASide_eta2")

#    mmGroup.defineHistogram('lb_mon,sector_lb_ASide_eta2;Occupancy_lb_ASide_eta2_PCB', type='TH2F', title="Occupancy lb ASide eta2 PCB; ; ;", path='Occupancy', xbins=100,xmin=-0.5,xmax=99.5,ymin=1, ymax=49,ylabels=thisLabely7,opt='kAddBinsDynamically')

#    mmGroup.defineHistogram('lb_ontrack,sector_lb_ASide_eta2_ontrack;Occupancy_lb_ASide_eta2_PCB_ontrack', type='TH2F', title="Occupancy lb ASide eta2 PCB ontrack; ; ;", path='Occupancy_ontrack', xbins=100,xmin=-0.5,xmax=99.5,ymin=1, ymax=49,ylabels=thisLabely7,opt='kAddBinsDynamically')

#    side = ["CSide","ASide"]
#    sector = ["S","L"]
#    etasector  = ["1","2"]
#    for iside in side:
#        if iside=="ASide":
#            thisLabelx11,thisLabely11=getsTGCLabel("x_lab_occ_ASide","y_lab_occ_ASide")
#        if iside=="CSide":
#            thisLabelx11,thisLabely11=getsTGCLabel("x_lab_occ_CSide","y_lab_occ_CSide")
#        sTGC_SideGroup="sTGC_sideGroup{0}".format(iside)
#        sTGCSideGroup=helper.addGroup(sTGCMonAlg, sTGC_SideGroup, "Muon/MuonRawDataMonitoring/sTGC/"+iside)
        # Histograms for each sector
#        phimax=8
#        multipletmin=1
#        multipletmax=2
#        for isector in sector:
#            for phi in range(1, phimax+1):
#                title_sTGCSummary="Number of strips per cluster,"+iside+" "+isector+" stPhi "+str(phi)   
#                var="sector_strip_"+iside+"_"+isector+"_phi"+str(phi)+",strip_number_"+iside+"_"+isector+"_phi"+str(phi)+";Strip_number_pergap_"+iside+"_"+isector+"stPhi"+str(phi)
#                mmSideGroup.defineHistogram(var, type='TH2F', title=title_sTGCSummary+"; ;Strip Number",      path='Number_of_strips_percluster_perPhiSector',   xbins=16, xmin=0, xmax=16, xlabels=thisLabelx11, ybins=5120, ymin=0., ymax=5120.)
#                for eta in etasector:
#                    for multi in range(multipletmin, multipletmax+1):
#                        for gas_gap in range(1,5):
                            # Histograms for each layer
#                            title_sTGCSummary_charge="Charge "+iside+" "+isector+" stPhi"+str(phi)+" stEta"+str(eta)+" multiplet"+str(multi)+" gap"+str(gas_gap)
#                            var1="charge_"+iside+"_sector_"+isector+"_phi"+str(phi)+"_stationEta"+str(eta)+"_multiplet"+str(multi)+"_gas_gap"+str(gas_gap)+";Charge_"+iside+"_"+isector+"_stPhi"+str(phi)+"_stEta"+str(eta)+"_multiplet"+str(multi)+"_gap"+str(gas_gap)
#                            mmSideGroup.defineHistogram(var1,  type='TH1F', title=title_sTGCSummary_charge+';Charge [fC];Number of Entries',path='Charge_perLayer',   xbins=120, xmin=0., xmax=1200.)
#                            title_sTGCSummary_angle="uTPC angle "+iside+" "+isector+" stPhi"+str(phi)+" stEta"+str(eta)+" multiplet"+str(multi)+" gap"+str(gas_gap)
#                            var3="mu_TPC_angle_"+iside+"_sector_"+isector+"_phi"+str(phi)+"_stationEta"+str(eta)+"_multiplet"+str(multi)+"_gas_gap"+str(gas_gap)+";uTPCangle_"+iside+"_"+isector+"_stPhi"+str(phi)+"_stEta"+str(eta)+"_multiplet"+str(multi)+"_gap"+str(gas_gap)
#                            mmSideGroup.defineHistogram(var3,  type='TH1F', title=title_sTGCSummary_angle+"; #muTPC angle [degrees];Number of Entries",path='uTPC_angle_perLayer',    xbins=2000, xmin=-100, xmax=100)
#                            var_residual="residuals_"+iside+"_phi"+str(phi)+"_stationEta"+str(eta)+"_multiplet"+str(multi)+"_gas_gap"+str(gas_gap)
#                            print(var_residual)
#                            title_residual = "residuals "+iside+" "+isector+" stPhi"+str(phi)+" stEta"+str(eta)+" multiplet"+str(multi)+" gap"+str(gas_gap)
#                            mmSideGroup.defineHistogram(var_residual,  type='TH1F', title=title_residual+"; res [mm];Number of Entries",path='Residuals',    xbins=200, xmin=-10, xmax=10)
#        for gas1 in range(1, 5):
#            for multi1 in range(1, 3):
#                title_ontrack="Posy vs Posx "+iside+" multiplet"+str(multi1)+" gap"+str(gas1)+" ontrack"
#                var_ontrack="x_"+iside+"_multiplet"+str(multi1)+"_gas_gap_"+str(gas1)+"_ontrack,y_"+iside+"_multiplet"+str(multi1)+"_gas_gap_"+str(gas1)+"_ontrack;Posy_vs_Posx_"+iside+"_multiplet"+str(multi1)+"_gas_gap_"+str(gas1)+"_ontrack"
#                mmSideGroup.defineHistogram(var_ontrack, type='TH2F', title=title_ontrack+";sTGC-GlobalX [mm];sTGC-GlobalY [mm];", path='PosY_vs_Posx_perLayer_ontrack',xbins=500, xmin=-5000, xmax=5000., ybins=500, ymin=-5000.,ymax=5000.)
#    mmMonAlg.TriggerChain = ''
    ####acc, seq = helper.result()
    acc = helper.result()
    result.merge(acc)
    return result
if __name__=='__main__':
    # Setup the Run III behavior
    from AthenaCommon.Configurable import Configurable
    Configurable.configurableRun3Behavior = 1
    #from AthenaCommon.AppMgr import ServiceMgr
    #ServiceMgr.Dump = False
    #from AthenaCommon.Constants import DEBUG
    #from AthenaCommon.Logging import log
    # Set the Athena configuration flags
    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    #ConfigFlags.Input.Files = ['/eos/home-s/sfuenzal/NSWSoftware_0804/MuonSpectrometerPackages_2504/InputData/group.det-muon/group.det-muon/group.det-muon.28531224.EXT0._000050.ESD.pool.root']
    

    """
    ConfigFlags.Input.Files = ['/eos/home-s/sfuenzal/NSWSoftware_0804/MuonSpectrometerPackages_0506/InputData/group.det-muon/group.det-muon.29033949.EXT1._001021.ESD.pool.root',
                               '/eos/home-s/sfuenzal/NSWSoftware_0804/MuonSpectrometerPackages_0506/InputData/group.det-muon/group.det-muon.29033949.EXT1._001455.ESD.pool.root'
    ]
    """
    
    ConfigFlags.Input.Files = ['/eos/home-s/sfuenzal/NSWSoftware_0804/MuonSpectrometerPackages_0806/InputData/mc21/ESD.29004502._000074.pool.root.1'
    ]
    
    #from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
    ConfigFlags.Output.HISTFileName = 'monitor_sTGC.root'

    ConfigFlags.Detector.GeometrysTGC=True
    #ConfigFlags.Muon.doMicromegas=True
    ConfigFlags.DQ.useTrigger=False

    ConfigFlags.lock()
    ConfigFlags.dump()
    # Initialize configuration object, add accumulator, merge, and run.
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg 
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    #from MuonConfig.MuonPrepDataConvConfig import MuonPrepDataConvCfg
    #MuonPrepDataConvCfgAcc = MuonPrepDataConvCfg(ConfigFlags)
    #MuonPrepDataConvCfgAcc.OutputLevel=2
    cfg = MainServicesCfg(ConfigFlags)
    cfg.merge(PoolReadCfg(ConfigFlags))

    #cfg.merge(MuonPrepDataConvCfgAcc)
    sTGCMonitorAcc  =  StgcMonitoringConfig(ConfigFlags)
    #sTGCMonitorAcc.OutputLevel=2
    sTGCMonitorAcc.OutputLevel=2
    #log.setLevel(DEBUG)
    cfg.merge(sTGCMonitorAcc)           
    #cfg.printConfig(withDetails=True, summariseProps = True) 
    # number of events selected in the ESD
    cfg.run(100)


