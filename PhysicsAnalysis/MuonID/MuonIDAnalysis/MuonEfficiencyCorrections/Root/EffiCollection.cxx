/*
 Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 */

#include <MuonEfficiencyCorrections/EffiCollection.h>
#include <MuonEfficiencyCorrections/EfficiencyScaleFactor.h>
#include <MuonEfficiencyCorrections/MuonEfficiencyScaleFactors.h>
#include <MuonEfficiencyCorrections/UtilFunctions.h>
#include <TTree.h>
namespace CP {
    namespace {
        static const std::vector<std::string> ToRemove { "GeV", "MeV", "[", "]", "{", "}", "(", ")", "#", " " };
        typedef std::pair<std::string,std::string> stringpair;
        static const std::vector<stringpair> ToReplace { stringpair("-","minus"), stringpair(".","p")};
        
    }
    EffiCollection::EffiCollection(const MuonEfficiencyScaleFactors& ref_tool) :
            m_ref_tool(ref_tool),
            m_central_eff(),
            m_calo_eff(),
            m_forward_eff(),
            m_lowpt_central_eff(),
            m_lowpt_calo_eff(),
            m_syst_set(std::make_unique<SystematicSet>()){

        m_central_eff = std::make_shared<CollectionContainer>(m_ref_tool, CollectionType::Central);
        if (m_ref_tool.filename_Calo() != m_ref_tool.filename_Central()){
            m_calo_eff = std::make_shared<CollectionContainer>(m_ref_tool,CollectionType::Calo);
        } else m_calo_eff = m_central_eff;
        
        if (m_ref_tool.filename_HighEta() != m_ref_tool.filename_Central()){
            m_forward_eff = std::make_shared<CollectionContainer>(m_ref_tool,CollectionType::Forward);
        } else m_forward_eff = m_central_eff;
        
        if (m_ref_tool.lowPtTransition() > 0  && m_ref_tool.filename_LowPt() != m_ref_tool.filename_Central()){
            m_lowpt_central_eff = std::make_shared<CollectionContainer>(m_ref_tool,CollectionType::CentralLowPt);
        } else m_lowpt_central_eff = m_central_eff;
       
         if (m_ref_tool.lowPtTransition() > 0  && m_ref_tool.filename_LowPtCalo() != m_ref_tool.filename_Central()){
            m_lowpt_calo_eff = std::make_shared<CollectionContainer>(m_ref_tool, CollectionType::CaloLowPt);
        } else m_lowpt_calo_eff = m_central_eff;
    }
    
    EffiCollection::EffiCollection(const EffiCollection* Nominal, const MuonEfficiencyScaleFactors& ref_tool, const std::string& syst, int syst_bit_map, bool is_up):
            m_ref_tool(ref_tool),
            m_central_eff(),
            m_calo_eff(),
            m_forward_eff(),
            m_lowpt_central_eff(),
            m_lowpt_calo_eff(),
            m_syst_set() {
    
        if (is_up) syst_bit_map |= EffiCollection::UpVariation;
        /// Use a lambda function to assign the maps easily
        std::function< std::shared_ptr<CollectionContainer>(CollectionType)>  make_variation = [this,&ref_tool, Nominal, syst_bit_map, syst](CollectionType type){
                if (syst_bit_map & type) {
                    return std::make_shared<CollectionContainer>(ref_tool,
                                                                 Nominal->retrieveContainer(type).get(),
                                                                 syst, syst_bit_map);                                                                             }
                /// Only the Z reconstruction analysis has different files for bulk / calo-tag / low-pt and forward eta
                if (ref_tool.measurement() != MuonEfficiencyType::Reco) return m_central_eff;
                return Nominal->retrieveContainer(type);
            };
        m_central_eff = make_variation(CollectionType::Central);
        m_forward_eff = make_variation(CollectionType::Forward);
        m_calo_eff = make_variation(CollectionType::Calo);
        m_lowpt_central_eff = make_variation(CollectionType::CentralLowPt);
        m_lowpt_calo_eff = make_variation(CollectionType::CaloLowPt);
        
        
    }
    
    std::shared_ptr<CollectionContainer> EffiCollection::retrieveContainer(CollectionType Type) const {
        if (Type == CollectionType::Central) return m_central_eff;
        if (Type == CollectionType::Forward) return m_forward_eff;
        if (Type == CollectionType::Calo) return m_calo_eff;
        if (Type == CollectionType::CentralLowPt) return m_lowpt_central_eff;
        if (Type == CollectionType::CaloLowPt) return m_lowpt_calo_eff;
        return std::shared_ptr<CollectionContainer>();
    }
    bool EffiCollection::CheckConsistency()  {
        if (!m_central_eff || !m_central_eff->CheckConsistency()) {
            Error("EffiCollection()", "Consistency check for central file failed");
            return false;
        }
        if (!m_calo_eff || !m_calo_eff->CheckConsistency()) {
            Error("EffiCollection()", "Consistency check for calo file failed");            
            return false;
        }
        if (!m_forward_eff || !m_forward_eff->CheckConsistency()) {
             Error("EffiCollection()", "Consistency check for forward file failed");            
            return false;
        }
        if (!m_lowpt_central_eff || !m_lowpt_central_eff->CheckConsistency()) {
            Error("EffiCollection()", "Consistency check for low-pt file failed");    
            return false;
        }
        if (!m_lowpt_calo_eff || !m_lowpt_calo_eff->CheckConsistency()) {
            Error("EffiCollection()", "Consistency check for low-pt calo file failed"); 
            return false;
        }
        /// At this stage we know that all efficiencies have been loaded 
        /// successfully. We need to now to order the maps to make global
        /// bin numbers
        unsigned int n = m_central_eff->nBins();
        std::function<void (CollectionContainer*)> assign_mapping =  [this, &n](CollectionContainer* container){
                if (container != m_central_eff.get() && container->separateBinSyst()){
                    container->SetGlobalOffSet(n);
                    n += container->nBins();
                };
        };
        assign_mapping(m_calo_eff.get());
        assign_mapping(m_calo_eff.get());
        assign_mapping(m_lowpt_central_eff.get());
        assign_mapping(m_lowpt_calo_eff.get());
     
        assign_mapping(m_forward_eff.get());
       
        /// Systematic constructor has been called. We can now assemble
        /// the systematic variations
        if (!m_syst_set){
            m_syst_set = std::make_unique<CP::SystematicSet>();
            size_t glob_sys = m_ref_tool.getPosition(this);
            if (glob_sys > m_ref_tool.getNCollections()){
                Error("EffiCollection()", "Invalid position in the list of systematics. It seems that I'm not part of the referring ScaleFactorTool.");
                return false;
            }
            for (const EffiCollection::CollectionType& file_type: {EffiCollection::Central, EffiCollection::Calo, EffiCollection::Forward,  
                                         EffiCollection::CentralLowPt, EffiCollection::CaloLowPt}){
                    
                std::shared_ptr<CollectionContainer> container = retrieveContainer(file_type);
                if (container->isNominal()) continue;
                if (container->separateBinSyst()){
                    /// Let the world implode... Yeaha register foreach bin
                    /// a systematic variation                    
                    for (unsigned int b = container->nBins() - 1; b > 0  ; --b){
                        unsigned int bin = b + container->globalOffSet();
                        if (container->isOverFlowBin(bin)) continue;
                        m_syst_set->insert(CP::SystematicVariation::makeToyVariation("MUON_EFF_" + container->sysname()  + GetBinName(bin) , bin, glob_sys));
                    }
                } else {
                    m_syst_set->insert( SystematicVariation("MUON_EFF_" + container->sysname(), container->isUpVariation() ? 1  : -1));
                }
             }
        }
        return true;
    }

    CollectionContainer* EffiCollection::FindContainer(const xAOD::Muon& mu) const {
        if (mu.pt() <  m_ref_tool.lowPtTransition()) {
            if (std::abs(mu.eta()) >= 2.5) {
                return m_forward_eff.get();
            }
            if (mu.muonType() == xAOD::Muon::CaloTagged) {
                return m_lowpt_calo_eff.get();
            }
            return m_lowpt_central_eff.get();
        }
        if (mu.muonType() == xAOD::Muon::CaloTagged) {
            return m_calo_eff.get();
        } else if (std::abs(mu.eta()) < 2.5) {
            return m_central_eff.get();
        } else {
            return m_forward_eff.get();
        }
    }
    CollectionContainer* EffiCollection::FindContainer(unsigned int bin) const{
        if (m_central_eff->isBinInMap(bin)) return m_central_eff.get();
        if (m_forward_eff->isBinInMap(bin)) return m_forward_eff.get();
        if (m_calo_eff->isBinInMap(bin)) return m_calo_eff.get();
        if (m_lowpt_central_eff->isBinInMap(bin)) return m_lowpt_central_eff.get();
        if (m_lowpt_calo_eff->isBinInMap(bin)) return m_lowpt_calo_eff.get();
        return nullptr;
    }
            
    EfficiencyScaleFactor* EffiCollection::retrieveSF(const xAOD::Muon& mu, unsigned int RunNumber) const {
        CollectionContainer* Cont = FindContainer(mu);
        if (Cont != nullptr) return Cont->retrieve(RunNumber);
        Warning("EffiCollection::retrieveSF()", "Invalid muon");
        return nullptr;
    }
    unsigned int EffiCollection::nBins() const {
        unsigned int Nbins = 0;
        if (m_central_eff) {
            Nbins += m_central_eff->nBins();
        }
        if (m_central_eff != m_calo_eff) {
            Nbins += m_calo_eff->nBins();
        }
        if (m_forward_eff != m_central_eff) {
            Nbins += m_forward_eff->nBins();
        }
        if (m_lowpt_central_eff != m_central_eff) {
            Nbins += m_lowpt_central_eff->nBins();
        }
        if (m_lowpt_calo_eff != m_central_eff) {
            Nbins += m_lowpt_calo_eff->nBins();
        }
        return Nbins;
    }
    bool EffiCollection::SetSystematicBin(unsigned int Bin) {
        CollectionContainer* Cont = FindContainer(Bin);
        if (!Cont) return false;
        return Cont->SetSystematicBin(Bin);
    }
    bool EffiCollection::IsLowPtBin(unsigned int Bin) const {
        return (m_central_eff != m_lowpt_central_eff && m_lowpt_central_eff->isBinInMap(Bin)) ||
               (m_central_eff != m_lowpt_calo_eff && m_lowpt_calo_eff->isBinInMap(Bin));
    }
    bool EffiCollection::IsForwardBin(unsigned int Bin) const {
        return m_forward_eff != m_central_eff && m_forward_eff->isBinInMap(Bin);
    }
    
    std::string EffiCollection::FileTypeName(EffiCollection::CollectionType T) {
        if (T == CollectionType::Central) return "Central";
        if (T == CollectionType::Calo) return "Calo";
        if (T == CollectionType::Forward) return "Forward";
        if (T == CollectionType::CentralLowPt) return "CentralLowPt";
        if (T == CollectionType::CaloLowPt) return "CaloLowPt";
        return "EffiCollection::FileTypeName() - WARNING: Unknown EffiCollection::CollectionType!";
    }

    std::string EffiCollection::GetBinName(unsigned int bin) const {
        CollectionContainer* Cont = FindContainer(bin);
        if (Cont) {
            std::string BinName = FileTypeName(Cont->type()) +"_"+ Cont->GetBinName(bin);
            for (const std::string& R : ToRemove) {
                BinName = ReplaceExpInString(BinName, R, "");
            }
            for (const stringpair& R : ToReplace) {
                BinName = ReplaceExpInString(BinName, R.first, R.second);
            }
            return BinName;
        }
        Warning("EffiCollection::GetBinName()", "Unknown bin %u", bin);
        
        return "UNKNOWN_BIN";
    }
    int EffiCollection::getUnCorrelatedSystBin(const xAOD::Muon& mu) const {
        CollectionContainer* container = FindContainer(mu);
        if (container) return container->FindBinSF(mu); 
        return -1;
    }
    SystematicSet* EffiCollection::getSystSet() const{
        return m_syst_set.get();
    }
    bool  EffiCollection::isAffectedBySystematic(const SystematicVariation& variation) const {
        return m_syst_set->find(variation) != m_syst_set->end();
    }
    bool  EffiCollection::isAffectedBySystematic(const SystematicSet& set) const{
        if (set.empty()) return m_syst_set->empty();
        for (const SystematicVariation& variation: set){
            if (isAffectedBySystematic(variation)) return true;
        }
        return false;
    }            
    
    
    //################################################################################
    //                               CollectionContainer
    //################################################################################
    CollectionContainer::CollectionContainer(const MuonEfficiencyScaleFactors& ref_tool, EffiCollection::CollectionType FileType):
                m_SF(),
                m_currentSF(nullptr),
                m_FileType(FileType),
                m_binOffSet(0){
        
            std::map<std::string, std::pair<unsigned int, unsigned int>> map = findPeriods(ref_tool);
            for (auto& period : map) {
             m_SF.push_back(std::make_shared<EfficiencyScaleFactor>(ref_tool, fileName(ref_tool), period.first));
             m_SF.back()->setFirstLastRun(period.second.first, period.second.second);
         }
    }
    CollectionContainer::CollectionContainer(const MuonEfficiencyScaleFactors& ref_tool, CollectionContainer* Nominal, const std::string& syst_name, unsigned int syst_bit_map):
                m_SF(),
                m_currentSF(nullptr),
                m_FileType(Nominal->type()),
                m_binOffSet(0) {
            std::map<std::string, std::pair<unsigned int, unsigned int>> map = findPeriods(ref_tool);
            for (auto& period : map) {
                std::vector< std::shared_ptr<EfficiencyScaleFactor>>::const_iterator itr = std::find_if( Nominal->m_SF.begin(), 
                                                                                                      Nominal->m_SF.end(), 
                                                                                                      [&period](const std::shared_ptr<EfficiencyScaleFactor>& a){
                                                                                                            return a->coversRunNumber(period.second.first);
                                                                                                    });
                m_SF.push_back(std::make_shared<EfficiencyScaleFactor>(*itr, ref_tool, fileName(ref_tool), period.first, syst_name, syst_bit_map));
                m_SF.back()->setFirstLastRun(period.second.first, period.second.second);
         }
    }
    std::map<std::string, std::pair<unsigned int, unsigned int>> CollectionContainer::findPeriods(const MuonEfficiencyScaleFactors& ref_tool) const{
        
        std::string file_name = fileName(ref_tool);
        std::map<std::string, std::pair<unsigned int, unsigned int>> map;
      
        std::unique_ptr<TFile> fin (TFile::Open(file_name.c_str()));
        if (!fin || !fin->IsOpen()) {
            Error("CollectionContainer","Unable to open file %s", file_name.c_str());
            return map;
        }
        TTree* intree = 0;
        fin->GetObject("LumiData", intree);
        
        // if no Lumi tree is found, we assume that the SF are not binned in time
        if (!intree) {
            map["All"] = std::pair<unsigned int,unsigned int>(1, 999999);
        } else {
            std::string *period = 0;
            unsigned int firstRun = 0;
            unsigned int lastRun = 0;
            intree->SetBranchAddress("Period", &period);
            intree->SetBranchAddress("FirstRun", &firstRun);
            intree->SetBranchAddress("LastRun", &lastRun);
            for (int t = 0; intree->GetEntry(t); t++) {
                map[*period] = std::pair<unsigned int, unsigned int>(firstRun, lastRun);
            }
        }
        return map;
    }
    std::string CollectionContainer::fileName(const MuonEfficiencyScaleFactors& ref_tool) const{
        if (type() == EffiCollection::CollectionType::Central) return ref_tool.filename_Central();
        if (type() == EffiCollection::CollectionType::Calo) return ref_tool.filename_Calo();
        if (type() == EffiCollection::CollectionType::CentralLowPt) return ref_tool.filename_LowPt();
        if (type() == EffiCollection::CollectionType::CaloLowPt) return ref_tool.filename_LowPtCalo();
        return ref_tool.filename_HighEta();
    }
                  
    bool CollectionContainer::CheckConsistency()  {
        if (m_SF.empty()) {
            Error("CollectionContainer", "Could not retrieve any SFs from the input file");
            return false;
        }
        /// Check that there are no overlapping run numbers in the periods        
        for (std::vector< std::shared_ptr<EfficiencyScaleFactor>>::const_iterator first_sf = m_SF.begin() ; first_sf != m_SF.end(); ++first_sf)  {
            if (!(*first_sf)->CheckConsistency()) return false;
            for (std::vector< std::shared_ptr<EfficiencyScaleFactor>>::const_iterator second_sf = m_SF.begin(); second_sf != first_sf; ++second_sf) {
                if ( (*first_sf)->coversRunNumber( (*second_sf)->firstRun()) || (*first_sf)->coversRunNumber((*second_sf)->lastRun()) || 
                     (*second_sf)->coversRunNumber( (*first_sf)->firstRun()) || (*second_sf)->coversRunNumber((*first_sf)->lastRun())){
                    Error("CollectionContainer", "Overlapping periods observed in file type %s. As run %i is in period %i - %i. Please check your SF file!",  
                           EffiCollection::FileTypeName(m_FileType).c_str(), (*first_sf)->firstRun(), (*second_sf)->firstRun(), (*second_sf)->lastRun());
                    return false;
                }
            }
        }
        
        std::sort(m_SF.begin(), m_SF.end(), [](const std::shared_ptr<EfficiencyScaleFactor>& a, const std::shared_ptr<EfficiencyScaleFactor>& b){
            return a->firstRun() < b->firstRun();});
        return true;
    }
    bool CollectionContainer::LoadPeriod(unsigned int RunNumber) {
        if (!m_currentSF || !m_currentSF->coversRunNumber(RunNumber)) {
            for (auto& period : m_SF) {
                if (period->coversRunNumber(RunNumber)) {
                    m_currentSF = period.get();
                    return true;
                }
            }
        } else return true;
        Warning("CollectionContainer", "Could not find any SF period in %s matching the run number %u", EffiCollection::FileTypeName(type()).c_str(), RunNumber);
        return false;
    }
    EfficiencyScaleFactor* CollectionContainer::retrieve(unsigned int RunNumber) {
        if (!LoadPeriod(RunNumber)) {
            return (*m_SF.begin()).get();
        }
        return m_currentSF;
    }
    bool CollectionContainer::SetSystematicBin(unsigned int Bin) {
        for ( std::shared_ptr<EfficiencyScaleFactor>& Period : m_SF) {
            if (!Period->SetSystematicBin(Bin- m_binOffSet)) {
                return false;
            }
        }
        return true;
    }
    unsigned int CollectionContainer::nBins() const {
        return m_SF.empty() ? 0 : (*m_SF.begin())->nBins();
    }
    unsigned int CollectionContainer::nOverFlowBins() const {
          return m_SF.empty() ? 0 : (*m_SF.begin())->nOverFlowBins();
    }
    bool CollectionContainer::isOverFlowBin(int b) const {
          return m_SF.empty() ? true : (*m_SF.begin())->isOverFlowBin(b-m_binOffSet);
    }
    bool CollectionContainer::isBinInMap(unsigned int bin) const{
        return  m_binOffSet <= bin && bin < m_binOffSet + nBins();
    }
    EffiCollection::CollectionType CollectionContainer::type() const {
        return m_FileType;
    }
    std::string CollectionContainer::GetBinName(unsigned int Bin) const {
        return (*m_SF.begin())->GetBinName(Bin- m_binOffSet);
    }
    int CollectionContainer::FindBinSF(const xAOD::Muon &mu) const {
        if (m_SF.empty() ) return -1;
        int bin =  (*m_SF.begin())->FindBinSF(mu);
        return bin > 0 ? m_binOffSet + bin : bin;
    }
    void CollectionContainer::SetGlobalOffSet(unsigned int OffSet){
        m_binOffSet = OffSet;
    }
    unsigned int CollectionContainer::globalOffSet() const{
        return m_binOffSet;
    }
    bool CollectionContainer::isNominal() const {
        if (m_SF.empty()) return false;
        return (*m_SF.begin())->sysname(false).empty();
    }
    bool CollectionContainer::isUpVariation() const {
        if (m_SF.empty()) return false;
        return (*m_SF.begin())->IsUpVariation();        
    }
    bool CollectionContainer::separateBinSyst() const {
        if (m_SF.empty()) return false;
        return (*m_SF.begin())->separateBinSyst();
    }
    std::string CollectionContainer::sysname() const{
        if (m_SF.empty()) return "UNKNOWN SYST";
        return (*m_SF.begin())->sysname(false);
    }
    
}
