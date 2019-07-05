// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

using namespace edm;
using namespace reco;
using namespace std;

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>


//
// Class declaration
//

class eventNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit eventNtupler(const edm::ParameterSet&);

      ~eventNtupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void fillGENInfo(const edm::Event &iEvent);
      virtual void fillLHEInfo(const edm::Event &iEvent);
      virtual void endJob() override;

      edm::EDGetTokenT<GenParticleCollection> mcLabel_;
      edm::EDGetTokenT<GenEventInfoProduct> GenEventInfoToken;
      edm::EDGetTokenT< LHEEventProduct > LHEEventProductToken;

      TTree *DYTree;

      bool theKeepAllGen;

      vector<double> gen_phi;
      vector<double> gen_eta;
      vector<double> gen_pt;
      vector<double> gen_mass;
      vector<double> gen_charge;
      vector<int> gen_mother_index;
      vector<int> gen_status;
      vector<int> gen_PID;
      vector<int> gen_isPrompt;
      vector<int> gen_isPromptFinalState;
      vector<int> gen_isTauDecayProduct;
      vector<int> gen_isPromptTauDecayProduct;
      vector<int> gen_isDirectPromptTauDecayProductFinalState;
      vector<int> gen_isHardProcess;
      vector<int> gen_isLastCopy;
      vector<int> gen_isLastCopyBeforeFSR;
      vector<int> gen_isPromptDecayed;
      vector<int> gen_isDecayedLeptonHadron;
      vector<int> gen_fromHardProcessBeforeFSR;
      vector<int> gen_fromHardProcessDecayed;
      vector<int> gen_fromHardProcessFinalState;
      vector<int> gen_isMostlyLikePythia6Status3;
      double gen_weight;
      double genWeight_Q;
      double genWeight_X1;
      double genWeight_X2;
      int genWeight_id1;
      int genWeight_id2;
      double genWeight_alphaQCD;
      double genWeight_alphaQED;

      // LHE info
      vector<double> lhe_phi;
      vector<double> lhe_eta;
      vector<double> lhe_pt;
      vector<double> lhe_mass;
      vector<double> lhe_PID;
      vector<double> lhe_status;

      // ----------member data ---------------------------
};
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
eventNtupler::eventNtupler(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  usesResource("TFileService");
  //mcLabel_ = consumes<GenParticleCollection>(edm::InputTag("genParticles")); // use genParticles as in the input root file
  mcLabel_ = ( consumes< std::vector<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("genSrc")));
  GenEventInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  LHEEventProductToken = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));

  theKeepAllGen                     = iConfig.getUntrackedParameter<bool>("KeepAllGen", true);
}


eventNtupler::~eventNtupler()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

////////////////////////
// -- Get GEN info -- //
////////////////////////
void eventNtupler::fillGENInfo(const edm::Event &iEvent)
{
  //cout << "fill pdf info" << endl;
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(mcLabel_,genParticles);
  
  int counter=0;
  for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it , ++counter) {      

    if(!theKeepAllGen && counter > 30) continue;
    
    gen_PID.push_back( it->pdgId() );
    gen_pt.push_back( it->pt() );
    gen_mass.push_back( it->mass() );
    gen_charge.push_back( it->charge() );
    gen_eta.push_back( it->eta() );
    gen_phi.push_back( it->phi() );
    gen_status.push_back( it->status() );
    
    //Flags (Ref: https://indico.cern.ch/event/402279/contribution/5/attachments/805964/1104514/mcaod-Jun17-2015.pdf)
    gen_isPrompt.push_back( it->statusFlags().isPrompt() ); //not from hadron, muon or tau decay 
    gen_isPromptFinalState.push_back( it->isPromptFinalState() ); //isPrompt && final state (status==1)
    gen_isTauDecayProduct.push_back( it->statusFlags().isTauDecayProduct() ); //is directly or indirectly from a tau decay
    gen_isPromptTauDecayProduct.push_back( it->statusFlags().isPromptTauDecayProduct() ); //is directly or indirectly from a tau decay, where the tau did not come from a hadron decay
    gen_isDirectPromptTauDecayProductFinalState.push_back( it->isDirectPromptTauDecayProductFinalState() ); // is the direct decay product from a tau decay (ie no intermediate hadron), where the tau did not come from a hadron decay && final state
    gen_isHardProcess.push_back( it->isHardProcess() );
    gen_isLastCopy.push_back( it->isLastCopy() );
    gen_isLastCopyBeforeFSR.push_back( it->isLastCopyBeforeFSR() );
    gen_isPromptDecayed.push_back( it->isPromptDecayed() );
    gen_isDecayedLeptonHadron.push_back( it->statusFlags().isDecayedLeptonHadron() );
    gen_fromHardProcessBeforeFSR.push_back( it->fromHardProcessBeforeFSR() );
    gen_fromHardProcessDecayed.push_back( it->fromHardProcessDecayed() );
    gen_fromHardProcessFinalState.push_back( it->fromHardProcessFinalState() );
    gen_isMostlyLikePythia6Status3.push_back( it->fromHardProcessBeforeFSR() );
    
    int idx = -1;
    for( reco::GenParticleCollection::const_iterator mit = genParticles->begin(); mit != genParticles->end(); ++mit ){
      if( it->mother()==&(*mit) ){
        idx = std::distance(genParticles->begin(),mit);
        break;
      }
    }
    
    gen_mother_index.push_back( idx );
    
  }
   
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByToken(GenEventInfoToken, genEvtInfo);
  gen_weight = genEvtInfo->weight();

  genWeight_Q = genEvtInfo->pdf()->scalePDF;
  genWeight_X1 = genEvtInfo->pdf()->x.first;
  genWeight_X2 = genEvtInfo->pdf()->x.second;
  genWeight_id1 = genEvtInfo->pdf()->id.first;
  genWeight_id2 = genEvtInfo->pdf()->id.second;
  genWeight_alphaQCD = genEvtInfo->alphaQCD();
  genWeight_alphaQED = genEvtInfo->alphaQED();
  
}

void eventNtupler::fillLHEInfo(const edm::Event &iEvent){

	Handle<LHEEventProduct> LHEInfo;
	iEvent.getByToken(LHEEventProductToken, LHEInfo);

	const lhef::HEPEUP& lheEvent = LHEInfo->hepeup();
	std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;

	Int_t _nLHEParticle = 0;
	for( size_t idxParticle = 0; idxParticle < lheParticles.size(); ++idxParticle )
	{
		Int_t id = lheEvent.IDUP[idxParticle];

		//if( fabs(id) == 13 || fabs(id) == 11 || fabs(id) == 15 )
		if( ( 0 < fabs(id) && fabs(id) < 7 ) || ( 10 < fabs(id) && fabs(id) < 17 ) || ( 20 < fabs(id) && fabs(id) < 25 ) )
		{
			Double_t Px = lheParticles[idxParticle][0];
			Double_t Py = lheParticles[idxParticle][1];
			Double_t Pz = lheParticles[idxParticle][2];
			Double_t E = lheParticles[idxParticle][3];
			// Double_t M = lheParticles[idxParticle][4];		
			Int_t status = lheEvent.ISTUP[idxParticle];
			TLorentzVector temp_lhe;
			temp_lhe.SetPxPyPzE(Px,Py,Pz,E);
			//std::cout << "id: " << id << " status: " << status << " Pt: " << temp_lhe.Pt() << " mass: " << temp_lhe.M() << std::endl;

  			lhe_phi.push_back(temp_lhe.Phi());
  			lhe_eta.push_back(temp_lhe.Eta());
  			lhe_pt.push_back(temp_lhe.Pt());
  			lhe_mass.push_back(temp_lhe.M());
  			lhe_PID.push_back(id);
  			lhe_status.push_back(status);

			_nLHEParticle++;
		}
	}
	//nLHEParticle = _nLHEParticle;

	//// -- PDf weights for theoretical uncertainties: scale, PDF replica and alphaS variation -- //
	//// -- ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW -- //
	//double OriginalWeight = LHEInfo->originalXWGTUP();
	//// std::cout << "OriginalWeight: " << OriginalWeight << endl;
	//int nWeight = (int)LHEInfo->weights().size();
	//// std::cout << "nWeight: " << nWeight << endl;

	//for(int i=0; i<nWeight; i++)
	//{
	//	double weight = LHEInfo->weights()[i].wgt;
	//	double ratio = weight / OriginalWeight;
	//	PDFWeights.push_back( ratio );

	//	// std::cout << i << "th weight = " << weight << "(ID=" << LHEInfo->weights()[i].id <<"), ratio w.r.t. original: " << ratio << endl;
	//}


}

// ------------ method called for each event  ------------
void
eventNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ////////////initialize/////////////


  gen_phi.clear();
  gen_eta.clear();
  gen_pt.clear();
  gen_mass.clear();
  gen_charge.clear();
  gen_mother_index.clear();
  gen_status.clear();
  gen_PID.clear();
  gen_isPrompt.clear();
  gen_isPromptFinalState.clear();
  gen_isTauDecayProduct.clear();
  gen_isPromptTauDecayProduct.clear();
  gen_isDirectPromptTauDecayProductFinalState.clear();
  gen_isHardProcess.clear();
  gen_isLastCopy.clear();
  gen_isLastCopyBeforeFSR.clear();
  gen_isPromptDecayed.clear();
  gen_isDecayedLeptonHadron.clear();
  gen_fromHardProcessBeforeFSR.clear();
  gen_fromHardProcessDecayed.clear();
  gen_fromHardProcessFinalState.clear();
  gen_isMostlyLikePythia6Status3.clear();
  gen_weight=-999;
  genWeight_Q=-999;
  genWeight_X1=-999;
  genWeight_X2=-999;
  genWeight_id1=-999;
  genWeight_id2=-999;
  genWeight_alphaQCD=-999;
  genWeight_alphaQED=-999;

  // clear LHE info
  lhe_phi.clear();
  lhe_eta.clear();
  lhe_pt.clear();
  lhe_mass.clear();
  lhe_PID.clear();
  lhe_status.clear();

  fillGENInfo(iEvent);
  fillLHEInfo(iEvent);

  DYTree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void 
eventNtupler::beginJob()
{

    edm::Service<TFileService> fs;
    DYTree = fs->make<TTree>("SKFlat","SKFlat");

    DYTree->Branch("gen_phi", "vector<double>", &gen_phi);
    DYTree->Branch("gen_eta", "vector<double>", &gen_eta);
    DYTree->Branch("gen_pt", "vector<double>", &gen_pt);
    DYTree->Branch("gen_mass", "vector<double>", &gen_mass);
    DYTree->Branch("gen_charge", "vector<double>", &gen_charge);
    DYTree->Branch("gen_mother_index", "vector<int>", &gen_mother_index);
    DYTree->Branch("gen_status", "vector<int>", &gen_status);
    DYTree->Branch("gen_PID", "vector<int>", &gen_PID);
    DYTree->Branch("gen_isPrompt", "vector<int>", &gen_isPrompt);
    DYTree->Branch("gen_isPromptFinalState", "vector<int>", &gen_isPromptFinalState);
    DYTree->Branch("gen_isTauDecayProduct", "vector<int>", &gen_isTauDecayProduct);
    DYTree->Branch("gen_isPromptTauDecayProduct", "vector<int>", &gen_isPromptTauDecayProduct);
    DYTree->Branch("gen_isDirectPromptTauDecayProductFinalState", "vector<int>", &gen_isDirectPromptTauDecayProductFinalState);
    DYTree->Branch("gen_isHardProcess", "vector<int>", &gen_isHardProcess);
    DYTree->Branch("gen_isLastCopy", "vector<int>", &gen_isLastCopy);
    DYTree->Branch("gen_isLastCopyBeforeFSR", "vector<int>", &gen_isLastCopyBeforeFSR);
    DYTree->Branch("gen_isPromptDecayed", "vector<int>", &gen_isPromptDecayed);
    DYTree->Branch("gen_isDecayedLeptonHadron", "vector<int>", &gen_isDecayedLeptonHadron);
    DYTree->Branch("gen_fromHardProcessBeforeFSR", "vector<int>", &gen_fromHardProcessBeforeFSR);
    DYTree->Branch("gen_fromHardProcessDecayed", "vector<int>", &gen_fromHardProcessDecayed);
    DYTree->Branch("gen_fromHardProcessFinalState", "vector<int>", &gen_fromHardProcessFinalState);
    DYTree->Branch("gen_isMostlyLikePythia6Status3", "vector<int>", &gen_isMostlyLikePythia6Status3);
    DYTree->Branch("gen_weight", &gen_weight, "gen_weight/D");
    DYTree->Branch("genWeight_Q", &genWeight_Q, "genWeight_Q/D");
    DYTree->Branch("genWeight_X1", &genWeight_X1, "genWeight_X1/D");
    DYTree->Branch("genWeight_X2", &genWeight_X2, "genWeight_X2/D");
    DYTree->Branch("genWeight_id1", &genWeight_id1, "genWeight_id1/I");
    DYTree->Branch("genWeight_id2", &genWeight_id2, "genWeight_id2/I");
    DYTree->Branch("genWeight_alphaQCD", &genWeight_alphaQCD, "genWeight_alphaQCD/D");
    DYTree->Branch("genWeight_alphaQED", &genWeight_alphaQED, "genWeight_alphaQED/D");

    DYTree->Branch("lhe_phi", "vector<double>", &lhe_phi);
    DYTree->Branch("lhe_eta", "vector<double>", &lhe_eta);
    DYTree->Branch("lhe_pt", "vector<double>", &lhe_pt);
    DYTree->Branch("lhe_mass", "vector<double>", &lhe_mass);
    DYTree->Branch("lhe_PID", "vector<double>", &lhe_PID);
    DYTree->Branch("lhe_status", "vector<double>", &lhe_status);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
eventNtupler::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
eventNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(eventNtupler);
