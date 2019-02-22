//ROOT Code to calculate trigger rates at each EventTime and make a tree with those values.
//
//N.Mast 01/2019

#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <algorithm>

//ROOT
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TObjArray.h"
#include "TPRegexp.h"
#include "TTree.h"
#include "TChain.h"
#include "TLeaf.h"
#include "TFile.h"
#include "Riostream.h"
#include "TTimeStamp.h"
#include "TCut.h"
#include "TEventList.h"

bool debug=0;

//Input params are path to a series's rq folder and output file name
//The parsing is stupid, so don't mess it up!
int main(int argc, char *argv[]){
if(argc!=3){
	cout << "buildRateTree Error: incorrect number of inputs\n";
	cout << "Usage is: ./buildRateTree [rqSeriesDir] [outputFile]\n";
	for(int iarg=0;iarg<argc;iarg++){
		cout<<iarg<<":"<<argv[iarg]<<endl;
	}
	return -1;
}

TString rqDir=TString(argv[1]);
TString outFile=TString(argv[2]);

if(debug){
	cout<<"rqDir:"<<rqDir<<endl;
	cout<<"outFile:"<<outFile<<endl;
}

//Open rqTree
TChain *z=new TChain("rqDir/zip1");
TChain *e=new TChain("rqDir/eventTree");
z->AddFriend(e);

z->Add(Form("%s/umn*.root",rqDir.Data()));
e->Add(Form("%s/umn*.root",rqDir.Data()));

if(z->GetEntries()==0){
	cout << "makeK100logTree Error: no events found\n";
	return -1;
}

if(debug){cout<<"z->GetEntries()="<<z->GetEntries()<<endl;}

//////////////////
//Initialize tree
double rate=0;
double rateTrigPTOF_50nA200nA=0, rateTrigPTOF_50nA200nA_Next=0, rateTrigPTOF_50nA200nA_Prev=0;
double rateTrigPTWKmax_50nA200nA=0, rateTrigPTWKmax_50nA200nA_Next=0, rateTrigPTWKmax_50nA200nA_Prev=0;
double rateTrigPTWKmax_200nA350nA=0, rateTrigPTWKmax_200nA350nA_Next=0, rateTrigPTWKmax_200nA350nA_Prev=0;

TTree *tree = new TTree("rateTree","Rate information");
//Rate of events saved in data in the current second given by EventTime. Includes triggered and randoms
tree->Branch("rate",&rate,"rate/D");

////////////////PTOFamps///////////////////////
//Rate of triggered events with 50e-9<PTOFamps<200e-9 in the current second
tree->Branch("rateTrigPTOF_50nA200nA",&rateTrigPTOF_50nA200nA,"rateTrigPTOF_50nA200nA/D");
//Rate of triggered events with 50e-9<PTOFamps<200e-9 in the NEXT second
tree->Branch("rateTrigPTOF_50nA200nA_Next",&rateTrigPTOF_50nA200nA_Next,"rateTrigPTOF_50nA200nA_Next/D");
//Rate of triggered events with 50e-9<PTOFamps<200e-9 in the PREVIOUS second
tree->Branch("rateTrigPTOF_50nA200nA_Prev",&rateTrigPTOF_50nA200nA_Prev,"rateTrigPTOF_50nA200nA_Prev/D");

////////////////PTWKmax///////////////////////
///////////
//50-200 nA

//Rate of triggered events with 50e-9<PTWKmax<200e-9 in the current second
tree->Branch("rateTrigPTWKmax_50nA200nA",&rateTrigPTWKmax_50nA200nA,"rateTrigPTWKmax_50nA200nA/D");
//Rate of triggered events with 50e-9<PTWKmax<200e-9 in the NEXT second
tree->Branch("rateTrigPTWKmax_50nA200nA_Next",&rateTrigPTWKmax_50nA200nA_Next,"rateTrigPTWKmax_50nA200nA_Next/D");
//Rate of triggered events with 50e-9<PTWKmax<200e-9 in the PREVIOUS second
tree->Branch("rateTrigPTWKmax_50nA200nA_Prev",&rateTrigPTWKmax_50nA200nA_Prev,"rateTrigPTWKmax_50nA200nA_Prev/D");

///////////
//200-350 nA

//Rate of triggered events with 200e-9<PTWKmax<350e-9 in the current second
tree->Branch("rateTrigPTWKmax_200nA350nA",&rateTrigPTWKmax_200nA350nA,"rateTrigPTWKmax_200nA350nA/D");
//Rate of triggered events with 200e-9<PTWKmax<350e-9 in the NEXT second
tree->Branch("rateTrigPTWKmax_200nA350nA_Next",&rateTrigPTWKmax_200nA350nA_Next,"rateTrigPTWKmax_200nA350nA_Next/D");
//Rate of triggered events with 200e-9<PTWKmax<350e-9 in the PREVIOUS second
tree->Branch("rateTrigPTWKmax_200nA350nA_Prev",&rateTrigPTWKmax_200nA350nA_Prev,"rateTrigPTWKmax_200nA350nA_Prev/D");


double EventTime=0;
z->SetBranchAddress("EventTime",&EventTime);

//Identify low energy events
TCut crand = "EventCategory";

TEventList* elist_OF=new TEventList("elist_OF");
z->Draw(">>elist_OF",!crand&&"PTOFamps>0.05e-6&&PTOFamps<0.2e-6");//Cut from 50nA to 200 nA

TEventList* elist_WK1=new TEventList("elist_WK1");
z->Draw(">>elist_WK1",!crand&&"PTWKmax>0.05e-6&&PTWKmax<0.2e-6");//Cut from 50nA to 200 nA

TEventList* elist_WK2=new TEventList("elist_WK2");
z->Draw(">>elist_WK2",!crand&&"PTWKmax>0.2e-6&&PTWKmax<0.35e-6");//Cut from 200nA to 350 nA


if(debug){cout<<"elist_OF->GetN():"<<elist_OF->GetN()<<endl;}
if(debug){cout<<"elist_WK1->GetN():"<<elist_WK1->GetN()<<endl;}
if(debug){cout<<"elist_WK2->GetN():"<<elist_WK2->GetN()<<endl;}

////////////
//Initialize
int Nentries=z->GetEntries();
int iEntry_tailPrev=0, iEntry_tail=0, iEntry_tailNext=0;//first entries of the three time segments
int iEntry=0;//head index
z->GetEntry(iEntry);
double ETcurr=EventTime;//head time
while(EventTime==ETcurr&&iEntry<Nentries){//Parse first time block
	if(elist_OF->Contains(iEntry)){rateTrigPTOF_50nA200nA++;}
	if(elist_WK1->Contains(iEntry)){rateTrigPTWKmax_50nA200nA++;}
	if(elist_WK2->Contains(iEntry)){rateTrigPTWKmax_200nA350nA++;}
	iEntry++;
	z->GetEntry(iEntry);
}
//Mark the start of the 2nd block
iEntry_tailNext=iEntry;
ETcurr=EventTime;
while(EventTime==ETcurr&&iEntry<Nentries){//Parse 2nd time block
	if(elist_OF->Contains(iEntry)){rateTrigPTOF_50nA200nA_Next++;}
	if(elist_WK1->Contains(iEntry)){rateTrigPTWKmax_50nA200nA_Next++;}
	if(elist_WK2->Contains(iEntry)){rateTrigPTWKmax_200nA350nA_Next++;}
	iEntry++;
	z->GetEntry(iEntry);
}
rate=iEntry_tailNext-iEntry_tail;//total number of events in 1st block
//Fill tree for curr block
for(;iEntry_tail<iEntry_tailNext;iEntry_tail++){tree->Fill();}
//Now iEntry_tailPrev==0,iEntry_tail==iEntry_tailNext
iEntry_tailNext=iEntry;
ETcurr=EventTime;
rateTrigPTOF_50nA200nA_Prev=rateTrigPTOF_50nA200nA;
rateTrigPTOF_50nA200nA=rateTrigPTOF_50nA200nA_Next;
rateTrigPTOF_50nA200nA_Next=0;

rateTrigPTWKmax_50nA200nA_Prev=rateTrigPTWKmax_50nA200nA;
rateTrigPTWKmax_50nA200nA=rateTrigPTWKmax_50nA200nA_Next;
rateTrigPTWKmax_50nA200nA_Next=0;

rateTrigPTWKmax_200nA350nA_Prev=rateTrigPTWKmax_200nA350nA;
rateTrigPTWKmax_200nA350nA=rateTrigPTWKmax_200nA350nA_Next;
rateTrigPTWKmax_200nA350nA_Next=0;

///////////
//Main Loop
for(;iEntry<Nentries;iEntry++){
	z->GetEntry(iEntry);
	if(debug){cout<<"iEntry:"<<iEntry<<" EventTime:"<<EventTime<<endl;}
	if(EventTime==ETcurr){
		if(elist_OF->Contains(iEntry)){rateTrigPTOF_50nA200nA_Next++;}
		if(elist_WK1->Contains(iEntry)){rateTrigPTWKmax_50nA200nA_Next++;}
		if(elist_WK2->Contains(iEntry)){rateTrigPTWKmax_200nA350nA_Next++;}
	}else{//hit a new time block
		rate=iEntry_tailNext-iEntry_tail;//total number of events in curr block
		for(;iEntry_tail<iEntry_tailNext;iEntry_tail++){tree->Fill();}
		iEntry_tailNext=iEntry;
		ETcurr=EventTime;

		rateTrigPTOF_50nA200nA_Prev=rateTrigPTOF_50nA200nA;
		rateTrigPTOF_50nA200nA=rateTrigPTOF_50nA200nA_Next;
		rateTrigPTOF_50nA200nA_Next=0;
		if(elist_OF->Contains(iEntry)){rateTrigPTOF_50nA200nA_Next++;}

		rateTrigPTWKmax_50nA200nA_Prev=rateTrigPTWKmax_50nA200nA;
		rateTrigPTWKmax_50nA200nA=rateTrigPTWKmax_50nA200nA_Next;
		rateTrigPTWKmax_50nA200nA_Next=0;
		if(elist_WK1->Contains(iEntry)){rateTrigPTWKmax_50nA200nA_Next++;}

		rateTrigPTWKmax_200nA350nA_Prev=rateTrigPTWKmax_200nA350nA;
		rateTrigPTWKmax_200nA350nA=rateTrigPTWKmax_200nA350nA_Next;
		rateTrigPTWKmax_200nA350nA_Next=0;
		if(elist_WK2->Contains(iEntry)){rateTrigPTWKmax_200nA350nA_Next++;}

	}
}
///////////////////////
//Last 2 time segments
rate=iEntry_tailNext-iEntry_tail;//total number of events in curr block
for(;iEntry_tail<iEntry_tailNext;iEntry_tail++){tree->Fill();}
iEntry_tailNext=iEntry;
ETcurr=EventTime;

rateTrigPTOF_50nA200nA_Prev=rateTrigPTOF_50nA200nA;
rateTrigPTOF_50nA200nA=rateTrigPTOF_50nA200nA_Next;
rateTrigPTOF_50nA200nA_Next=0;

rateTrigPTWKmax_50nA200nA_Prev=rateTrigPTWKmax_50nA200nA;
rateTrigPTWKmax_50nA200nA=rateTrigPTWKmax_50nA200nA_Next;
rateTrigPTWKmax_50nA200nA_Next=0;

rateTrigPTWKmax_200nA350nA_Prev=rateTrigPTWKmax_200nA350nA;
rateTrigPTWKmax_200nA350nA=rateTrigPTWKmax_200nA350nA_Next;
rateTrigPTWKmax_200nA350nA_Next=0;

rate=iEntry_tailNext-iEntry_tail;//total number of events in final block
for(;iEntry_tail<iEntry_tailNext;iEntry_tail++){tree->Fill();}


//////////////////////////////////////////
//A second loop to find number of hits
//  in blocks with fixed number of events

TEventList* elist_WKn0=new TEventList("elist_WKn0");
z->Draw(">>elist_WKn0","PTOFamps>0.05e-6&&PTOFamps<0.2e-6");//Cut from 50nA to 200 nA
TEventList* elist_WKn1=new TEventList("elist_WKn1");
z->Draw(">>elist_WKn1","PTOFamps>0.2e-6&&PTOFamps<0.5e-6");//Cut from 200nA to 500 nA
TEventList* elist_WKn2=new TEventList("elist_WKn2");
z->Draw(">>elist_WKn2","PTOFamps>0.5e-6&&PTOFamps<1e-6");//Cut from 500nA to 1 uA

//N=50
double nPTWKmax_50nA200nA_N50=0;//number of ROI0 events in current 50 event block
double nPTWKmax_200nA500nA_N50=0;//number of ROI1 events in current 50 event block
double nPTWKmax_500nA1uA_N50=0;//number of ROI2 events in current 50 event block
double nPTWKmax_N_N50=0;//actual number of events in current 50 event block. 
//This should be 50 for all but the last ~50 events in the series
//N=100
double nPTWKmax_50nA200nA_N100=0;//number of ROI0 events in current 100 event block
double nPTWKmax_200nA500nA_N100=0;//number of ROI1 events in current 100 event block
double nPTWKmax_500nA1uA_N100=0;//number of ROI2 events in current 100 event block
double nPTWKmax_N_N100=0;//actual number of events in current 100 event block. 
//This should be 100 for all but the last ~100 events in the series
//N=500
double nPTWKmax_50nA200nA_N500=0;//number of ROI0 events in current 500 event block
double nPTWKmax_200nA500nA_N500=0;//number of ROI1 events in current 500 event block
double nPTWKmax_500nA1uA_N500=0;//number of ROI2 events in current 500 event block
double nPTWKmax_N_N500=0;//actual number of events in current 500 event block. 
//This should be 500 for all but the last ~500 events in the series


TBranch *bnPTWKmax_50nA200nA_N50 = tree->Branch("nPTWKmax_50nA200nA_N50",&nPTWKmax_50nA200nA_N50,"nPTWKmax_50nA200nA_N50/D");
TBranch *bnPTWKmax_200nA500nA_N50 = tree->Branch("nPTWKmax_200nA500nA_N50",&nPTWKmax_200nA500nA_N50,"nPTWKmax_200nA500nA_N50/D");
TBranch *bnPTWKmax_500nA1uA_N50 = tree->Branch("nPTWKmax_500nA1uA_N50",&nPTWKmax_500nA1uA_N50,"nPTWKmax_500nA1uA_N50/D");
TBranch *bnPTWKmax_N_N50 = tree->Branch("nPTWKmax_N_N50",&nPTWKmax_N_N50,"nPTWKmax_N_N50/D");

TBranch *bnPTWKmax_50nA200nA_N100 = tree->Branch("nPTWKmax_50nA200nA_N100",&nPTWKmax_50nA200nA_N100,"nPTWKmax_50nA200nA_N100/D");
TBranch *bnPTWKmax_200nA500nA_N100 = tree->Branch("nPTWKmax_200nA500nA_N100",&nPTWKmax_200nA500nA_N100,"nPTWKmax_200nA500nA_N100/D");
TBranch *bnPTWKmax_500nA1uA_N100 = tree->Branch("nPTWKmax_500nA1uA_N100",&nPTWKmax_500nA1uA_N100,"nPTWKmax_500nA1uA_N100/D");
TBranch *bnPTWKmax_N_N100 = tree->Branch("nPTWKmax_N_N100",&nPTWKmax_N_N100,"nPTWKmax_N_N100/D");

TBranch *bnPTWKmax_50nA200nA_N500 = tree->Branch("nPTWKmax_50nA200nA_N500",&nPTWKmax_50nA200nA_N500,"nPTWKmax_50nA200nA_N500/D");
TBranch *bnPTWKmax_200nA500nA_N500 = tree->Branch("nPTWKmax_200nA500nA_N500",&nPTWKmax_200nA500nA_N500,"nPTWKmax_200nA500nA_N500/D");
TBranch *bnPTWKmax_500nA1uA_N500 = tree->Branch("nPTWKmax_500nA1uA_N500",&nPTWKmax_500nA1uA_N500,"nPTWKmax_500nA1uA_N500/D");
TBranch *bnPTWKmax_N_N500 = tree->Branch("nPTWKmax_N_N500",&nPTWKmax_N_N500,"nPTWKmax_N_N500/D");

if(debug){cout<<"block tree"<<endl;}
for(iEntry=0;iEntry<Nentries;iEntry++){
	tree->GetEntry(iEntry);
	//Increment block counters
	nPTWKmax_N_N50++;
	nPTWKmax_N_N100++;
	nPTWKmax_N_N500++;
	//Increment ROI counters
	if(elist_WKn0->Contains(iEntry)){
		nPTWKmax_50nA200nA_N50++;
		nPTWKmax_50nA200nA_N100++;
		nPTWKmax_50nA200nA_N500++;
	}
	if(elist_WKn1->Contains(iEntry)){
		nPTWKmax_200nA500nA_N50++;
		nPTWKmax_200nA500nA_N100++;
		nPTWKmax_200nA500nA_N500++;
	}
	if(elist_WKn2->Contains(iEntry)){
		nPTWKmax_500nA1uA_N50++;
		nPTWKmax_500nA1uA_N100++;
		nPTWKmax_500nA1uA_N500++;
	}

	if(debug){
		cout<<"iEntry:"<<iEntry<<", "<<nPTWKmax_200nA500nA_N50<<"/"<<nPTWKmax_N_N50;
		cout<<", "<<nPTWKmax_200nA500nA_N100<<"/"<<nPTWKmax_N_N100;
		cout<<", "<<nPTWKmax_200nA500nA_N500<<"/"<<nPTWKmax_N_N500<<endl;
	}

	//Handle the end of each block

	//NOTE! These blocks must be handled in order of increasing size
	//because moving back and forth in the tree can overwrite variables 
	//in memory with their values somewhere in the tree
	//That's also why the counts only get reset after handling the writing of all blocks

	//End of a 50 block
	if(nPTWKmax_N_N50==50){
		//At the end of a 50 block, go back and fill entries for the block
		for(int jEntry=iEntry-49;jEntry<=iEntry;jEntry++){
			tree->GetEntry(jEntry);
			bnPTWKmax_50nA200nA_N50->Fill();
			bnPTWKmax_200nA500nA_N50->Fill();
			bnPTWKmax_500nA1uA_N50->Fill();
			bnPTWKmax_N_N50->Fill();
		}
		//Counters reset below
	}

	//End of a 100 block
	if(nPTWKmax_N_N100==100){
		for(int jEntry=iEntry-99;jEntry<=iEntry;jEntry++){
			tree->GetEntry(jEntry);
			bnPTWKmax_50nA200nA_N100->Fill();
			bnPTWKmax_200nA500nA_N100->Fill();
			bnPTWKmax_500nA1uA_N100->Fill();
			bnPTWKmax_N_N100->Fill();
		}
		//Counters reset below
	}

	//End of a 500 block
	if(nPTWKmax_N_N500==500){
		//At the end of a 500 block, go back and fill entries for the block
		for(int jEntry=iEntry-499;jEntry<=iEntry;jEntry++){
			tree->GetEntry(jEntry);
			bnPTWKmax_50nA200nA_N500->Fill();
			bnPTWKmax_200nA500nA_N500->Fill();
			bnPTWKmax_500nA1uA_N500->Fill();
			bnPTWKmax_N_N500->Fill();
		}
	}



	//Reset counters
	if(nPTWKmax_N_N50==50){
		nPTWKmax_50nA200nA_N50=0;
		nPTWKmax_200nA500nA_N50=0;
		nPTWKmax_500nA1uA_N50=0;
		nPTWKmax_N_N50=0;
	}
	if(nPTWKmax_N_N100==100){
		nPTWKmax_50nA200nA_N100=0;
		nPTWKmax_200nA500nA_N100=0;
		nPTWKmax_500nA1uA_N100=0;
		nPTWKmax_N_N100=0;
	}
	if(nPTWKmax_N_N500==500){
		nPTWKmax_50nA200nA_N500=0;
		nPTWKmax_200nA500nA_N500=0;
		nPTWKmax_500nA1uA_N500=0;
		nPTWKmax_N_N500=0;
	}

}
if(debug){cout<<"iEntry:"<<iEntry<<", "<<nPTWKmax_200nA500nA_N50<<"/"<<nPTWKmax_N_N50<<endl;}

//Still need to do the last <50 events if Nentries%50!=0 
//last N<50
if(Nentries%50!=0){
	for(int jEntry=Nentries-Nentries%50;jEntry<Nentries;jEntry++){
		if(debug){cout<<"jEntry:"<<jEntry<<", "<<nPTWKmax_200nA500nA_N50<<"/"<<nPTWKmax_N_N50<<endl;}
		tree->GetEntry(jEntry);
		bnPTWKmax_50nA200nA_N50->Fill();
		bnPTWKmax_200nA500nA_N50->Fill();
		bnPTWKmax_500nA1uA_N50->Fill();
		bnPTWKmax_N_N50->Fill();
	}
}
//last N<100
if(Nentries%100!=0){
	for(int jEntry=Nentries-Nentries%100;jEntry<Nentries;jEntry++){
		tree->GetEntry(jEntry);
		bnPTWKmax_50nA200nA_N100->Fill();
		bnPTWKmax_200nA500nA_N100->Fill();
		bnPTWKmax_500nA1uA_N100->Fill();
		bnPTWKmax_N_N100->Fill();
	}
}
//last N<500
if(Nentries%500!=0){
	for(int jEntry=Nentries-Nentries%500;jEntry<Nentries;jEntry++){
		tree->GetEntry(jEntry);
		bnPTWKmax_50nA200nA_N500->Fill();
		bnPTWKmax_200nA500nA_N500->Fill();
		bnPTWKmax_500nA1uA_N500->Fill();
		bnPTWKmax_N_N500->Fill();
	}
}

////////////////
//Save and Close
TFile hFile(outFile, "RECREATE");
tree->Write();
hFile.Close();

return 0;

}
