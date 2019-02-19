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

bool debug=1;

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

TEventList* elist_WKn=new TEventList("elist_WKn");
z->Draw(">>elist_WKn","PTOFamps>0.05e-6&&PTOFamps<0.5e-6");//Cut from 50nA to 500 nA

double nPTWKmax_50nA500nA_N50=0;//number of ROI events in current 50 event block
double nPTWKmax_N_N50=0;//actual number of events in current 50 event block. 
//This should be 50 for all but the last ~50 events in the series

TBranch *bnPTWKmax_50nA500nA_N50 = tree->Branch("nPTWKmax_50nA500nA_N50",&nPTWKmax_50nA500nA_N50,"nPTWKmax_50nA500nA_N50/D");
TBranch *bnPTWKmax_N_N50 = tree->Branch("nPTWKmax_N_N50",&nPTWKmax_N_N50,"nPTWKmax_N_N50/D");

if(debug){cout<<"block tree"<<endl;}
for(iEntry=0;iEntry<Nentries;iEntry++){
	tree->GetEntry(iEntry);
	nPTWKmax_N_N50++;
	if(elist_WKn->Contains(iEntry)){nPTWKmax_50nA500nA_N50++;}

	if(debug){cout<<"iEntry:"<<iEntry<<", "<<nPTWKmax_50nA500nA_N50<<"/"<<nPTWKmax_N_N50<<endl;}

	if((iEntry+1)%50==0){
		//At the end of a block, go back and fill entries for the block
		for(int jEntry=iEntry-49;jEntry<=iEntry;jEntry++){
			tree->GetEntry(jEntry);
			bnPTWKmax_50nA500nA_N50->Fill();
			bnPTWKmax_N_N50->Fill();
		}
		nPTWKmax_50nA500nA_N50=0;
		nPTWKmax_N_N50=0;
	}
}
if(debug){cout<<"iEntry:"<<iEntry<<", "<<nPTWKmax_50nA500nA_N50<<"/"<<nPTWKmax_N_N50<<endl;}
//Still need to do the last <50 events if Nentries%50!=0 
if(Nentries%50!=0){
	for(int jEntry=Nentries-Nentries%50;jEntry<Nentries;jEntry++){
		if(debug){cout<<"jEntry:"<<jEntry<<", "<<nPTWKmax_50nA500nA_N50<<"/"<<nPTWKmax_N_N50<<endl;}
		tree->GetEntry(jEntry);
		bnPTWKmax_50nA500nA_N50->Fill();
		bnPTWKmax_N_N50->Fill();
	}
}

////////////////
//Save and Close
TFile hFile(outFile, "RECREATE");
tree->Write();
hFile.Close();

return 0;

}
