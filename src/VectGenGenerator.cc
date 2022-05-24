/**
 * @file VectGenGenerator.cc
 *
 * @date 2022-03-09
 * @author Y.Koshio
 *
 * @history
 *  21 APR 2022 S.Izumiyama: integrated flux table class
 */

#include <math.h>

#include "VectGenSnConst.hh"
#include "VectGenGenerator.hh"
#include "VectGenSetBin.hh"

#include "VectGenUtil.hh"

#include "FluxCalculation.hh"

void VectGenGenerator::convDirection( const double eTheta, const double ePhi,
				      double * eDir )
{
	double origVec[3];
	origVec[0] = sin( eTheta ) * cos( ePhi );
	origVec[1] = sin( eTheta ) * sin( ePhi );
	origVec[2] = cos( eTheta );

	eDir[0] = 0.;
	eDir[1] = 0.;
	eDir[2] = 0.;

	for( int i = 0; i < 3; i++ ){
		for( int j = 0; j < 3; j++ ){
			eDir[i] += Rmat[i][j] * origVec[j];
		}
	}

	return;
}

void VectGenGenerator::determineAngleNuebarP( const double nuEne,
					      double & eEne, double & eTheta, double & ePhi )
{

	double nuEnergy = nuEne;
	double cost, eMom, p;

	//find maximum values, which depends on nuEne
	double maxP = 0.;
	for( int iCost = 0; iCost < costNBins; iCost++ ){
		cost = costMin + costBinSize * ( double(iCost) + 0.5 );
		nucrs->DcsNuebP_SV(nuEnergy, cost, eEne, p);

		if( p > maxP ){ maxP = p; }
	}

	while( 1 ){
		cost = getRandomReal( costMin, costMax, generator );
		nucrs->DcsNuebP_SV(nuEnergy, cost, eEne, p);

		double x = getRandomReal( 0., maxP, generator );
		if( x < p ){
			eTheta = acos( cost );
			ePhi = getRandomReal( -M_PI,  M_PI, generator );
			break;
		}
	}

	return;

}

void VectGenGenerator::determineAngleElastic( const int nReact, const double nuEne
		, double & eEne, double & eTheta, double & ePhi )
{

	double nuEnergy = nuEne;
	double cost, eEnergy;
	double p=0., x;

	//we know the maximum prob. happens at cost=1
	cost = 1. - zero_precision;
	eEnergy = nucrs->calcElectronTotEnergyElastic( nuEnergy, cost );

	double maxP = 0.;
	if( nReact == 1 ) maxP = sl_nue_dif_rad_( & nuEnergy, & eEnergy);
	if( nReact == 2 ) maxP = sl_neb_dif_rad_( & nuEnergy, & eEnergy);
	if( nReact == 3 ) maxP = sl_num_dif_rad_( & nuEnergy, & eEnergy);
	if( nReact == 4 ) maxP = sl_nmb_dif_rad_( & nuEnergy, & eEnergy);
	maxP *= nucrs->calcDeEneDcostElastic( nuEnergy, cost );

	/*
	//take into account the electron energy threshold
	//double costTh = nucrs->calcCosTthElastic( nuEnergy, eEneThr );
	*/

	// avoid too low energy event
	double eEneThrElastic = Me + zero_precision ;
	double costTh = nucrs->calcCosTthElastic( nuEnergy, eEneThrElastic );

	if(fabs(costTh) > 1.){ // neutrino energy is too low
	  std::cerr << "Neutrino energy is too low : neutrino energy " << nuEnergy << std::endl;
	  iSkip = 1;
	  return;
	}

	while( 1 ){
		cost = getRandomReal( costTh, 1., generator);
		eEnergy = nucrs->calcElectronTotEnergyElastic( nuEnergy, cost );

		if( nReact == 1 ) p = sl_nue_dif_rad_( & nuEnergy, & eEnergy);
		if( nReact == 2 ) p = sl_neb_dif_rad_( & nuEnergy, & eEnergy);
		if( nReact == 3 ) p = sl_num_dif_rad_( & nuEnergy, & eEnergy);
		if( nReact == 4 ) p = sl_nmb_dif_rad_( & nuEnergy, & eEnergy);
		p *= nucrs->calcDeEneDcostElastic( nuEnergy, cost );

		x = getRandomReal( 0., maxP, generator );
		if( x < p ){
			eEne = eEnergy;
			eTheta = acos( cost );
			ePhi = getRandomReal( -M_PI, M_PI, generator );
			break;
		}
	}

}

void VectGenGenerator::determineKinematics( const int nReact, const double nuEne, double * snDir, MCInfo * mc )
{

	if( nReact == 0 ){ // nuebar + p -> e+ + n
		mc->nvc = 4;

		// Original neutrino
		mc->ipvc[0] = -12; // anti-electron neutrino
		mc->energy[0] = nuEne; // ENERGY ( MEV )
		mc->pvc[0][0] = nuEne * snDir[0]; // MOMENTUM OF I-TH PARTICLE ( MEV/C )
		mc->pvc[0][1] = nuEne * snDir[1];
		mc->pvc[0][2] = nuEne * snDir[2];
		mc->iorgvc[0] = 0;  // ID OF ORIGIN PARTICLE  PARENT PARTICLE
		mc->ivtivc[0] = 0;  // VERTEX # ( INITIAL )
		mc->iflgvc[0] = -1; // FINAL STATE FLAG
		mc->icrnvc[0] = 0;  // CHERENKOV FLAG
		mc->ivtfvc[0] = 1;  // VERTEX # ( FINAL )

		// Original proton
		mc->ipvc[1] = 2212; // anti-electron neutrino
		mc->energy[1] = Mp; // ENERGY ( MEV )
		mc->pvc[1][0] = 0.; // MOMENTUM OF I-TH PARTICLE ( MEV/C )
		mc->pvc[1][1] = 0.;
		mc->pvc[1][2] = 0.;
		mc->iorgvc[1] = 0;  // ID OF ORIGIN PARTICLE  PARENT PARTICLE
		mc->ivtivc[1] = 0;  // VERTEX # ( INITIAL )
		mc->iflgvc[1] = -1; // FINAL STATE FLAG
		mc->icrnvc[1] = 0;  // CHERENKOV FLAG
		mc->ivtfvc[1] = 1;  // VERTEX # ( FINAL )

		// Positron
		double eEne, eTheta, ePhi, eDir[3];
		determineAngleNuebarP( nuEne, eEne, eTheta, ePhi );
		double amom = sqrt(SQ( eEne ) - SQ( Me ));
		convDirection( eTheta, ePhi, eDir );
		//std::cout << nuEne << " " << eEne << " " << eDir[0] << " " << eDir[1] << " " << eDir[2] << std::endl;

		mc->ipvc[2] = -11; // positron
		mc->energy[2] = eEne; // total ENERGY ( MEV )
		mc->pvc[2][0] = amom * eDir[0]; // MOMENTUM OF I-TH PARTICLE ( MEV/C )
		mc->pvc[2][1] = amom * eDir[1];
		mc->pvc[2][2] = amom * eDir[2];
		mc->iorgvc[2] = 1;  // ID OF ORIGIN PARTICLE  PARENT PARTICLE
		mc->ivtivc[2] = 1;  // VERTEX # ( INITIAL )
		mc->iflgvc[2] = 0; // FINAL STATE FLAG
		mc->icrnvc[2] = 1;  // CHERENKOV FLAG
		mc->ivtfvc[2] = 1;  // VERTEX # ( FINAL )

		// Neutron
		mc->ipvc[3] = 2112; // neutron
		mc->pvc[3][0] = (mc->pvc[0][0]) - (mc->pvc[1][0]); // MOMENTUM OF I-TH PARTICLE ( MEV/C )
		mc->pvc[3][1] = (mc->pvc[0][1]) - (mc->pvc[1][1]);
		mc->pvc[3][2] = (mc->pvc[0][2]) - (mc->pvc[1][2]);
		mc->energy[3] = sqrt(SQ( mc->pvc[2][0] ) + SQ( mc->pvc[2][1] ) + SQ( mc->pvc[2][2] )  + SQ( Mn )); // ENERGY ( MEV )
		mc->iorgvc[3] = 1;  // ID OF ORIGIN PARTICLE  PARENT PARTICLE
		mc->ivtivc[3] = 1;  // VERTEX # ( INITIAL )
		mc->iflgvc[3] = 0; // FINAL STATE FLAG
		mc->icrnvc[3] = 1;  // CHERENKOV FLAG
		mc->ivtfvc[3] = 1;  // VERTEX # ( FINAL )

	}
	else if( nReact == 1 || nReact == 2 || nReact == 3 || nReact == 4 ){ //nu + e Elastic

		mc->nvc = 2;

		// Original neutrino
		if( nReact == 1) mc->ipvc[0] =  12;
		if( nReact == 2) mc->ipvc[0] = -12;
		if( nReact == 3) mc->ipvc[0] =  14;
		if( nReact == 4) mc->ipvc[0] = -14;
		mc->energy[0] = nuEne; // ENERGY ( MEV )
		mc->pvc[0][0] = nuEne * snDir[0]; // MOMENTUM OF I-TH PARTICLE ( MEV/C )
		mc->pvc[0][1] = nuEne * snDir[1];
		mc->pvc[0][2] = nuEne * snDir[2];
		mc->iorgvc[0] = 0;  // ID OF ORIGIN PARTICLE  PARENT PARTICLE
		mc->ivtivc[0] = 0;  // VERTEX # ( INITIAL )
		mc->iflgvc[0] = -1; // FINAL STATE FLAG
		mc->icrnvc[0] = 0;  // CHERENKOV FLAG
		mc->ivtfvc[0] = 1;  // VERTEX # ( FINAL )

		// Recoil electron
		double eEne, eTheta, ePhi, eDir[3];
		determineAngleElastic( nReact, nuEne, eEne, eTheta, ePhi );
		double amom = sqrt(SQ( eEne ) - SQ( Me ));
		convDirection( eTheta, ePhi, eDir );

		mc->ipvc[1] = 11; // electron
		mc->energy[1] = eEne; // total ENERGY ( MEV )
		mc->pvc[1][0] = amom * eDir[0]; // MOMENTUM OF I-TH PARTICLE ( MEV/C )
		mc->pvc[1][1] = amom * eDir[1];
		mc->pvc[1][2] = amom * eDir[2];
		mc->iorgvc[1] = 1;  // ID OF ORIGIN PARTICLE  PARENT PARTICLE
		mc->ivtivc[1] = 1;  // VERTEX # ( INITIAL )
		mc->iflgvc[1] = 0; // FINAL STATE FLAG
		mc->icrnvc[1] = 1;  // CHERENKOV FLAG
		mc->ivtfvc[1] = 1;  // VERTEX # ( FINAL )

		double costh = snDir[0] * eDir[0] + snDir[1] * eDir[1] + snDir[2] * eDir[2];
		//std::cout << "Elastic " << costh << std::endl;


	}
	else {
		std::cout << "Not Supported yet" << std::endl;
		exit(-1);
	}
	/*
	   }else if( nReact == 5 || nReact == 6 ){ //nue + O or nuebar +O
	   determineAngleNueO( nReact, nuEne, eEne, eTheta, ePhi );
	   pType = 11;
	   }else if( nReact == 7 || nReact == 8 || nReact == 9 || nReact == 10){// NC 
	   determineAngleOfGamma( nReact, eEne, eTheta, ePhi );
	   pType = 22;
	   }
	   */

	return;
}

void VectGenGenerator::determinePosition( double &x, double &y, double &z )
{
	//random inside the full tank (32.5kton)
	double r2 = getRandomReal( 0., 1., generator ) * RINTK * RINTK ;
	double r = sqrt( r2 );
	double phi = getRandomReal( 0., 1., generator ) * 2. * M_PI;
	x = r * cos( phi );
	y = r * sin( phi );
	z = ZMINTK + getRandomReal( 0., 1., generator ) * ( ZPINTK - ZMINTK );
	return;
}

void VectGenGenerator::FillEvent()
{

	/*---- Time sorting ----*/
	sort( vEvtInfo.begin(), vEvtInfo.end(), evtInfoTSort );

	/*---- Generate particle kinematics and save into rootfile ----*/

	// define class
	MCInfo *mc = new MCInfo;
	mc->Clear();
	mc->SetName("MC");

	SNEvtInfo *sngen = new SNEvtInfo;
	sngen->Clear();
	sngen->SetName("SN");

	// define Branch
	TList *TopBranch = new TList;
	TopBranch->Add(mc);
	TopBranch->Add(sngen);

	// make output root file

	int nsub = 0;

	std::ostringstream sout;
	sout << std::setfill('0') << std::setw(6) << nsub;
	std::stringstream ssname;
	ssname << OutDir << sout.str() << ".root";
	std::string fname = ssname.str();
	std::cout << "file name " << fname << std::endl;

	TFile *fout = new TFile(fname.c_str(), "RECREATE");
	//------------------------------------------------------------------------
	// set write cache to 40MB 
	TFileCacheWrite *cachew = new TFileCacheWrite(fout, 40*1024*1024);
	//------------------------------------------------------------------------

	// define tree
	TTree *theOTree = new TTree("data","SK 3 tree");
	// new MF
	theOTree->SetCacheSize(40*1024*1024);
	int bufsize = 8*1024*1024;      // may be this is the best 15-OCT-2007 Y.T.
	theOTree->Branch(TopBranch,bufsize);

	int totGenNuebarp=0;
	int totGenNueElastic=0, totGenNuebarElastic=0, totGenNuxElastic=0, totGenNuxbarElastic=0;
	int totGenNueO=0, totGenNuebarO=0;
	int totGenNcNup=0, totGenNcNun=0, totGenNcNubarp=0, totGenNcNubarn=0;

	for( uint iEvt = 0; iEvt < vEvtInfo.size(); iEvt++ ){

		iSkip = 0;

		if((iEvt != 0) && (iEvt%NeventFile == 0) ) {

		  fout->cd();
		  theOTree->Write();
		  theOTree->Reset();
		  fout->Close();
		  delete fout;

		  // make output root file

		  nsub++;

		  std::ostringstream sout;
		  sout << std::setfill('0') << std::setw(6) << nsub;
		  std::stringstream ssname;
		  ssname << OutDir << sout.str() << ".root";
		  std::string fname = ssname.str();
		  std::cout << "file name " << fname << std::endl;

		  fout = new TFile(fname.c_str(), "RECREATE");
		  //------------------------------------------------------------------------
		  // set write cache to 40MB 
		  TFileCacheWrite *cachew = new TFileCacheWrite(fout, 40*1024*1024);
		  //------------------------------------------------------------------------

		  theOTree = new TTree("data","SK 3 tree");
		  theOTree->SetCacheSize(40*1024*1024);
		  theOTree->Branch(TopBranch,bufsize);

		}


		SNEvtInfo & p = vEvtInfo[iEvt];

		// fill SNEvtInfo (see $SKOFL_ROOT/include/lowe/snevtinfo.h )

		sngen->iEvt = iEvt;
		sngen->rType = p.rType;
		sngen->rTime = p.rTime;
		sngen->nuType = p.nuType;
		sngen->nuEne = p.nuEne;
		sngen->nuDir[0] = p.nuDir[0];
		sngen->nuDir[1] = p.nuDir[1];
		sngen->nuDir[2] = p.nuDir[2];
		sngen->rVtx[0] = p.rVtx[0];
		sngen->rVtx[1] = p.rVtx[1];
		sngen->rVtx[2] = p.rVtx[2];

		// MCVERTEX (see $SKOFL_ROOT/inc/vcvrtx.h )

		mc->nvtxvc = 1;
		mc->pvtxvc[0][0] = p.rVtx[0];
		mc->pvtxvc[0][1] = p.rVtx[1];
		mc->pvtxvc[0][2] = p.rVtx[2];
		mc->iflvvc[0] = 1;
		mc->iparvc[0] = 0;
		mc->timvvc[0] = 0.; // impossible store here because it is float and no enough precision for SN time, instead of this, fill it into sngen->rTime above

		//std::cout << iEvt << " t=" << p.rTime << " " << p.rType << " " << p.nuType << " E=" << p.nuEne << " x=" << p.rVtx[0] << " y=" << p.rVtx[1] << " z=" << p.rVtx[2] << std::endl;

		// Calculate neutrino interaction vector and save into MCVECT
		determineKinematics( p.rType, p.nuEne, p.nuDir, mc );

		if(iSkip == 0) {
		  theOTree->Fill();

		  if(p.rType == 0) totGenNuebarp++;
		  if(p.rType == 1) totGenNueElastic++;
		  if(p.rType == 2) totGenNuebarElastic++;
		  if(p.rType == 3) totGenNuxElastic++;
		  if(p.rType == 4) totGenNuxbarElastic++;
		}
	}

	int totalNumOfGenEvts = ( totGenNuebarp + totGenNueElastic
				  + totGenNuebarElastic + totGenNuxElastic + totGenNuxbarElastic
				  + totGenNueO + totGenNuebarO 
				  + totGenNcNup + totGenNcNun + totGenNcNubarp + totGenNcNubarn
				  );

	fprintf( stdout, "------------------------------------\n" );
	fprintf( stdout, "total generated number of events %d\n", totalNumOfGenEvts );
	fprintf( stdout, "   nuebar + p = %d\n", totGenNuebarp );
	fprintf( stdout, "   nue + e = %d\n", totGenNueElastic );
	fprintf( stdout, "   nuebar + e = %d\n", totGenNuebarElastic );
	fprintf( stdout, "   nux + e = %d\n", totGenNuxElastic );
	fprintf( stdout, "   nuxbar + e = %d\n", totGenNuxbarElastic );
	fprintf( stdout, "------------------------------------\n" );

	fout->cd();
	theOTree->Write();
	fout->Close();
	delete fout;
}

void VectGenGenerator::MakeEvent(double time, double nu_energy, int nReact, int nuType, double rate){

	double totcrsIBD[nuEneNBins] = {0.};
	double dRandTotEvts = generator->Poisson(rate);

	if(dRandTotEvts > 0){
		for(int i=0; i<dRandTotEvts; i++){
			double ene_s = nu_energy - nuEneBinSize/2., ene_e = nu_energy + nuEneBinSize/2.;
			double nuEne = getRandomReal( ene_s, ene_e, generator );
			double time_s = time - tBinSize/2., time_e = time + tBinSize/2.;
			double tReact = getRandomReal( time_s, time_e , generator );

			double ver_x = 9999., ver_y = 9999., ver_z = 9999.;
			determinePosition( ver_x, ver_y, ver_z );

			SNEvtInfo evtInfo;
			//evtInfo.iEvt = dRandTotEvts;
			evtInfo.rType = nReact;
			evtInfo.rTime = tReact;
			evtInfo.nuType = nuType;
			evtInfo.nuEne = nuEne;
			evtInfo.nuDir[0] = snDir[0];
			evtInfo.nuDir[1] = snDir[1];
			evtInfo.nuDir[2] = snDir[2];
			evtInfo.rVtx[0] = ver_x;
			evtInfo.rVtx[1] = ver_y;
			evtInfo.rVtx[2] = ver_z;
			      
			vEvtInfo.push_back( evtInfo );
		}
	}
	return;
}

void VectGenGenerator::Process(){

	/*---- Fill total cross section into array to avoid repeating calculation ----*/
	double totcrsIBD[nuEneNBins] = {0.};
	double totcrsNue[nuEneNBins] = {0.}, totcrsNueb[nuEneNBins] = {0.}, totcrsNux[nuEneNBins] = {0.}, totcrsNuxb[nuEneNBins] = {0.};

	std::cout << "calculate cross section and fill to array" << std::endl;
	for(int i_nu_ene =0; i_nu_ene < nuEneNBins; i_nu_ene++) {
		double nu_energy = nuEneMin + ( double(i_nu_ene) + 0.5 ) * nuEneBinSize;
		//if(i_nu_ene % 10 == 0) std::cout << "Neutrino Energy  " << nu_energy << std::endl;
		/*----- inverse beta decay -----*/
		if(nu_energy > eEneThr + DeltaM) totcrsIBD[i_nu_ene] = nucrs->CsNuebP_SV(nu_energy);

		/*----- electron elastic -----*/

		nuType = 12; // Nue
		totcrsNue[i_nu_ene] = nucrs->CsNuElastic(nuType, nu_energy, flag_event);

		nuType = -12; // Nuebar
		totcrsNueb[i_nu_ene] = nucrs->CsNuElastic(nuType, nu_energy, flag_event);

		nuType = 14; // Nux
		totcrsNux[i_nu_ene] = nucrs->CsNuElastic(nuType, nu_energy, flag_event);

		nuType = -14; // Nuxbar
		totcrsNuxb[i_nu_ene] = nucrs->CsNuElastic(nuType, nu_energy, flag_event);
	}

	/*---- loop ----*/
	double time;
	double nuEne;
	for(Int_t i_time =0; i_time < tNBins; i_time++) {

		time = tStart + (double(i_time)+0.5)*tBinSize; //center value of each bin[s]
		int itime_sn = int(time);

		if(itime_sn > (int)(tEnd * 1000.)){
			exit(0);
		}
		double rate;
		for(int i_nu_ene =0; i_nu_ene < nuEneNBins; i_nu_ene++) {

			double nu_energy = nuEneMin + ( double(i_nu_ene) + 0.5 ) * nuEneBinSize;

			double nspcne = nuflux->VectGenSnNspeNue(time, nu_energy); //Nue
			double nspcneb = nuflux->VectGenSnNspeNueb(time, nu_energy); //Nuebar
			double nspcnx = nuflux->VectGenSnNspeNux(time, nu_energy); //Nux or Nuxbar

			/*----- inverse beta decay -----*/

			rate = Const_p * (oscneb1*nspcneb + oscneb2*nspcnx) * totcrsIBD[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
			totNuebarp += rate;
			if(flag_event == 1) {
			  nReact = 0;
			  nuType = -12;
			  MakeEvent(time, nu_energy, nReact, nuType, rate);
			}


			/*----- electron elastic -----*/

			rate = Const_e * (oscnue1*nspcne + oscnue2*nspcnx) * totcrsNue[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
			totNueElastic += rate;
			if(flag_event == 1) {
			  nReact = 1;
			  nuType = 12;
			  MakeEvent(time, nu_energy, nReact, nuType, rate);
			}

			rate = Const_e * (oscneb1*nspcneb + oscneb2*nspcnx) * totcrsNueb[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
			totNuebarElastic += rate;
			if(flag_event == 1) {
			  nReact = 2;
			  nuType = -12;
			  MakeEvent(time, nu_energy, nReact, nuType, rate);
			}

			rate = Const_e * (oscnux1*nspcnx + oscnux2*nspcne) * totcrsNux[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
			totNuxElastic += rate;
			if(flag_event == 1) {
			  nReact = 3;
			  nuType = 14;
			  MakeEvent(time, nu_energy, nReact, nuType, rate);
			}

			rate = Const_e * (oscnxb1*nspcnx + oscnxb2*nspcneb) * totcrsNuxb[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
			totNuxbarElastic += rate;
			if(flag_event == 1) {
			  nReact = 4;
			  nuType = -14;
			  MakeEvent(time, nu_energy, nReact, nuType, rate);
			}
		}

		//std::cout << time << " " << totNuebarp << " " << totNueElastic << std::endl;

	}

	double totalNumOfEvts = ( totNuebarp + totNueElastic
			+ totNuebarElastic + totNuxElastic + totNuxbarElastic
			+ totNueO + totNuebarO 
			+ totNcNup + totNcNun + totNcNubarp + totNcNubarn
			);


	fprintf( stdout, "------------------------------------\n" );
	fprintf( stdout, "total expected number of events %e\n", totalNumOfEvts );
	fprintf( stdout, "   nuebar + p = %e\n", totNuebarp );
	fprintf( stdout, "   nue + e = %e\n", totNueElastic );
	fprintf( stdout, "   nuebar + e = %e\n", totNuebarElastic );
	fprintf( stdout, "   nux + e = %e\n", totNuxElastic );
	fprintf( stdout, "   nuxbar + e = %e\n", totNuxbarElastic );
	fprintf( stdout, "------------------------------------\n" );

	/*
	fprintf( stdout, "------------------------------------\n" );
	fprintf( stdout, "totreaction:\n" );
	fprintf( stdout, "   nuebar + p = %e\n", totNuebarp );
	fprintf( stdout, "   nue + e = %e\n", totNueElastic );
	fprintf( stdout, "   nuebar + e = %e\n", totNuebarElastic );
	fprintf( stdout, "   nux + e = %e\n", totNuxElastic );
	fprintf( stdout, "   nuxbar + e = %e\n", totNuxbarElastic );
	fprintf( stdout, "   nue + O = %e\n", totNueO );
	fprintf( stdout, "   nuebar + O = %e\n", totNuebarO );
	fprintf( stdout, "   nu + O = nu + p + N%e\n", totNcNup );
	fprintf( stdout, "   nu + O = nu + n + O%e\n", totNcNun );
	fprintf( stdout, "   nubar + O = nubar + p + N%e\n", totNcNubarp );
	fprintf( stdout, "   nubar + O = nubar + n + O%e\n", totNcNubarn );
	fprintf( stdout, "------------------------------------\n" );
	*/

	if(flag_event == 1) FillEvent();

}

void VectGenGenerator::Process(int NumEv){

	/*-----input file name-----*/
  FluxCalculation &nuflux = *nuflux_dsnb;
//  nuflux.dumpFlux();
//  double nuEne_min = nuEneMin;
//  double nuEne_max = nuEneMax;
//  if( nuEne_min < nuflux.getFluxLimit(true /* true: lower limit, false: higher limit*/) ){
//    nuEne_min = nuflux.getFluxLimit(true);
//    std::cerr << "Flux lower limit is higher than specified value: reset to nu_ene_min = " << nuEne_min << std::endl;
//  }
//  if( nuEne_max > nuflux.getFluxLimit(false) ) {
//    nuEne_max = nuflux.getFluxLimit(false);
//    std::cerr << "Flux upper limit is lower than specified value: reset to nu_ene_max = " << nuEne_max << std::endl;
//  }

//
//	/*-----calculate maximum value-----*/
//	double maxP = 0.;
//	for( int j = 0; j < int(nuflux.getNBins()); j++ ){
//	  double nuEne = nuflux.getBinnedEnergy(j);
//
//	  if((nuEne < nuEne_min) || (nuEne > nuEne_max)) continue;
//
//	  for( int iCost = 0; iCost < costNBins; iCost++ ){
//	    double cost = costMin + costBinSize * ( double(iCost) + 0.5 );
//	    double eEne, sigm;
//	    nucrs->DcsNuebP_SV(nuEne, cost, eEne, sigm );
//
//	    double p = nuflux.getBinnedFlux(j) * sigm;
//	    if( p > maxP ){ maxP = p; }
//	    //std::cout << nuEne << " " << Flux[j] << " " << sigm << " " << p << " " << maxP << std::endl;                          
//	  }
//	}

//	/*-----define class-----*/
//	MCInfo *mc = new MCInfo;
//	mc->Clear();
//	mc->SetName("MC");
//
//	/*-----define branch-----*/
//	TList *TopBranch = new TList;
//	TopBranch->Add(mc);
//
//	/*-----make output root file-----*/
//
//	int nsub = 0;
//
//  auto GenerateOutputFilePath = [](std::string odir, int nsub){
//    std::ostringstream sout;
//    sout << odir << "/" <<  std::setfill('0') << std::setw(6) << nsub << ".root";
//    return sout.str();
//  };
//	std::string fname = GenerateOutputFilePath(OutDir, nsub);
//	std::cout << "file name " << fname << std::endl;
//
//	TFile *fout = new TFile(fname.c_str(), "RECREATE");
//	//------------------------------------------------------------------------
//	// set write cache to 40MB 
//	TFileCacheWrite *cachew = new TFileCacheWrite(fout, 40*1024*1024);
//	//------------------------------------------------------------------------
//
//	// define tree
//	TTree *theOTree = new TTree("data","SK 3 tree");
//	// new MF
//	theOTree->SetCacheSize(40*1024*1024);
//	int bufsize = 8*1024*1024;      // may be this is the best 15-OCT-2007 Y.T.
//	theOTree->Branch(TopBranch,bufsize);
  

	/*---- loop ----*/
	for( int iEvt = 0; iEvt < NumEv; iEvt++ ){
	  std::cout << iEvt << std::endl;

//		if((iEvt != 0) && (iEvt%NeventFile == 0) ) {
//
//		  fout->cd();
//		  theOTree->Write();
//		  theOTree->Reset();
//		  fout->Close();
//		  delete fout;
//
//		  // make output root file
//
//		  nsub++;
//
//		  std::ostringstream sout;
//		  sout << std::setfill('0') << std::setw(6) << nsub;
//		  std::stringstream ssname;
//		  ssname << OutDir << "/" << sout.str() << ".root";
//		  std::string fname = ssname.str();
//		  std::cout << "file name " << fname << std::endl;
//
//		  fout = new TFile(fname.c_str(), "RECREATE");
//		  //------------------------------------------------------------------------
//		  // set write cache to 40MB 
//		  TFileCacheWrite *cachew = new TFileCacheWrite(fout, 40*1024*1024);
//		  //------------------------------------------------------------------------
//
//		  theOTree = new TTree("data","SK 3 tree");
//		  theOTree->SetCacheSize(40*1024*1024);
//		  theOTree->Branch(TopBranch,bufsize);
//
//		}

		// determine neutrino and positron energy, and its direction                                                              
		double nuEne, cost, eEne;
		while( 1 ){
		  nuEne = getRandomReal( nuEne_min, nuEne_max, generator );

      const double nuFlux = nuflux.getFlux(nuEne);

		  double sigm;
		  cost = getRandomReal( costMin, costMax, generator );
		  nucrs->DcsNuebP_SV(nuEne, cost, eEne, sigm );

		  double p = nuFlux * sigm;
		  double x = getRandomReal( 0., maxProb, generator );
		  if( x < p ) break;
		}

		// determine neutrino direction
		double nuDir[3];
		double theta = getRandomReal( 0., 1., generator ) * M_PI;
		double phi = getRandomReal( 0., 1., generator ) * 2. * M_PI;
		nuDir[0] = sin( theta ) * cos( phi );
		nuDir[1] = sin( theta ) * sin( phi );
		nuDir[2] = cos( theta );

		// Rotation matrix of neutrino direction
		double Rmat[3][3];
		Rmat[0][0] = cos( theta ) * cos( phi );
		Rmat[0][1] = -sin( phi );
		Rmat[0][2] = sin( theta ) * cos( phi );
		Rmat[1][0] = cos( theta ) * sin( phi );
		Rmat[1][1] = cos( phi );
		Rmat[1][2] = sin( theta ) * sin( phi );
		Rmat[2][0] = -sin( theta );
		Rmat[2][1] = 0.;
		Rmat[2][2] = cos( theta );

		// interaction point
		double ver_x, ver_y, ver_z;
		determinePosition( ver_x, ver_y, ver_z );

		// Fill into class
		// MCVERTEX (see $SKOFL_ROOT/inc/vcvrtx.h )                                                                               
		fMC->nvtxvc = 1;
		fMC->pvtxvc[0][0] = float(ver_x);
		fMC->pvtxvc[0][1] = float(ver_y);
		fMC->pvtxvc[0][2] = float(ver_z);
		fMC->iflvvc[0] = 1;
		fMC->iparvc[0] = 0;
		fMC->timvvc[0] = 0.;

		// IBD interaction

		fMC->nvc = 4;

		// Original neutrino
		fMC->ipvc[0] = -12; // anti-electron neutrino
		fMC->energy[0] = nuEne; // ENERGY ( MEV )
		fMC->pvc[0][0] = nuEne * nuDir[0]; // MOMENTUM OF I-TH PARTICLE ( MEV/C )
		fMC->pvc[0][1] = nuEne * nuDir[1];
		fMC->pvc[0][2] = nuEne * nuDir[2];
		fMC->iorgvc[0] = 0;  // ID OF ORIGIN PARTICLE  PARENT PARTICLE
		fMC->ivtivc[0] = 0;  // VERTEX # ( INITIAL )
		fMC->iflgvc[0] = -1; // FINAL STATE FLAG
		fMC->icrnvc[0] = 0;  // CHERENKOV FLAG
		fMC->ivtfvc[0] = 1;  // VERTEX # ( FINAL )

		// Original proton
		fMC->ipvc[1] = 2212; // proton
		fMC->energy[1] = Mp; // ENERGY ( MEV )
		fMC->pvc[1][0] = 0.; // MOMENTUM OF I-TH PARTICLE ( MEV/C )
		fMC->pvc[1][1] = 0.;
		fMC->pvc[1][2] = 0.;
		fMC->iorgvc[1] = 0;  // ID OF ORIGIN PARTICLE  PARENT PARTICLE
		fMC->ivtivc[1] = 0;  // VERTEX # ( INITIAL )
		fMC->iflgvc[1] = -1; // FINAL STATE FLAG
		fMC->icrnvc[1] = 0;  // CHERENKOV FLAG
		fMC->ivtfvc[1] = 1;  // VERTEX # ( FINAL )

		// Positron

		double amom = sqrt(SQ( eEne ) - SQ( Me ));
		double eTheta = acos( cost );
		double ePhi = getRandomReal( -M_PI,  M_PI, generator );

		// conversion the positron direction along the neutrino direction
		double origVec[3], eDir[3];

		origVec[0] = sin( eTheta ) * cos( ePhi );
		origVec[1] = sin( eTheta ) * sin( ePhi );
		origVec[2] = cos( eTheta );

		eDir[0] = 0.;
		eDir[1] = 0.;
		eDir[2] = 0.;

		for( int i = 0; i < 3; i++ ){
		  for( int j = 0; j < 3; j++ ){
		    eDir[i] += Rmat[i][j] * origVec[j];
		  }
		}
		//------------------------------------------------------------------------

		// Positron
		fMC->ipvc[2] = -11; // positron
		fMC->energy[2] = eEne; // total ENERGY ( MEV )
		fMC->pvc[2][0] = amom * eDir[0]; // MOMENTUM OF I-TH PARTICLE ( MEV/C )
		fMC->pvc[2][1] = amom * eDir[1];
		fMC->pvc[2][2] = amom * eDir[2];
		fMC->iorgvc[2] = 1;  // ID OF ORIGIN PARTICLE  PARENT PARTICLE
		fMC->ivtivc[2] = 1;  // VERTEX # ( INITIAL )
		fMC->iflgvc[2] = 0; // FINAL STATE FLAG
		fMC->icrnvc[2] = 1;  // CHERENKOV FLAG
		fMC->ivtfvc[2] = 1;  // VERTEX # ( FINAL )

		// Neutron
		fMC->ipvc[3] = 2112; // neutron
		fMC->pvc[3][0] = (fMC->pvc[0][0]) - (fMC->pvc[1][0]); // MOMENTUM OF I-TH PARTICLE ( MEV/C )
		fMC->pvc[3][1] = (fMC->pvc[0][1]) - (fMC->pvc[1][1]);
		fMC->pvc[3][2] = (fMC->pvc[0][2]) - (fMC->pvc[1][2]);
		fMC->energy[3] = sqrt(SQ( fMC->pvc[2][0] ) + SQ( fMC->pvc[2][1] ) + SQ( fMC->pvc[2][2] )  + SQ( Mn )); // ENERGY ( MEV )
		fMC->iorgvc[3] = 1;  // ID OF ORIGIN PARTICLE  PARENT PARTICLE
		fMC->ivtivc[3] = 1;  // VERTEX # ( INITIAL )
		fMC->iflgvc[3] = 0; // FINAL STATE FLAG
		fMC->icrnvc[3] = 1;  // CHERENKOV FLAG
		fMC->ivtfvc[3] = 1;  // VERTEX # ( FINAL )

    fMC->mcinfo[0] = fRefRunNum;

		theOTree->Fill();
	}
}
