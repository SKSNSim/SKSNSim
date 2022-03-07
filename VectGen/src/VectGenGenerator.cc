
#include "VectGenGenerator.hh"

void VectGenGenerator::Process(){

	// SN variables in the input root file
	int nEbin = nuEneNBins;
	int nuType;
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

			double nu_energy = nuEneMin + ( i_nu_ene + 0.5 ) * nuEneBinSize;
			double totcrs = 0.;

			double nspcne = nuflux->VectGenSnNspeNue(time, nu_energy); //Nue
			double nspcneb = nuflux->VectGenSnNspeNueb(time, nu_energy); //Nuebar
			double nspcnx = nuflux->VectGenSnNspeNux(time, nu_energy); //Nux or Nuxbar

			totcrs = nucrs->CsNuebP_SV(nu_energy, eEneThr);
			//std::cout<<i_nu_ene<<" "<<"totcrs= "<<totcrs<<std::endl;
			//rate = Const_p * nspcneb * totcrs * nuEneBinSize * tBinSize;
			rate = Const_p * (oscneb1*nspcneb + oscneb2*nspcnx) * totcrs * nuEneBinSize * tBinSize;
			totNuebarp += rate;

			nuType = 12; // Nue
			totcrs = nucrs->CsNuElastic(nuType, nu_energy);
			//rate = Const_e * nspcne * totcrs * nuEneBinSize * tBinSize;
			rate = Const_e * (oscnue1*nspcne + oscnue2*nspcnx) * totcrs * nuEneBinSize * tBinSize;
			totNueElastic += rate;

			nuType = -12; // Nuebar
			totcrs = nucrs->CsNuElastic(nuType, nu_energy);
			//rate = Const_e * nspcneb * totcrs * nuEneBinSize * tBinSize;
			rate = Const_e * (oscneb1*nspcneb + oscneb2*nspcnx) * totcrs * nuEneBinSize * tBinSize;
			totNuebarElastic += rate;

			nuType = 14; // Nux
			totcrs = nucrs->CsNuElastic(nuType, nu_energy);
			//rate = Const_e * nspcnx * totcrs * nuEneBinSize * tBinSize * 2.;
			rate = Const_e * (oscnux1*nspcnx + oscnux2*nspcne) * totcrs * nuEneBinSize * tBinSize * 2.;
			totNuxElastic += rate;

			nuType = -14; // Nuxbar
			totcrs = nucrs->CsNuElastic(nuType, nu_energy);
			//rate = Const_e * nspcnx * totcrs * nuEneBinSize * tBinSize * 2.;
			rate = Const_e * (oscnxb1*nspcnx + oscnxb2*nspcneb) * totcrs * nuEneBinSize * tBinSize * 2.;
			totNuxbarElastic += rate;

		}

		/*
		NuEnergy.push_back( ene );
		intRate00.push_back( rate00 );
		intRate01.push_back( rate01 );
		intRate02.push_back( rate02 );
		intRate03.push_back( rate03 );
		intRate04.push_back( rate04 );
		*/

		//if(i_time % 10000 == 0) std::cout << i_time/10000. << std::endl;
		if(i_time % 100 == 0) std::cout << i_time/1000. << " " << totNuebarp << " " << totNueElastic << std::endl;

	}

	double totalNumOfEvts = ( totNuebarp + totNueElastic
			+ totNuebarElastic + totNuxElastic + totNuxbarElastic
			+ totNueO + totNuebarO 
			+ totNcNup + totNcNun + totNcNubarp + totNcNubarn
			);


	fprintf( stdout, "------------------------------------\n" );
	fprintf( stdout, "total number of events with random %e\n", totalNumOfEvts );
	fprintf( stdout, "totreaction:\n" );
	fprintf( stdout, "   nuebar + p = %e\n", totNuebarp );
	fprintf( stdout, "   nue + e = %e\n", totNueElastic );
	fprintf( stdout, "   nuebar + e = %e\n", totNuebarElastic );
	fprintf( stdout, "   nux + e = %e\n", totNuxElastic );
	fprintf( stdout, "   nuxbar + e = %e\n", totNuxbarElastic );

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

}
