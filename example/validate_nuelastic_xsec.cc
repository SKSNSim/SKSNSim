#include <SKSNSimCrosssection.hh>

int main(int argc, char **argv){

    std::unique_ptr<SKSNSimXSecNuElastic> xsec = std::make_unique<SKSNSimXSecNuElastic>();

    constexpr double emin = 1.0;
    constexpr double emax = 100.0;
    constexpr int nstep_per_one = 10.0;
    constexpr int nstep = (emax - emin) * nstep_per_one;

    for( int i = 0; i < nstep; i ++){
        double e = emin + 1.0 / nstep_per_one * i;
        double E = SKSNSimXSecNuElastic::CalcElectronTotEnergy(e, 0.99);
        printf("%10.5f %10.5g\n", e, xsec->GetDiffCrosssection(e, 0.99, PDG_ELECTRON_NEUTRINO).first); 
    }

    return EXIT_SUCCESS;
}
