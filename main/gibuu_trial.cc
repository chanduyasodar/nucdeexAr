#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TText.h"

#include "NucDeExUtils.hh"
#include "NucDeExRandom.hh"
#include "NucDeExDeexcitation.hh"
#include "NucDeExEventInfo.hh"

using namespace std;
namespace fs = std::filesystem;


// ============================================================
// Masses
// ============================================================

const double mass_neutron = 939.566; // MeV
const double mass_proton  = 938.272; // MeV


// ============================================================
// Parse folder name
//
// Example:
// 40Ar_CCQE_e_100000
//
// Returns:
// target   = 40Ar
// process  = CCQE
// flavour  = e
// nevents  = 100000
// ============================================================

bool ParseFolderName(
    const string& folderName,
    string& target,
    string& process,
    string& flavour,
    long long& requestedEvents
)
{
    vector<string> tokens;

    stringstream ss(folderName);
    string token;

    while (getline(ss, token, '_')) {
        tokens.push_back(token);
    }

    if (tokens.size() != 4) {
        cerr << "ERROR: Folder name must have format:\n";
        cerr << "target_process_flavour_events\n";
        cerr << "Example:\n";
        cerr << "40Ar_CCQE_e.n_100000\n";
        return false;
    }

    target = tokens[0];
    process = tokens[1];
    flavour = tokens[2];

    try {
        requestedEvents = stoll(tokens[3]);
    }
    catch (...) {
        cerr << "ERROR: Cannot read number of events from folder name.\n";
        return false;
    }

    return true;
}


// ============================================================
// Main
// ============================================================

int main(int argc, char* argv[])
{

    if (argc != 2) {

        cerr << "\nUsage:\n";

        cerr << "  " << argv[0]
             << " <GiBUU-folder-name>\n\n";

        cerr << "Example:\n";

        cerr << "  " << argv[0]
             << " 40Ar_CCQE_e.n_100000\n\n";

        return 1;
    }


    // ========================================================
    // Input
    // ========================================================

    string folderName = argv[1];

    string target;
    string process;
    string flavour;

    long long requestedEvents = 0;


    if (!ParseFolderName(
            folderName,
            target,
            process,
            flavour,
            requestedEvents
        ))
    {
        return 1;
    }


    cout << "\n========================================\n";

    cout << "GiBUU analysis\n";

    cout << "========================================\n";

    cout << "Target           : " << target << endl;

    cout << "Process          : " << process << endl;

    cout << "Flavour          : " << flavour << endl;

    cout << "Requested events : "
         << requestedEvents << endl;


    // ========================================================
    // Folder paths
    // ========================================================

    string inputDir =
        "output_gibuu/" + folderName;


    cout << "\nInput directory:\n";

    cout << inputDir << endl;


    if (!fs::exists(inputDir)) {

        cerr << "\nERROR: Directory does not exist:\n";

        cerr << inputDir << endl;

        return 1;
    }


    // ========================================================
    // Find ROOT files
    // ========================================================

    vector<string> rootFiles;


    for (const auto& entry :
         fs::directory_iterator(inputDir))
    {

        if (!entry.is_regular_file())
            continue;


        string filename =
            entry.path().string();


        if (entry.path().extension() == ".root") {

            rootFiles.push_back(filename);
        }
    }


    sort(
        rootFiles.begin(),
        rootFiles.end()
    );


    if (rootFiles.empty()) {

        cerr << "\nERROR: No ROOT files found.\n";

        return 1;
    }


    cout << "\nROOT files found:\n";


    for (const auto& file : rootFiles) {

        cout << "  "
             << file
             << endl;
    }


    // ========================================================
    // Create TChain
    //
    // All ROOT files contain RootTuple
    // ========================================================

    TChain chain("RootTuple");


    for (const auto& file : rootFiles) {

        chain.Add(file.c_str());
    }


    long long totalEntries =
        chain.GetEntries();


    cout << "\nTotal entries available:\n";

    cout << totalEntries << endl;


    // ========================================================
    // NucDeEx initialization
    // ========================================================

    int seed = 1;

    int verbose = 1;


    cout << "\nVERBOSE = "
         << verbose
         << endl;


    cout << "SEED = "
         << seed
         << endl;


    NucDeEx::Utils::fVerbose = verbose;


    NucDeEx::Random::SetSeed(seed);


    NucDeExDeexcitation* deex =
        new NucDeExDeexcitation();


    deex->Init();


    // ========================================================
    // Target information
    // ========================================================

    int Zt = 0;

    int Nt = 0;

    double S = 0.0;

    const char* nuc_name = nullptr;


    if (target == "12C") {

        Zt = 6;

        Nt = 6;

        nuc_name = "12C";
    }

    else if (target == "16O") {

        Zt = 8;

        Nt = 8;

        nuc_name = "16O";
    }

    else if (target == "40Ar") {

        Zt = 18;

        Nt = 22;

        nuc_name = "40Ar";
    }

    else {

        cerr << "\nERROR: Unsupported target:\n";

        cerr << target << endl;

        return 1;
    }


    cout << "\nTarget nucleus:\n";

    cout << "Z = "<< Zt<< endl;

    cout << "N = "<< Nt<< endl;


    // ========================================================
    // ROOT branches
    // ========================================================

    vector<int>* barcode = nullptr;


    vector<double>* Px = nullptr;

    vector<double>* Py = nullptr;

    vector<double>* Pz = nullptr;

    vector<double>* E = nullptr;


    double weight = 0.0;


    int process_ID = 0;

    int flavor_ID = 0;

    int evType = 0;


    double lepIn_E  = 0.0;

    double lepIn_Px = 0.0;

    double lepIn_Py = 0.0;

    double lepIn_Pz = 0.0;


    double lepOut_E  = 0.0;

    double lepOut_Px = 0.0;

    double lepOut_Py = 0.0;

    double lepOut_Pz = 0.0;


    double nuc_E  = 0.0;

    double nuc_Px = 0.0;

    double nuc_Py = 0.0;

    double nuc_Pz = 0.0;


    int nuc_charge = 0;


    // ========================================================
    // Connect branches
    // ========================================================

    chain.SetBranchAddress(
        "weight",
        &weight
    );


    chain.SetBranchAddress(
        "barcode",
        &barcode
    );


    chain.SetBranchAddress(
        "Px",
        &Px
    );


    chain.SetBranchAddress(
        "Py",
        &Py
    );


    chain.SetBranchAddress(
        "Pz",
        &Pz
    );


    chain.SetBranchAddress(
        "E",
        &E
    );


    chain.SetBranchAddress(
        "process_ID",
        &process_ID
    );


    chain.SetBranchAddress(
        "flavor_ID",
        &flavor_ID
    );


    chain.SetBranchAddress(
        "evType",
        &evType
    );


    // Incoming neutrino

    chain.SetBranchAddress(
        "lepIn_E",
        &lepIn_E
    );


    chain.SetBranchAddress(
        "lepIn_Px",
        &lepIn_Px
    );


    chain.SetBranchAddress(
        "lepIn_Py",
        &lepIn_Py
    );


    chain.SetBranchAddress(
        "lepIn_Pz",
        &lepIn_Pz
    );


    // Outgoing lepton

    chain.SetBranchAddress(
        "lepOut_E",
        &lepOut_E
    );


    chain.SetBranchAddress(
        "lepOut_Px",
        &lepOut_Px
    );


    chain.SetBranchAddress(
        "lepOut_Py",
        &lepOut_Py
    );


    chain.SetBranchAddress(
        "lepOut_Pz",
        &lepOut_Pz
    );


    // Initial struck nucleon

    chain.SetBranchAddress(
        "nuc_E",
        &nuc_E
    );


    chain.SetBranchAddress(
        "nuc_Px",
        &nuc_Px
    );


    chain.SetBranchAddress(
        "nuc_Py",
        &nuc_Py
    );


    chain.SetBranchAddress(
        "nuc_Pz",
        &nuc_Pz
    );


    chain.SetBranchAddress(
        "nuc_charge",
        &nuc_charge
    );


    // ========================================================
    // Histograms
    // ========================================================


    TH1D* h_Pinit =
        new TH1D(
            "h_Pinit",
            "",
            100,
            0,
            500
        );


    TH1D* h_nmulti_postFSI =
        new TH1D(
            "h_nmulti_postFSI",
            "",
            10,
            -0.5,
            9.5
        );


    TH1D* h_nmulti_postdeex =
        new TH1D(
            "h_nmulti_postdeex",
            "",
            10,
            -0.5,
            9.5
        );


    TH1D* h_MissE =
        new TH1D(
            "h_MissE",
            "",
            400,
            0,
            200
        );


    TH1D* h_Ex[4];


    for (int i = 0; i < 4; i++) {

        stringstream os;

        os << "h_Ex_" << i;


        h_Ex[i] =
            new TH1D(
                os.str().c_str(),
                "",
                500,
                -100,
                400
            );
    }


    TH1D* h_Ex_multi =
        new TH1D(
            "h_Ex_multi",
            "",
            500,
            -100,
            400
        );


    TH2D* h_MissE_Pinit =
        new TH2D(
            "h_MissE_Pinit",
            "",
            100,
            0,
            500,
            400,
            0,
            200
        );


    // ========================================================
    // Number of events to process
    // ========================================================

    long long nToProcess =
        min(
            requestedEvents,
            totalEntries
        );


    cout << "\nProcessing events:\n";

    cout << nToProcess
         << endl;


    // ========================================================
    // Event loop
    // ========================================================

    long long processedEvents = 0;


    for (
        long long i = 0;
        i < totalEntries;
        i++
    )
    {

        if (processedEvents >= nToProcess)
            break;


        chain.GetEntry(i);


        // ----------------------------------------------------
        // IMPORTANT:
        //
        // Each ROOT entry may contain multiple instances
        // because barcode, Px, Py, Pz and E are vectors.
        //
        // We only use Instance 0.
        // ----------------------------------------------------


        if (
            !barcode ||
            !Px ||
            !Py ||
            !Pz ||
            !E
        )
        {
            continue;
        }


        if (
            barcode->empty() ||
            Px->empty() ||
            Py->empty() ||
            Pz->empty() ||
            E->empty()
        )
        {
            continue;
        }


        // ====================================================
        // Instance 0
        // ====================================================

        int pdg =
            barcode->at(0);


        double px =
            Px->at(0);


        double py =
            Py->at(0);


        double pz =
            Pz->at(0);


        double energy =
            E->at(0);


        // ====================================================
        // Incoming neutrino energy
        // ====================================================

        double Ev =
            lepIn_E;


        // ====================================================
        // Outgoing lepton energy
        // ====================================================

        double El =
            lepOut_E;


        // ====================================================
        // Initial struck nucleon momentum
        //
        // nuc_Px, nuc_Py, nuc_Pz
        // ====================================================

        double Pinit =
            sqrt(
                nuc_Px * nuc_Px +
                nuc_Py * nuc_Py +
                nuc_Pz * nuc_Pz
            ) * 1e3;


        h_Pinit->Fill(
            Pinit
        );


        // ====================================================
        // Initial struck nucleon energy
        //
        // nuc_E
        //
        // GiBUU unit: GeV
        // ====================================================

        double Enucc =nuc_E;


        // ====================================================
        // Determine nucleon mass
        //
        // nuc_charge = 0 -> neutron
        // nuc_charge = 1 -> proton
        // ====================================================

        double massnuc = 0.0;


        if (nuc_charge == 0) {

            massnuc =
                mass_neutron * 1e-3;
        }

        else if (nuc_charge == 1) {

            massnuc =
                mass_proton * 1e-3;
        }

        else {

            cerr << "\nWARNING: Unknown nucleon charge: "
                 << nuc_charge
                 << endl;

            continue;
        }


        // ====================================================
        // Missing energy
        //
        // Emiss =
        //
        // Ev - El - Enucc + m_N
        //
        // All GiBUU quantities are GeV.
        //
        // Convert result to MeV.
        // ====================================================

        double MissE =
            (Ev - El - Enucc + massnuc  )*1e3;


        h_MissE->Fill(MissE);


        h_MissE_Pinit->Fill(Pinit,MissE);


        // ====================================================
        // Separation energy
        //
        // Same logic as genie.cc
        // ====================================================

        int Z = Zt;

        int N = Nt;


        if (process == "CCQE") {


            // ------------------------------------------------
            // Neutrino
            //
            // neutrino + n -> l- + p
            //
            // Initial struck nucleon is neutron
            // ------------------------------------------------

            if (
                nuc_charge == 0
            )
            {

                Z++;

                N--;


                S =
                    NucDeEx::Utils::
                    NucleusTable->GetNucleusPtr(nuc_name)->S[1];
}


            // ------------------------------------------------
            // Antineutrino
            //
            // antineutrino + p -> l+ + n
            //
            // Initial struck nucleon is proton
            // ------------------------------------------------

            else if (
                nuc_charge == 1
            )
            {

                Z--;

                N++;


                S =
                    NucDeEx::Utils::
                    NucleusTable
                    ->
                    GetNucleusPtr(
                        nuc_name
                    )
                    ->
                    S[2];
            }
        }


        // ====================================================
        // Excitation energy
        //
        // Ex = Emiss - Separation energy
        // ====================================================

        double Ex =MissE-S;


        h_Ex[0]->Fill( Ex);


        // ====================================================
        // Neutron multiplicity
        //
        // Count neutrons in final GiBUU particle list.
        //
        // barcode is PDG code.
        //
        // PDG:
        // neutron = 2112
        //
        // This is the final-state particle multiplicity
        // available in RootTuple.
        // ====================================================

        int nmulti = 0;


        for (
            size_t k = 0;
            k < barcode->size();
            k++
        )
        {

            if (
                barcode->at(k)
                ==
                2112
            )
            {

                nmulti++;
            }
        }


        // ----------------------------------------------------
        // Keep GENIE-compatible histogram name
        // ----------------------------------------------------

        h_nmulti_postFSI->Fill(
            nmulti
        );


        // ====================================================
        // NucDeEx simulation
        // ====================================================

        NucDeExEventInfo result =
            deex->DoDeex(
                Zt,
                Nt,
                Z,
                N,
                Ex,
                TVector3(
                    0,
                    0,
                    0
                )
            );


        // ====================================================
        // Shell information
        // ====================================================

        int shell =
            result.fShell;


        // ====================================================
        // Add neutrons produced by de-excitation
        // ====================================================

        vector<NucDeExParticle> particle =
            result.ParticleVector;


        int size =
            particle.size();


        for (
            int k = 0;
            k < size;
            k++
        )
        {

            NucDeExParticle p =
                particle.at(k);


            if (
                p._PDG
                ==
                2112
            )
            {

                nmulti++;
            }
        }


        // ====================================================
        // Fill shell excitation histogram
        // ====================================================

        if (
            shell >= 0
            &&
            shell < 4
        )
        {

            h_Ex[shell]->Fill(
                Ex
            );
        }


        // ====================================================
        // Neutron multiplicity after de-excitation
        // ====================================================

        h_nmulti_postdeex->Fill(
            nmulti
        );


        processedEvents++;


        // ====================================================
        // Progress
        // ====================================================

        if (
            processedEvents % 10000
            ==
            0
        )
        {

            cout
                << "Processed "
                << processedEvents
                << " / "
                << nToProcess
                << endl;
        }
    }


    // ========================================================
    // Final event count
    // ========================================================

    cout << "\n========================================\n";

    cout << "Processed events:\n";

    cout << processedEvents
         << endl;

    cout << "========================================\n";


    // ========================================================
    // Output directories
    // ========================================================

    fs::create_directories(
        "fig_gibuu"
    );


    fs::create_directories(
        "output_gibuu"
    );


    // ========================================================
    // ROOT style
    // ========================================================

    gStyle->SetTextFont(
        132
    );


    gStyle->SetTextSize(
        0.08
    );


    gStyle->SetTitleSize(
        0.05,
        "XYZ"
    );


    gStyle->SetTitleFont(
        132,
        "XYZ"
    );


    gStyle->SetLabelSize(
        0.05,
        "XYZ"
    );


    gStyle->SetLabelFont(
        132,
        "XYZ"
    );


    gStyle->SetLegendFont(
        132
    );


    gStyle->SetTitleYOffset(
        0.95
    );


    // ========================================================
    // Plot:
    // Initial nucleon momentum
    // ========================================================

    TCanvas* c_Pinit =
        new TCanvas(
            "c_Pinit",
            "c_Pinit",
            0,
            0,
            800,
            600
        );


    h_Pinit
        ->
        GetXaxis()
        ->
        SetTitle(
            "Momentum of target nucleon (MeV)"
        );


    h_Pinit
        ->
        GetYaxis()
        ->
        SetTitle(
            "Events/bin"
        );


    h_Pinit
        ->
        Draw(
            "HIST"
        );


    string outputName =
        "fig_gibuu/fig_Pinit_"
        +
        folderName
        +
        ".pdf";


    c_Pinit
        ->
        Print(
            outputName.c_str()
        );


    // ========================================================
    // Plot:
    // Neutron multiplicity
    // ========================================================

    TCanvas* c_nmulti =
        new TCanvas(
            "c_nmulti",
            "c_nmulti",
            0,
            0,
            800,
            600
        );


    h_nmulti_postFSI
        ->
        GetXaxis()
        ->
        SetTitle(
            "Neutron multiplicity"
        );


    h_nmulti_postFSI
        ->
        GetYaxis()
        ->
        SetTitle(
            "Fraction of events"
        );


    h_nmulti_postFSI
        ->
        SetStats(
            0
        );


    // --------------------------------------------------------
    // Normalize
    // --------------------------------------------------------

    if (
        h_nmulti_postFSI
        ->
        GetEntries()
        >
        0
    )
    {

        h_nmulti_postFSI
            ->
            Scale(
                1.0
                /
                h_nmulti_postFSI
                ->
                GetEntries()
            );
    }


    if (
h_nmulti_postdeex->GetEntries()> 0)
    { h_nmulti_postdeex->Scale( 1.0/h_nmulti_postdeex->GetEntries());
    }


    // --------------------------------------------------------
    // Draw GiBUU neutrons
    // --------------------------------------------------------

    h_nmulti_postFSI->Draw("HIST");


    // --------------------------------------------------------
    // Draw after NucDeEx
    // --------------------------------------------------------

    h_nmulti_postdeex->SetLineColor(kRed);
    h_nmulti_postdeex->Draw("HIST SAME");


    // --------------------------------------------------------
    // Mean multiplicity
    // --------------------------------------------------------

    stringstream meanText;


    meanText
        <<
        "Mean n multi. = "
        <<
        fixed
        <<
        setprecision(3)
        <<
        h_nmulti_postFSI
        ->
        GetMean();


    TLatex* l_mean =
        new TLatex(
            3,
            h_nmulti_postFSI
            ->
            GetMaximum()
            *
            0.8,
            meanText.str().c_str()
        );


    l_mean
        ->
        Draw(
            "same"
        );


    meanText.str(
        ""
    );


    meanText
        <<
        "Mean n multi. = "
        <<
        fixed
        <<
        setprecision(3)
        <<
        h_nmulti_postdeex
        ->
        GetMean();


    TLatex* l_mean_deex =
        new TLatex(
            3,
            h_nmulti_postFSI
            ->
            GetMaximum()
            *
            0.7,
            meanText.str().c_str()
        );


    l_mean_deex
        ->
        SetTextColor(
            kRed
        );


    l_mean_deex
        ->
        Draw(
            "same"
        );


    outputName =
        "fig_gibuu/fig_nmulti_"
        +
        folderName
        +
        ".pdf";


    c_nmulti
        ->
        Print(
            outputName.c_str()
        );


    // ========================================================
    // Plot:
    // Missing energy
    // ========================================================

    TCanvas* c_MissE =
        new TCanvas(
            "c_MissE",
            "c_MissE",
            0,
            0,
            800,
            600
        );


    h_MissE
        ->
        GetXaxis()
        ->
        SetTitle(
            "Missing energy (MeV)"
        );


    h_MissE
        ->
        GetYaxis()
        ->
        SetTitle(
            "Events/bin"
        );


    h_MissE
        ->
        Draw(
            "HIST"
        );


    outputName =
        "fig_gibuu/fig_MissE_"
        +
        folderName
        +
        ".pdf";


    c_MissE
        ->
        Print(
            outputName.c_str()
        );


    // ========================================================
    // Plot:
    // Excitation energy
    // ========================================================

    TCanvas* c_Ex =
        new TCanvas(
            "c_Ex",
            "c_Ex",
            0,
            0,
            800,
            600
        );


    h_Ex[0]
        ->
        GetXaxis()
        ->
        SetTitle(
            "Excitation energy (MeV)"
        );


    h_Ex[0]
        ->
        GetYaxis()
        ->
        SetTitle(
            "Events/bin"
        );


    h_Ex[0]
        ->
        SetStats(
            0
        );


    h_Ex[0]
        ->
        GetXaxis()
        ->
        SetRangeUser(
            -10,
            100
        );


    h_Ex[0]
        ->
        Draw(
            "HIST"
        );


    outputName =
        "fig_gibuu/fig_Ex_"
        +
        folderName
        +
        ".pdf";


    c_Ex
        ->
        Print(
            outputName.c_str()
        );


    // ========================================================
    // Plot:
    // Missing energy vs initial momentum
    // ========================================================

    TCanvas* c_MissE_Pinit =
        new TCanvas(
            "c_MissE_Pinit",
            "c_MissE_Pinit",
            0,
            0,
            800,
            600
        );


    gPad
        ->
        SetLogz();


    h_MissE_Pinit
        ->
        GetXaxis()
        ->
        SetTitle(
            "Momentum of target nucleon (MeV)"
        );


    h_MissE_Pinit
        ->
        GetYaxis()
        ->
        SetTitle(
            "Missing energy (MeV)"
        );


    h_MissE_Pinit
        ->
        GetYaxis()
        ->
        SetRangeUser(
            0,
            100
        );


    h_MissE_Pinit
        ->
        SetStats(
            0
        );


    h_MissE_Pinit
        ->
        Draw(
            "COLZ"
        );


    outputName =
        "fig_gibuu/fig_MissE_Pinit_"
        +
        folderName
        +
        ".pdf";


    c_MissE_Pinit
        ->
        Print(
            outputName.c_str()
        );


    // ========================================================
    // Save ROOT histograms
    // ========================================================

    string rootOutput =
        "output_gibuu/histogram_deex_"
        +
        folderName
        +
        ".root";


    TFile* outf =
        new TFile(
            rootOutput.c_str(),
            "RECREATE"
        );


    h_Pinit
        ->
        Write();


    h_nmulti_postFSI
        ->
        Write();


    h_nmulti_postdeex
        ->
        Write();


    h_MissE
        ->
        Write();


    for (
        int i = 0;
        i < 4;
        i++
    )
    {

        h_Ex[i]
            ->
            Write();
    }


    h_Ex_multi
        ->
        Write();


    h_MissE_Pinit
        ->
        Write();


    outf
        ->
        Close();


    delete outf;


    delete deex;


    cout << "\n========================================\n";

    cout << "Output ROOT file:\n";

    cout << rootOutput << endl;

    cout << "\nDone.\n";

    cout << "========================================\n";


    return 0;
}
