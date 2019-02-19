/**
 *  makeTree.cc
 *     author: Cameron Bravo <bravo@slac.stanford.edu>
 *     Based on baseline.cxx by Omar Moreno
 *     created: October 24, 2018
 */

//--- C++ ---//
//-----------//
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>

//--- DAQ ---//
//-----------//
#include <TrackerEvent.h>
#include <TrackerSample.h>
#include <DataRead.h>
#include <Data.h>

//--- ROOT ---//
//------------//
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH2.h>
#include <TError.h>
#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>

//--- Utils ---//
//-------------//
//#include <//PlotUtils.h>
#include <Apv25Utils.h>

using namespace std; 

void displayUsage(); 

int main(int argc, char **argv)
{

    string input_file_name; 
    string output_path = ".";

    int fpga = -1;
    int hybrid = -1; 
    int flip = 0;

    int half_module_id = -1; 
    int option_char; 
    int number_events = 0; 
    int run = -1;

    // Parse any command line arguments. If there are no valid command line
    // arguments given, print the usage.
    while((option_char = getopt(argc, argv, "i:d:u:n:r:p:m:f:")) != -1){
        switch(option_char){
            case 'i':
                input_file_name = optarg;
                break;
            case 'm': 
                half_module_id = atoi(optarg);
                break;	
            case 'd': 
                hybrid = atoi(optarg); 
                break; 
            case 'n': 
                number_events = atoi(optarg); 
                break;
            case 'r':
                run = atoi(optarg);
                break;
            case 'p':
                output_path = optarg;
                break;
            case 'f':
                flip = atoi(optarg);
                break;
            case 'u':
                displayUsage(); 
            default:
                displayUsage(); 
        }
    }

    // If an input file (binary or EVIO) was not specified, exit the program
    if(input_file_name.length() == 0){
        cout << "\n Please specify an input file to process, without the extension." 
            << "\n Use '-h' flag for usage.\n" << endl;
        return EXIT_FAILURE; 
    }

    // If a valid FPGA address was not specified, exit the program
    /*if(fpga == uint(-1) || fpga < uint(0) || fpga > uint(6)){ 
      cout << "\nPlease specify a valid FPGA address."
      << "\nUse '-h' flag for usage.\n" << endl;
      return EXIT_FAILURE;  
      }*/

    // If a valid hybrid or half-module id was not specified, exit the program
    if(hybrid <= 0 && half_module_id <= 0){
        cerr << "\n\tPlease specify the hybrid or half-module ID of the device under test"
            << "\n\tUse '-h' option to see usage.\n" << endl;
        return EXIT_FAILURE;
    } else if(hybrid > 0 && half_module_id > 0){
        cerr << "\n\tCannot specify both a hybrid and half-module ID!"
            << "\n\tUse '-h' option to see usage.\n" << endl;
        return EXIT_FAILURE;
    }

    // If a run number has not been set, exit the program
    if(run <= 0){ 
        cerr << "\n\tPlease specify a run number."
            << "\n\tUse '-h' option for usage.\n" << endl;
        return EXIT_FAILURE; 
    }

    // 
    TrackerEvent *event = new TrackerEvent(); 
    TrackerSample *samples = NULL;
    DataRead *dataRead = new DataRead(); 

    // Open the input file.  If the input file cannot be opened exit.
    if(!dataRead->open(input_file_name + ".bin")){
        cerr << "\n Error! File " << input_file_name << " cannot be opened. Be sure to NOT include the .bin extension."
            << endl;
        return EXIT_FAILURE;
    }

    cout << "Processing file: " << input_file_name + ".bin" << endl;

    //Setup TFile and TH2Ds
    TFile *oFile = new TFile(Form("%s_HD.root", input_file_name.c_str()),"RECREATE");
    TH2D  smData_hh[6]; 
    for(int ss = 0; ss < 6; ss++) smData_hh[ss] = TH2D( Form("smData%i_hh",ss) , Form("Sample %i Data;Physical Channel # ; ADC Value",ss) , 640, -0.5, 639.5, 16384, -0.5, 16383.5);
    TH2D *smDataSum_hh = new TH2D( "smDataSum_hh" , "Sample Summed Data;Physical Channel # ;ADC Value" , 640, -0.5, 639.5, 16384, -0.5, 16383.5);
    //for(int ss = 0; ss < 6; ss++) smData_hh[ss] = TH2D( Form("smData%i_hh",ss) , Form("Sample %i Data;Physical Channel # ; ADC Value",ss) , 640, -0.5, 639.5, 16384*2, -16383.5, 16383.5);
    //TH2D *smDataSum_hh = new TH2D( "smDataSum_hh" , "Sample Summed Data;Physical Channel # ;ADC Value" , 640, -0.5, 639.5, 16384*2, -16383.5, 16383.5);

    //--- Fill TTree ---//
    //------------------//
    int channel;
    int event_number = 0; 
    double noise, baseline;
    int old_apv; 

    // Loop through all of the events until the end of file is reached 
    while(dataRead->next(event))
    {

        // Skip the event if the event doesn't contain data from the FPGA 
        // of interest
        //if(fpga != event->fpgaAddress()) continue;  

        // Once the number of desired events has been reached, skip the
        // rest of the events
        ++event_number; 
        if(event_number == number_events) break; 
        if(event_number%500 == 0)
        {
            cout << "Event: " << event_number << endl; 
        }

        for(uint x = 0; x < event->count(); x++){

            // Get the samples
            samples = event->sample(x); 

            // Skip the sample if it doesn't come from the hybrid of interest
            //if(hybrid != samples->hybrid()) continue; 

            // Get the physical channel 
            //old_apv = Apv25Utils::getOldApv(samples->apv()); 
            old_apv = samples->apv(); 
            channel = Apv25Utils::getPhysicalChannel(old_apv, samples->channel());

            // For source test, comment out otherwise
            /*if((samples->value(1) > samples->value(0) && samples->value(2) > samples->value(1)) ||
              (samples->value(2) > samples->value(1) && samples->value(3) > samples->value(2)) ||
              (samples->value(3) > samples->value(2) && samples->value(4) > samples->value(3)) ||
              (samples->value(4) > samples->value(3) && samples->value(5) > samples->value(4))) continue; 
              */
            for(int ss = 0; ss < 6; ss++){
                //smData_hh[ss].Fill( channel, double(samples->value(ss)) );
                //smDataSum_hh->Fill( channel, double(samples->value(ss)) );
                //cout << "EN: " << event_number << " modifier: " << (2*(event_number%2)-1) << " sample: " << samples->value(ss) << " result: " << double((2*(event_number%2)-1)*int(samples->value(ss))) << endl;
                smData_hh[ss].Fill( channel, double((2*((event_number+flip)%2)-1)*int(samples->value(ss))) );
                smDataSum_hh->Fill( channel, double((2*((event_number+flip)%2)-1)*int(samples->value(ss))) );
            }
        }
    }

    cout << "Finish processing " << event_number << Form(" events. Saving %s",Form("%s_HD.root", input_file_name.c_str())) << endl;

    // Write the TFile
    oFile->cd();
    smDataSum_hh->Write();
    for(int ss = 0; ss < 6; ss++) smData_hh[ss].Write();
    oFile->Close();

    return EXIT_SUCCESS; 
}

void displayUsage()
{
    cout << "Usage: makeHD [OPTIONS] ..."                 << endl;
    cout << "Example: makeHD -i input_file"           << endl;
    cout << "\n\t -i  Input file to be processed"         << endl;
    cout << "\t -h    Hybrid ID of the device under test"   << endl;
    cout << "\t -r    Run number. Used if ROOT output is enabled" << endl;
    cout << "\t -p    Path to save all output files to"     << endl;
    cout << "\t -m 	  Half-module ID of the device under test" << endl;
    cout << "\t -u    Show this usage \n"                   << endl;
    exit(EXIT_FAILURE);
}
