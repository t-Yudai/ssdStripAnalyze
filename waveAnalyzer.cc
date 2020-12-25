#include <iostream>
#include <vector>
#include <fstream>

#include "Math/Math.h"
#include "Math/SpecFuncMathCore.h"

//const string input_filename = "../../tmpdatalink/tmprootfile/run00385.root.txt";
//const string input_filename = "../../tmpdatalink/tmprootfile/run00918.root.txt";
const string input_filename = "../../tmpdatalink/tmprootfile/run00956.root.txt";
const string input_param_filename = "../paramTable/fitParam.dat";
//TString input_rootfilename = "../../tmpdatalink/tmprootfile/run00385.root";
//TString input_rootfilename = "../../tmpdatalink/tmprootfile/run00918.root";
TString input_rootfilename = "../../tmpdatalink/tmprootfile/run00956.root";
//TString output_rootfilename = "../rootfile_analyzed/run00385_analyzed.root";
TString output_rootfilename = "../rootfile_analyzed/run00956_analyzed_test.root";
TString output_pdf_name_a4 = "../pdf/hoge_a4_test.pdf";
TString output_pdf_name_b5 = "../pdf/hoge_b5_test.pdf";
const string input_filename_calib[2]
= {"../../tmpdatalink/Calibdata/SsdCalibration/SSD1/MeanRms_202006152204.txt",
   "../../tmpdatalink/Calibdata/SsdCalibration/SSD2/MeanRms_202006152204.txt"};



const int n_sample = 8;
const int n_sample_spline = 200;
const int n_plane = 6;
const int n_chip = 6;
const int n_strip = 128;
const int n_stripOnPlane = n_chip * n_strip;
const int PLANE = 2;
const int CHIP = 1;
const int n_parameter = 3;
const int maxEventInEntry = 100;



struct _waveEvent{
   int cluster_num; //clusterize
   int entry;
   int plane;
   int chip;
   int strip;
   int tdc;
   int flag;
   double waveHeight;
   double adc_max;
   int sample_max;
   int sample_min;
   double fit_range_min;
   double fit_range_max;
   double fitted_peak_time;
   double fitted_peak_adc;
   double fitted_hitTiming;
   double fitted_peak_adc_Error;
   double fitted_hitTiming_Error;
   double timing_weighted_mean;
   double chiSquare;
   int    ndf;
   double chiSquare_over_ndf;
   double a_fitParam[n_parameter];
   double a_fitParam_Error[n_parameter];
   vector<double> v_waveform;
   vector<double> v_timing;
   vector<double> v_waveform_original;
   vector<double> v_timing_original;
   double a_waveform[n_sample];
   double a_timing[n_sample];
   double a_waveform_original[n_sample];
   double a_timing_original[n_sample];
   vector<double> v_waveform_spline;
   vector<double> v_timing_spline;
   double a_waveform_spline[n_sample_spline];
   double a_timing_spline[n_sample_spline];
   TGraph *graph_wave;
   TSpline3 *spline;
   TGraph *graph_spline;
   TF1 *fit_func_chip;
   TGraphErrors *graph_wave_Error;

   int waveType;

   _waveEvent(void) : flag(1) {}
};


struct _cluster{
   vector<_waveEvent> cl_waveEvent;
   int w_num;
   int entry_cluster;
   int c_cluster_num;
   double strip_cluster;//adc -> a_fitParam_Error[0]/e
   double simple_strip_cluster;//adc -> adc_max
   double adc_max_cluster;
   double tdc_cluster;
   double chiSquareSum;
   double ndfSum;
};


struct _slice{
   TH1D *h_projection;
   vector<double> v_histdata;
   double mean;
   double stdDev;
};


struct _wave{
   int entry;
   int tdc;
   std::vector<double> v_waveform;
   std::vector<double> v_timing;
};

struct _loop{
   vector<_slice> c_slice;
   TH2D *h_waves_piled;
   TF1  *fit_func;
   Double_t fit_parameter[n_parameter];

};


struct _strip{
   vector<double> v_calib_mean;
   vector<double> v_calib_stdDev;
   vector<double> v_timing_resolution;
   double a_calib_mean[n_sample];
   double a_calib_stdDev[n_sample];
   double a_timing_resolution[n_sample];
};

struct _chip{
   vector<_waveEvent> c_waveEvent;
   vector<_cluster> v_cluster;
   vector<_loop>      c_loop;
   double a_fitParam[n_parameter];
   TH1D *h_chi_square;
   TH1D *h_peak_time;
   double m_calib_mean[n_strip][n_sample];
   double m_calib_stdDev[n_strip][n_sample];
   _strip a_strip_calibdata[n_strip];
};


struct _allchip{
   TH1D *h_chi_square_over_ndf;
   TH1D *h_chi_square;
   TF1  *chi_square_fit;
   TH1D *h_peak_time;
   TH1D *h_timeZero;
   TH1D *h_timeZero_Max3;
   TH1D *h_timeZero_Max4;
   TH1D *h_timeZero_Max5;
   THStack *hs_timeZero;
   TH2D *h_chiSquare_vs_waveHeight;
   TH1D *h_waveHeight;
   int n_wave;
   _allchip(void) : n_wave(0) {}
};


struct _stripID{
   int plane;
   int chip;
   int strip;
   int orderInChip;
};

struct _eventHitID{
   vector< _stripID > v_stripID;
};

struct _allEventID{
   vector< _eventHitID > v_eventHitID;
};

struct _plane{
   TGraphErrors *hit_map;
   vector<double> v_stripNo;
   vector<double> v_stripNo_Error;
   vector<double> v_numOfHit;
   vector<double> v_numOfHit_Error;
   double  hitSumCount;
   double  hitAverage;
   double  a_stripNo[n_stripOnPlane];
   double  a_stripNo_Error[n_stripOnPlane];
   double  a_numOfHit[n_stripOnPlane];
   double  a_numOfHit_Error[n_stripOnPlane];

};

struct _planes{
   TGraphErrors *hit_map_planes;
};

struct _entry{
   _chip m_chip[n_plane][n_chip];
};


bool isMaxAt3(vector<double>);
bool isMaxAt4(vector<double>);
bool isMaxAt5(vector<double>);
bool isBiggerAt8Than0(vector<double>);
bool isTailLow(vector<double>);
bool isNotNoisyStrip(int, int ,int);
double  getMax(vector<double>);
double  getMin(vector<double>);
double  getMean(vector<double>);
double  getStdDev(vector<double>);
double  getWaveHeight(vector<double>);
int     getMaxSampleNum(vector<double>);
double  getTimingWeightedMean(vector<double>, vector<double> ,int);
double  getMinSampleNum(vector<double>);
bool waveThreshold(vector<double>, _chip m_chip[n_plane][n_chip], int, int, int);
int     classifyWaveType(vector<double>);
TVector2 calcPosDecalt(double, double, double);

void waveAnalyzer_2999(){

   vector<_entry> v_entry;
   _allchip s_allchip;
   _allEventID s_allEventID;
   _chip m_chip[6][6];
   ifstream ifs;
   ifs.open(input_filename);
   if(ifs.fail()){
      cerr << "cannot open .txt file" << endl;
      return 0;
   }
   vector< _wave > c_wave;
   int tdc;

   ifstream ifs_calib[2];
   for(int data=0; data<2; data++){
      ifs_calib[data].open(input_filename_calib[data]);
      if(ifs_calib[data].fail()){
         cerr << "cannot open calibration data" << endl;
         return 0;
      }
   }

   TTree *tree;
   TFile *file;
   file = new TFile(input_rootfilename, "read");
   tree = (TTree*)file->Get("ssd_rawdata");
   tree->SetBranchAddress("tdc", &tdc);
   //int n_entry = tree->GetEntries();
   int n_entry = 1000;
   s_allEventID.v_eventHitID.resize(n_entry);
   vector<int> v_signal_info;
   v_signal_info.resize(3);
   vector< vector<int> > m_signal_info;
   int entry_input;
   int entry_multiplicity;
   int damp;
   int numOfWave = 0;
   vector<double> v_wave;
   v_wave.resize(n_sample);
   std::string s;

   //----------

   for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
         for(int strip=0; strip<n_strip; strip++){
            _strip s_strip;
            for(int sample=0; sample<n_sample; sample++){
               s_strip.v_calib_mean.resize(n_sample);
               s_strip.v_calib_stdDev.resize(n_sample);
               s_strip.v_timing_resolution.resize(n_sample);
               s_strip.v_timing_resolution.at(sample) = 8.0 / sqrt(12);
               s_strip.a_timing_resolution[sample] = 8.0 / sqrt(12);
               double mean;
               double stdDev;
               if(plane<3){
                  ifs_calib[0] >> mean;
                  s_strip.v_calib_mean.at(sample) = mean;
                  s_strip.a_calib_mean[sample] = mean;
                  ifs_calib[0] >> stdDev;
                  s_strip.v_calib_stdDev.at(sample) = stdDev;
                  s_strip.a_calib_stdDev[sample] = stdDev;
               }//if plane
               else{
                  ifs_calib[1] >>  mean;
                  s_strip.v_calib_mean.at(sample) = mean;
                  s_strip.a_calib_mean[sample] = mean;
                  ifs_calib[1] >>  stdDev;
                  s_strip.v_calib_stdDev.at(sample) = stdDev;
                  s_strip.a_calib_stdDev[sample] = stdDev;
               }//if plane
               m_chip[plane][chip].a_strip_calibdata[strip] = s_strip;
            }//for sample
         }//for strip
      }//for chip
   }//for plane

   //-----------

   for(int entry=0; entry<n_entry; entry++){
      int plane, chip, strip;
      _eventHitID s_eventHitID;
      tree->GetEntry(entry);
      _wave s_wave;

      while(true){
         ifs >> s;
         if(s[0] == '='){
            break;
         }
      }

      ifs >> entry_input;
      if(entry_input % 1000 == 0){
         cerr << entry_input << "/" << n_entry << endl;
      }
      if(entry_input != entry){
         cerr << "entry_input is not consistent with entry, entry_input = " << entry_input << ", entry = " << entry << endl;
         return 0;
      }

      while(true){
         ifs >> s;
         if(s[0] == '='){
            break;
         }
      }
      ifs >> entry_multiplicity;

      for(int multi=0; multi<entry_multiplicity; multi++){
         ifs >> plane >> chip >> strip;
         for(int sample=0; sample<n_sample; sample++){
            int num;
            ifs >> num;
            v_wave.at(sample) = (double)num;
         }//for sample
         double min = getMin(v_wave);
         double max = getMax(v_wave);
         double range = max - min;
         if( (isMaxAt5(v_wave) || isMaxAt4(v_wave) || isMaxAt3(v_wave))
               && isBiggerAt8Than0(v_wave)
               //&& isTailLow(v_wave)
               && isNotNoisyStrip(plane, chip, strip)
               && waveThreshold(v_wave, m_chip, plane, chip, strip)
           ){
            s_allchip.n_wave++;
            _stripID s_stripID;
            s_stripID.plane = plane;
            s_stripID.chip = chip;
            s_stripID.strip = strip;
            s_stripID.orderInChip = m_chip[plane][chip].c_waveEvent.size();
            s_eventHitID.v_stripID.push_back(s_stripID);
            _waveEvent s_waveEvent;
            s_waveEvent.v_waveform.resize(n_sample);
            s_waveEvent.v_timing.resize(n_sample);
            s_waveEvent.v_waveform_original.resize(n_sample);
            s_waveEvent.v_timing_original.resize(n_sample);
            for(int sample=0; sample<n_sample; sample++){
               s_waveEvent.v_waveform_original.at(sample) = v_wave.at(sample);
               s_waveEvent.a_waveform_original[sample] = v_wave.at(sample);
               s_waveEvent.v_timing_original.at(sample) = sample*25.0;
               s_waveEvent.v_waveform.at(sample) =  v_wave.at(sample) - min;
               s_waveEvent.a_waveform[sample] =  v_wave.at(sample) - min;
               s_waveEvent.v_timing.at(sample) = sample*25.0 + (tdc-1000)*0.035;
               s_waveEvent.a_timing[sample] = s_waveEvent.v_timing.at(sample);
            }//for sample
            s_waveEvent.waveHeight = range;
            s_waveEvent.plane = plane;
            s_waveEvent.chip = chip;
            s_waveEvent.strip = strip;
            s_waveEvent.entry = entry;
            s_waveEvent.tdc = tdc;
            s_waveEvent.waveType = classifyWaveType(v_wave);
            m_chip[plane][chip].c_waveEvent.push_back(s_waveEvent);
         }//if v_wave
      }//for multi
      for(int hitNo=0; hitNo<s_eventHitID.v_stripID.size(); hitNo++){
         s_allEventID.v_eventHitID.at(entry).v_stripID.push_back(s_eventHitID.v_stripID.at(hitNo));
      }//for hitNo

      _entry s_entry;
      v_entry.push_back(s_entry);
   }//for entry


   //---------------------------------------------------------------

   //---------------------------------------------------------------


   ifs.close();
   ifs.open(input_param_filename);
   if(ifs.fail()){
      cerr << "cannnot read parameter file" << endl;
      return 0;
   }
   while(!(ifs.eof())){
      int plane, chip;
      double a_param[n_parameter];
      ifs >> plane >> chip >> a_param[0] >> a_param[1] >> a_param[2];
      for(int data=0; data<n_parameter; data++){
         m_chip[plane][chip].a_fitParam[data] = a_param[data];
      }//for data
   };

   for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
         for(int data=0; data<n_parameter; data++){
            cout << m_chip[plane][chip].a_fitParam[data] << " ";
         }//for data
         cout << endl;
      }//for chip
   }//for plane



   //-----------------------------



   int count=0;
   for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
         for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){

            m_chip[plane][chip].c_waveEvent.at(data).sample_max
               = getMaxSampleNum(m_chip[plane][chip].c_waveEvent.at(data).v_waveform_original);

            m_chip[plane][chip].c_waveEvent.at(data).sample_min
               = getMinSampleNum(m_chip[plane][chip].c_waveEvent.at(data).v_waveform_original);

            m_chip[plane][chip].c_waveEvent.at(data).timing_weighted_mean
               = getTimingWeightedMean(
                     m_chip[plane][chip].c_waveEvent.at(data).v_waveform,
                     m_chip[plane][chip].c_waveEvent.at(data).v_timing,
                     m_chip[plane][chip].c_waveEvent.at(data).sample_max
                     );

            int sampleMax
               = m_chip[plane][chip].c_waveEvent.at(data).sample_max;


            m_chip[plane][chip].c_waveEvent.at(data).adc_max
               = getMax(m_chip[plane][chip].c_waveEvent.at(data).v_waveform_original);

            m_chip[plane][chip].c_waveEvent.at(data).graph_wave
               = new TGraph(n_sample,
                     &m_chip[plane][chip].c_waveEvent.at(data).v_timing[0],
                     &m_chip[plane][chip].c_waveEvent.at(data).v_waveform_original[0]
                     );

            double fit_range_min, fit_range_max; 
            int waveType = m_chip[plane][chip].c_waveEvent.at(data).waveType;
            switch(m_chip[plane][chip].c_waveEvent.at(data).waveType){
               case 41:
                  //cout << waveType << endl;
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax-2);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+3);
                  break;
               case 411:
                  //cout << waveType << endl;
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax-1);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+3);
                  break;
               case 42:
                  //cout << waveType << endl;
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax-2);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+3);
                  break;
               case 421:
                  //cout << waveType << endl;
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax-1);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+3);
                  break;
               case 43:
                  //cout << waveType << endl;
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+3);
                  break;
               case 31: 
                  //cout << waveType << endl;
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax-2);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+3);
                  break;
               case 32:
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax-1);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+4);
                  break;
               case 33:
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax-2);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+4);
                  break;
               case 331:
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax-1);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+4);
                  break;
               case 51:
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax-2);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+2);
                  break;
               case 511:
                  fit_range_min
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax-1);
                  fit_range_max
                     = m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sampleMax+2);
                  break;
            }
            m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip
               = new TF1(Form("fit_p%d_c%d_d%d", plane, chip, data),
                     "[1]*(x-[2])/[0]*exp(-(x-[2])/[0])", fit_range_min-0.001, fit_range_max+0.001 
                     );

            double param_fix[n_parameter];
            for(int p_data=0; p_data<n_parameter; p_data++){
               param_fix[p_data] = m_chip[plane][chip].a_fitParam[p_data];
            }//for p_data

            m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip
               -> SetParameter(0, param_fix[2]*3.5);
            int minSample
               = m_chip[plane][chip].c_waveEvent.at(data).sample_min;
            m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip
               -> SetParameter(2, 
                     m_chip[plane][chip].c_waveEvent.at(data).timing_weighted_mean
                     - 140);

            m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip
               -> SetParameter(1,
                     m_chip[plane][chip].c_waveEvent.at(data).waveHeight * 2.718);


            int strip_num
               = m_chip[plane][chip].c_waveEvent.at(data).strip;


            m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
               = new TGraphErrors(
                     n_sample,
                     m_chip[plane][chip].c_waveEvent.at(data).a_timing,
                     m_chip[plane][chip].c_waveEvent.at(data).a_waveform,
                     m_chip[plane][chip].a_strip_calibdata[strip_num].a_timing_resolution,
                     m_chip[plane][chip].a_strip_calibdata[strip_num].a_calib_stdDev
                     );

            m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
               -> Fit(m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip, 
                     "OQMEX0", "", fit_range_min-0.001, fit_range_max+0.001);

            double chiSquare =
               m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip -> GetChisquare();
            int ndf =
               m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip -> GetNDF();
            double chiSquare_over_ndf
               = chiSquare/(double)ndf;


            m_chip[plane][chip].c_waveEvent.at(data).chiSquare = chiSquare;
            m_chip[plane][chip].c_waveEvent.at(data).ndf = ndf;
            m_chip[plane][chip].c_waveEvent.at(data).chiSquare_over_ndf = chiSquare_over_ndf;
            double param[n_parameter];
            double param_Error[n_parameter];
            for(int p_data=0; p_data<n_parameter; p_data++){
               param[p_data]
                  = m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip -> GetParameter(p_data);
               m_chip[plane][chip].c_waveEvent.at(data).a_fitParam[p_data] = param[p_data];
               param_Error[p_data]
                  = m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip -> GetParError(p_data);
               m_chip[plane][chip].c_waveEvent.at(data).a_fitParam_Error[p_data]
                  = param_Error[p_data];
            }//for p_data
            double fit_peak_time;
            double fit_peak_adc;
            double fit_hitTiming;
            double fit_peak_adc_Error;
            double fit_hitTiming_Error;
            fit_peak_adc = param[1] / 2.714;
            fit_peak_time = param[0] + param[2];
            fit_hitTiming = param[2];
            fit_peak_adc_Error = param_Error[1] / 2.714;
            fit_hitTiming_Error = param_Error[2];

            m_chip[plane][chip].c_waveEvent.at(data).fitted_peak_time = fit_peak_time;
            m_chip[plane][chip].c_waveEvent.at(data).fitted_peak_adc  = fit_peak_adc;
            m_chip[plane][chip].c_waveEvent.at(data).fitted_hitTiming  = fit_hitTiming;


            count++;

            if(count%1000 == 0){
               cerr << count << "/" << s_allchip.n_wave << endl;
            }

            //--------------------------
            //pushing back to entry structure
            _waveEvent s_waveEvent
               = m_chip[plane][chip].c_waveEvent.at(data);
            Int_t entry = s_waveEvent.entry;

            if(fit_hitTiming > 0 && fit_peak_adc<2000){
               v_entry.at(entry).m_chip[plane][chip].c_waveEvent.push_back(s_waveEvent);
            }
            //--------------------------

         }//for data
      }//for chip
   }//for plane


   //---------------------------------------
   //adding parameter to entry structure
   for(int entry=0; entry<n_entry; entry++){
      for(int chip=0; chip<n_chip; chip++){
         for(int plane=0; plane<n_plane; plane++){
            memcpy(v_entry.at(entry).m_chip[plane][chip].a_fitParam,
                  m_chip[plane][chip].a_fitParam,
                  sizeof(m_chip[plane][chip].a_fitParam));
            memcpy(v_entry.at(entry).m_chip[plane][chip].m_calib_mean,
                  m_chip[plane][chip].m_calib_mean,
                  sizeof(m_chip[plane][chip].m_calib_mean));
            memcpy(v_entry.at(entry).m_chip[plane][chip].m_calib_stdDev,
                  m_chip[plane][chip].m_calib_stdDev,
                  sizeof(m_chip[plane][chip].m_calib_stdDev));
         }//for plane
      }//for chip
   }//for entry
   //------------------------------------------





   //---------------------------------------------------------------
   //ADD CODE FOR CLUSTERING



   //add cluster num
   int arientry,ariplane,arichip,ariwave,ariwavesize; //plane,chip,wave counter
   //int arientry;
   int aristrip,aristrip_next; //strip number
   int c_num; //cluster number
   for(arientry=0; arientry<n_entry; arientry++){
      for(ariplane=0;ariplane<6;ariplane++){
         for(arichip=0;arichip<6;arichip++){
            ariwavesize = v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.size();
            c_num = 0;

            for(ariwave=0;ariwave<ariwavesize-1;ariwave++){
               //arientry = m_chip[ariplane][arichip].c_waveEvent.at(ariwave).entry;
               aristrip = v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(ariwave).strip;
               aristrip_next = v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(ariwave+1).strip;

               if(!(ariwave == ariwavesize-2)){
                  v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(ariwave).cluster_num = c_num;
                  if(!((aristrip+1) == aristrip_next)){
                     c_num++;
                  }
               }
               else if(ariwave == ariwavesize-2){
                  v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(ariwave).cluster_num = c_num;
                  if((aristrip+1) == aristrip_next){
                     v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(ariwave+1).cluster_num = c_num;
                  }
                  else{
                     v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(ariwave+1).cluster_num = c_num+1;
                  }
               }


            }//for ariwave
         }//for arichip
      }//for ariplane
   }//for arientry


   //_cluster
   int wave_num;
   int cl_wave_num=0;
   int wave_num_max=0;
   int cluster_num_max=-1;
   double adc_maxC=0;
   double adc_max_fit=0;
   int numerator_tdcAvr=0;
   int simplenumerator_stripAvr=0;
   double numerator_stripAvr=0;
   double numerator_chiSquareSum=0;
   double numerator_ndfSum=0;
   int cl_num = 0;
   //   std::cout << "entry " << "cnum " << "plane " << "chip " << std::endl;

   for(arientry=0; arientry<n_entry; arientry++){
      for(ariplane=0;ariplane<6;ariplane++){
         for(arichip=0;arichip<6;arichip++){
            wave_num_max = v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.size();
            if(wave_num_max==0){
               cluster_num_max=-1;
            }
            else{
               cluster_num_max = v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num_max-1).cluster_num;

               //adc,tdc,strip clusterize
               v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.resize(cluster_num_max+1);

               for(cl_num=0;cl_num<cluster_num_max+1;cl_num++){
                  for(wave_num=0;wave_num<wave_num_max;wave_num++){
                     if(v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).cluster_num==cl_num){
                        v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).w_num++;
                        v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).c_cluster_num = c_num;
                        v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).entry_cluster = v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).entry;
                        //adc_maxC += v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).adc_max;
                        adc_maxC += v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).fitted_peak_adc;
                        adc_max_fit += v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).fitted_peak_adc;
                        //numerator_tdcAvr += v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).tdc * v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).adc_max;
                        numerator_chiSquareSum += v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).chiSquare;
                        numerator_ndfSum += v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).ndf;
                        numerator_tdcAvr += v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).fitted_hitTiming * v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).adc_max;
                        //cout << "1 " << v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).adc_max << endl;
                        simplenumerator_stripAvr += v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).strip * v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).adc_max;
                        numerator_stripAvr += v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).strip * v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(wave_num).fitted_peak_adc;
                        cl_wave_num++;
                        v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).cl_waveEvent.resize(cl_wave_num);
                        v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).cl_waveEvent.at(cl_wave_num-1)=m_chip[ariplane][arichip].c_waveEvent.at(wave_num);
                     }//if cl_num
                     cl_wave_num=0;//initialize
                  }//for wave_num
                  v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).adc_max_cluster = adc_maxC;
                  if(!(adc_maxC==0)){
                     //			   std::cout << "e,p,c=" << arientry << ariplane << arichip << " " << "tdc_numerator=" << numerator_tdcAvr << " " << "adc_sum=" << adc_maxC << std::endl;
                     v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).tdc_cluster = numerator_tdcAvr / (Double_t)adc_maxC;
                     v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).chiSquareSum = numerator_chiSquareSum;
                     v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).ndfSum = numerator_ndfSum;
                     
                     //cout << v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).tdc_cluster << endl;
                     v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).simple_strip_cluster = simplenumerator_stripAvr / adc_maxC;
                     v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(cl_num).strip_cluster = numerator_stripAvr / adc_max_fit;
                  }//if adcmax !=0
                  else{
                     std::cout << "error! adc = 0" << std::endl;
                  }
                  adc_maxC=0;//initialize adc_maxC
                  numerator_tdcAvr=0;//initialize numerator
                  simplenumerator_stripAvr=0;//initialize numerator
                  numerator_stripAvr=0;//initialize numerator
                  numerator_chiSquareSum=0;
                  numerator_ndfSum=0;
                  adc_max_fit=0;//initialize

                  //std::cout << arientry << " " << v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(c_num).c_cluster_num << " " << ariplane << " " << arichip << std::endl;


               }//for cl_num
            }//if wave num !=0(else)
         }//arichip
      }//ariplane
   }//for arientry



   //---------------------------------------------------------------








   _plane c_plane[n_plane];
   for(int plane=0; plane<n_plane; plane++){
      int stripSize = n_chip * n_strip;
      c_plane[plane].v_stripNo.resize(stripSize);
      c_plane[plane].v_stripNo_Error.resize(stripSize);
      c_plane[plane].v_numOfHit.resize(stripSize);
      c_plane[plane].v_numOfHit_Error.resize(stripSize);
      for(int chip=0; chip<n_chip; chip++){
         for(int strip=0; strip<n_strip; strip++){
            int stripNo = strip + chip*n_strip;
            c_plane[plane].v_stripNo.at(stripNo) = stripNo;
            c_plane[plane].v_stripNo_Error.at(stripNo) = 0;
            c_plane[plane].v_numOfHit.at(stripNo) = 0;
            c_plane[plane].v_numOfHit_Error.at(stripNo) = 0;
            c_plane[plane].a_stripNo[stripNo] = stripNo;
            c_plane[plane].a_stripNo_Error[stripNo] = 0;
            c_plane[plane].a_numOfHit[stripNo] = 0;
            c_plane[plane].a_numOfHit_Error[stripNo] = 0;
         }//for strip
      }//for chip
   }//for plane

   for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
         for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
            double chiSquare = m_chip[plane][chip].c_waveEvent.at(data).chiSquare;
            int strip = m_chip[plane][chip].c_waveEvent.at(data).strip;
            int stripNo = strip + chip*n_strip;

            c_plane[plane].v_numOfHit.at(stripNo) += 1.0;
            c_plane[plane].a_numOfHit[stripNo] += 1,0;
            c_plane[plane].hitSumCount += 1.0;


         }//for data
      }//for chip
   }//for data

   for(int plane=0; plane<n_plane; plane++){
      c_plane[plane].hitAverage = c_plane[plane].hitSumCount/(double)n_stripOnPlane;
      for(int chip=0; chip<n_chip; chip++){
         for(int strip=0; strip<n_strip; strip++){
            int stripNo = strip + chip*n_strip;
            int numOfHit = c_plane[plane].v_numOfHit.at(stripNo);
            //double numOfHit_Error = sqrt(numOfHit);
            double numOfHit_Error = 0;
            c_plane[plane].v_numOfHit_Error.at(stripNo) = numOfHit_Error;
            c_plane[plane].a_numOfHit_Error[stripNo] = numOfHit_Error;
         }//for strip
      }//for chip
   }//for plane

   for(int plane=0; plane<n_plane; plane++){
      c_plane[plane].hit_map
         = new TGraphErrors(
               n_stripOnPlane,
               c_plane[plane].a_stripNo,
               c_plane[plane].a_numOfHit,
               c_plane[plane].a_stripNo_Error,
               c_plane[plane].a_numOfHit_Error
               );

   }//for plane



   //TCanvas *c1 = new TCanvas("c1", "My Canvas", 0, 0, 596, 804);
   //TCanvas *c_a4v = new TCanvas("c_a4v", "My Canvas a4v", 0, 0, 596, 804);
   //TCanvas *c_a4h = new TCanvas("c_a4h", "My Canvas a4h", 0, 0, 804, 596);
   //TCanvas *c_b5v = new TCanvas("c_b5v", "My Canvas b5v", 0, 0, 512, 728);
   //TCanvas *c_b5h = new TCanvas("c_b5h", "My Canvas b5h", 0, 0, 728, 512);


   /*
      int drawCount=0;
      for(int data=0; data<m_chip[0][0].c_waveEvent.size(); data++){
      if(drawCount==0){
      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> SetTitle("SSD signal waveform");
      m_chip[0][0].c_waveEvent.at(data).graph_wave
      -> GetXaxis() -> SetTitle("timing[#musec]");
      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> GetYaxis() -> SetTitle("adc");

      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> SetMarkerColor(4);

      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> SetMarkerStyle(21);


      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> SetMarkerSize(1);

      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error->Draw("AP");
      drawCount++;
      }//
      }//for data
      */



   //c1->Divide(6,6);

   /*
      s_allchip.h_chi_square_over_ndf
      = new TH1D("chiSquare/ndf", "chiSquare/ndf", 40, 0, 40);
      s_allchip.h_chi_square
      = new TH1D("chiSquare", "chiSquare", 40, 0, 40);

      for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
      for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
      double waveHeight
      = m_chip[plane][chip].c_waveEvent.at(data).waveHeight;
      s_allchip.h_chi_square_over_ndf
      -> Fill(m_chip[plane][chip].c_waveEvent.at(data).chiSquare_over_ndf);
      if(30 < waveHeight && waveHeight < 40){
      s_allchip.h_chi_square
      -> Fill(m_chip[plane][chip].c_waveEvent.at(data).chiSquare);
      }//if waveHeight
      }//for data
      }//for chip
      }//for plane


      s_allchip.h_chi_square -> SetTitle("Fitting #chi^{2}");
      s_allchip.h_chi_square  -> GetXaxis() -> SetTitle("#chi^{2}");
      s_allchip.h_chi_square  -> GetYaxis() -> SetTitle("NumOfEvents");
      s_allchip.h_chi_square->Draw();



      cout << "ndf = " << m_chip[0][0].c_waveEvent.at(0).ndf << endl;
      */

   /*
      s_allchip.h_peak_time
      = new TH1D("peaktime", "peaktime", 200, 50, 150);

      for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
      for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
      if(m_chip[plane][chip].c_waveEvent.at(data).chiSquare < 7.815){
      s_allchip.h_peak_time
      -> Fill(m_chip[plane][chip].c_waveEvent.at(data).fitted_peak_time);
      }//if chiSquare
      }//for data
      }//for chip
      }//for plane

      s_allchip.h_peak_time -> SetTitle("Fitted peaktime");
      s_allchip.h_peak_time  -> GetXaxis() -> SetTitle("peaktime[nsec]");
      s_allchip.h_peak_time  -> GetYaxis() -> SetTitle("num of events");
      s_allchip.h_peak_time -> Draw();

      cout << "ndf = " << m_chip[0][0].c_waveEvent.at(0).ndf << endl;
      */

   /*
      s_allchip.h_raising_time
      = new TH1D("raising_time", "rising_time", 100, 5, 65);

      for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
      for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
      s_allchip.h_raising_time
      -> Fill(m_chip[plane][chip].c_waveEvent.at(data).a_fitParam[2]);
      }//for data
      }//for chip
      }//for plane

      s_allchip.h_raising_time -> SetTitle("RisingTime");
      s_allchip.h_raising_time  -> GetXaxis() -> SetTitle("risingTime[nsec]");
      s_allchip.h_raising_time  -> GetYaxis() -> SetTitle("num of events");
      s_allchip.h_raising_time -> Draw();
      */




   /*
      for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
      m_chip[plane][chip].h_peak_time
      = new TH1D(Form("peakTime_p%d_c%d", plane, chip),
      Form("peakTime_p%d_c%d", plane, chip),
      60, 60, 120);
      for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
      if(m_chip[plane][chip].c_waveEvent.at(data).chiSquare < 16.268){
      m_chip[plane][chip].h_peak_time
      -> Fill(m_chip[plane][chip].c_waveEvent.at(data).fitted_peak_time);
      }//if chiSquare
      }//for data
   //c1->cd(plane*n_plane+chip+1);
   //m_chip[plane][chip].h_peak_time -> Draw();
   }//for chip
   }//for plane
   m_chip[plane][chip].h_peak_time->SetTitle("Fitted peaktime");
   m_chip[plane][chip].h_peak_time-> GetXaxis() -> SetTitle("peaktime[nsec]");
   m_chip[plane][chip].h_peak_time-> GetYaxis() -> SetTitle("NumOfEvents");
   m_chip[plane][chip].h_peak_time->Draw();
   */

   /*
      int drawCount=0;
      for(int data=0; data<m_chip[0][0].c_waveEvent.size(); data++){
      if(drawCount==0){
      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> SetTitle("SSD signal waveform");
      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> GetXaxis() -> SetTitle("timing[nsec]");
      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> GetYaxis() -> SetTitle("adc incremented");

      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> SetMarkerColor(4);

      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> SetMarkerStyle(21);


      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error
      -> SetMarkerSize(1.5);


      TLegend *legend = new TLegend(0.65, 0.20, 0.95, 0.38);
      legend -> SetTextSize(0.035);
      legend->AddEntry(m_chip[0][0].c_waveEvent.at(data).graph_wave_Error, "dataplot");
      legend->AddEntry(m_chip[0][0].c_waveEvent.at(data).fit_func_chip, "FitCurve");

      m_chip[0][0].c_waveEvent.at(data).graph_wave_Error->Draw("AP");
      gStyle -> SetOptStat(111111111);
      m_chip[0][0].c_waveEvent.at(data).fit_func_chip->Draw("same");
      legend->Draw("same");
      drawCount++;
      }//if drawCount
      }//for data
      */

   /*
      c1->Divide(8,8);
      int canvNum=0;
      for(int data=0; data<m_chip[0][0].c_waveEvent.size(); data++){
      if(canvNum<64){
      if(m_chip[0][0].c_waveEvent.at(data).chiSquare > 30){
      c1->cd(canvNum+1);

      m_chip[0][0].c_waveEvent.at(data).graph_wave
      -> SetMarkerColor(4);

      m_chip[0][0].c_waveEvent.at(data).graph_wave
      -> SetMarkerStyle(21);


      m_chip[0][0].c_waveEvent.at(data).graph_wave
      -> SetMarkerSize(0.1);


      m_chip[0][0].c_waveEvent.at(data).graph_wave->Draw("AP");
      m_chip[0][0].c_waveEvent.at(data).fit_func_chip->Draw("same");
      gStyle -> SetOptStat(111111111);
      canvNum++;
      }//if ChiSquare
      }//if canvNum
      }//for data
      */

   /*
      c1->Divide(5,5);
      int canvNum=0;
      for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
      for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
      if(canvNum<25){
      c1->cd(canvNum+1);
      m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
      -> SetMarkerColor(4);
      m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
      -> SetMarkerStyle(21);
      m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
      -> SetMarkerSize(0.5);
      m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
      -> GetYaxis() -> SetLabelSize(0.1);
      m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
      -> GetXaxis() -> SetLabelSize(0.08);
   //m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
   //-> GetYaxis() -> SetTitleSize(0.1);
   m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error->Draw("AP");
   m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip->Draw("same");
   canvNum++;
   }//for canvNum
   else{
   break;
   }
   }//for data
   }//for chip
   }//for plane
   */


   /*
      for(int entry=0; entry<s_allEventID.v_eventHitID.size(); entry++){
      for(int hitNo=0; hitNo<s_allEventID.v_eventHitID.at(entry).v_stripID.size(); hitNo++){
      int plane, chip, strip, orderInChip;
      plane = s_allEventID.v_eventHitID.at(entry).v_stripID.at(hitNo).plane;
      chip = s_allEventID.v_eventHitID.at(entry).v_stripID.at(hitNo).chip;
      strip = s_allEventID.v_eventHitID.at(entry).v_stripID.at(hitNo).strip;
      orderInChip = s_allEventID.v_eventHitID.at(entry).v_stripID.at(hitNo).orderInChip;
      cout << entry << " " << hitNo << " " << plane << " " << chip << " " << strip << " " << orderInChip <<  endl;
      }//for hitNo
      }//for entry
      */
   /*
      vector<double> chiSquare_Thr = {0, 15, 30, 45, 60, 75, 90, 120, 999};


      c_a4v->SaveAs(output_pdf_name_a4+"[", "pdf");
      c_b5h->SaveAs(output_pdf_name_b5+"[", "pdf");



      c_b5h->cd();
      s_allchip.h_chi_square
      = new TH1D("chiSquare", "chiSquare", 60, 0, 30);
      for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
      for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
      s_allchip.h_chi_square
      -> Fill(m_chip[plane][chip].c_waveEvent.at(data).chiSquare);
      }//for data
      }//for chip
      }//for plane

      s_allchip.h_chi_square -> SetTitle("Fitting #chi^{2}");
      s_allchip.h_chi_square -> GetXaxis() -> SetTitle("#chi^{2}");
      s_allchip.h_chi_square -> GetYaxis() -> SetTitle("NumOfEvents");
      s_allchip.h_chi_square -> Draw();
      c_b5h->SaveAs(output_pdf_name_b5, "pdf");
      c_b5h->Clear();



      c_b5h->cd();

      s_allchip.h_chi_square_over_ndf
      = new TH1D("chiSquare/ndf", "chiSquare/ndf", 45, 0, 15);
      for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
      for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
      s_allchip.h_chi_square_over_ndf
      -> Fill(m_chip[plane][chip].c_waveEvent.at(data).chiSquare_over_ndf);
      }//for data
      }//for chip
      }//for plane

      s_allchip.h_chi_square_over_ndf -> SetTitle("Fitting #chi^{2}/ndf");
      s_allchip.h_chi_square_over_ndf -> GetXaxis() -> SetTitle("#chi^{2}/ndf");
      s_allchip.h_chi_square_over_ndf -> GetYaxis() -> SetTitle("NumOfEvents");
      s_allchip.h_chi_square_over_ndf -> Draw();
      c_b5h->SaveAs(output_pdf_name_b5, "pdf");
      c_b5h->Clear();




      c_b5h->cd();
      s_allchip.h_chiSquare_vs_waveHeight
      = new TH2D("chiSquare vs waveHeight", "chiSquare vs waveHeight", 24, 0, 12, 30, 20, 150);

      for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
      for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
      double chiSquare
      = m_chip[plane][chip].c_waveEvent.at(data).chiSquare;
      double waveHeight
      = m_chip[plane][chip].c_waveEvent.at(data).waveHeight;
      s_allchip.h_chiSquare_vs_waveHeight
      -> Fill(chiSquare, waveHeight);
      }//for data
      }//for chip
      }//for plane
      s_allchip.h_chiSquare_vs_waveHeight
      -> GetXaxis() -> SetTitle("#chi^{2}");
   s_allchip.h_chiSquare_vs_waveHeight
      -> GetYaxis() -> SetTitle("waveHeight(adc)");

   s_allchip.h_chiSquare_vs_waveHeight->Draw("colz");

   c_b5h->SaveAs(output_pdf_name_b5, "pdf");
   c_b5h->Clear();


   c_b5h->cd();
   s_allchip.h_waveHeight
      = new TH1D("waveHeight", "waveHeight", 80, 20, 100);
   for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
         for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
            double waveHeight
               = m_chip[plane][chip].c_waveEvent.at(data).waveHeight;
            s_allchip.h_waveHeight -> Fill(waveHeight);
         }//for data
      }//for chip
   }//for plane
   s_allchip.h_waveHeight
      -> GetXaxis() -> SetTitle("waveHeight(adc)");
   s_allchip.h_waveHeight
      -> GetYaxis() -> SetTitle("numOfEvents");
   s_allchip.h_waveHeight->Draw();
   c_b5h->SaveAs(output_pdf_name_b5, "pdf");
   c_b5h->Clear();




   c_b5h->cd();
   s_allchip.h_timeZero
      = new TH1D("timeZero", "timeZero", 120, 0, 120);

   for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
         for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
            double chiSquare
               = m_chip[plane][chip].c_waveEvent.at(data).chiSquare;
            if(chiSquare<15){
               s_allchip.h_timeZero
                  -> Fill(m_chip[plane][chip].c_waveEvent.at(data).a_fitParam[2]);
            }  
         }//for data
      }//for chip
   }//for plane

   s_allchip.h_timeZero -> SetTitle("timeZero (#chi^{2}<15)");
   s_allchip.h_timeZero  -> GetXaxis() -> SetTitle("T_{0}[nsec]");
   s_allchip.h_timeZero  -> GetYaxis() -> SetTitle("num of events");
   s_allchip.h_timeZero -> Draw();

   c_b5h->SaveAs(output_pdf_name_b5, "pdf");
   c_b5h->Clear();




   c_b5h->cd();

   s_allchip.hs_timeZero = new THStack("timeZeroStack", "timeZero");
   s_allchip.h_timeZero_Max3 = new TH1D("timeZero_Max3", "timeZero_Max3", 120, 0, 120);
   s_allchip.h_timeZero_Max4 = new TH1D("timeZero_Max4", "timeZero_Max4", 120, 0, 120);
   s_allchip.h_timeZero_Max5 = new TH1D("timeZero_Max5", "timeZero_Max5", 120, 0, 120);

   for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
         for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
            double chiSquare
               = m_chip[plane][chip].c_waveEvent.at(data).chiSquare;
            int maxSampleNum
               = m_chip[plane][chip].c_waveEvent.at(data).sample_max;
            double timeZero =
               m_chip[plane][chip].c_waveEvent.at(data).a_fitParam[2];
            if(chiSquare < 15){
               if(maxSampleNum ==3){
                  s_allchip.h_timeZero_Max3->Fill(timeZero);
               }
               else if(maxSampleNum == 4){
                  s_allchip.h_timeZero_Max4->Fill(timeZero);
               }
               else if(maxSampleNum == 5){
                  s_allchip.h_timeZero_Max5->Fill(timeZero);
               }
               else{
                  cerr << "maxSampleNum is maybe invalid" << endl;
                  return 0;
               }
            }
         }//for data
      }//for chip
   }//for plane

   cout << "test1" << endl;

   s_allchip.h_timeZero_Max3->SetFillColor(kRed);
   s_allchip.h_timeZero_Max4->SetFillColor(kBlue);
   s_allchip.h_timeZero_Max5->SetFillColor(kGreen);

   s_allchip.hs_timeZero -> Add(s_allchip.h_timeZero_Max3);
   s_allchip.hs_timeZero -> Add(s_allchip.h_timeZero_Max4);
   s_allchip.hs_timeZero -> Add(s_allchip.h_timeZero_Max5);

   s_allchip.hs_timeZero->Draw();

   s_allchip.hs_timeZero -> SetTitle("timeZero (#chi^{2}<15)");
   s_allchip.hs_timeZero  -> GetXaxis() -> SetTitle("T_{0}[nsec]");
   s_allchip.hs_timeZero  -> GetYaxis() -> SetTitle("num of events");

   c_b5h->Modified();

   c_b5h->SaveAs(output_pdf_name_b5, "pdf");
   c_b5h->Clear();



   c_b5h->Divide(3,2);
   int canvNum=0;
   for(int plane=0; plane<n_plane; plane++){
      c_b5h->cd(canvNum+1);
      if(plane<3){
         c_plane[plane].hit_map
            -> SetTitle(Form("sensorRIGHT%d", plane%3));
      }
      else{
         c_plane[plane].hit_map
            -> SetTitle(Form("sensorLEFT%d", plane%3));
      }
      c_plane[plane].hit_map
         -> GetXaxis() -> SetTitle("strip");
      c_plane[plane].hit_map
         -> GetYaxis() -> SetTitle("hitCount");
      c_plane[plane].hit_map
         -> SetMarkerColor(4);
      c_plane[plane].hit_map
         -> SetMarkerStyle(21);
      c_plane[plane].hit_map
         -> SetMarkerSize(0.3);
      c_plane[plane].hit_map->Draw("AP");
      canvNum++;
   }//for plane
   c_b5h->SaveAs(output_pdf_name_b5, "pdf");
   c_b5h->Clear();



   c_a4v->cd();
   for(int t_data=0; t_data < chiSquare_Thr.size()-1; t_data++){

      c_a4v->Divide(5,5);
      int canvNum=0;
      for(int plane=0; plane<n_plane; plane++){
         for(int chip=0; chip<n_chip; chip++){
            for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
               double chiSquare
                  = m_chip[plane][chip].c_waveEvent.at(data).chiSquare;
               if(canvNum<25 && chiSquare_Thr.at(t_data) < chiSquare
                     && chiSquare < chiSquare_Thr.at(t_data + 1)
                 ){
                  c_a4v->cd(canvNum+1);
                  m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
                     -> SetMarkerColor(4);
                  m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
                     -> SetMarkerStyle(21);
                  m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
                     -> SetMarkerSize(0.5);
                  m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
                     -> GetYaxis() -> SetLabelSize(0.06);
                  m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
                     -> GetXaxis() -> SetLabelSize(0.04);
                  //m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error
                  //-> GetYaxis() -> SetTitleSize(0.1);
                  m_chip[plane][chip].c_waveEvent.at(data).graph_wave_Error->Draw("AP");
                  m_chip[plane][chip].c_waveEvent.at(data).fit_func_chip->Draw("same");
                  canvNum++;
               }//if canvNum
            }//for data
         }//for chip
      }//for plane
      c_a4v->SaveAs(output_pdf_name_a4, "pdf");
      c_a4v->Clear();
   }//fot t_data


   c_a4v->SaveAs(output_pdf_name_a4+"]", "pdf");
   c_b5h->SaveAs(output_pdf_name_b5+"]", "pdf");


   */

      //--- output file -----------------------------------------------------
      /*
      //for debug
      int ent,clu,wav,pla,cp,str,ad,td;

      ofstream outputfile("rawdata.dat");
      outputfile << "entry" << " " << "cluster_num" << " " << "wave_num" << " " << "plane" << " " << "chip" << " " << "strip" << " " << "adc" << " " << "tdc" << std::endl;
      for(arientry=0; arientry<300; arientry++){
      for(ariplane=0; ariplane<6; ariplane++){
      for(arichip=0; arichip<6; arichip++){
      for(ariwave=0; ariwave<v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.size(); ariwave++){
      ent=arientry;
      clu=v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(ariwave).cluster_num;
      wav=ariwave;
      pla=ariplane;
      cp=arichip;
      str=v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(ariwave).strip;
      ad=v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(ariwave).adc_max;
      td=v_entry.at(arientry).m_chip[ariplane][arichip].c_waveEvent.at(ariwave).tdc;

      outputfile << ent << " " << clu << " " << wav <<  " " << pla << " " << cp << " " << str << " " << ad << " " << td << std::endl;
      }
      }
      }
      }
      */


      TFile *fout = new TFile("SSDdata2999.root","recreate");
   TTree *out_tree = new TTree("tree","Cluster Tree");
   Int_t Entry;
   std::vector<Int_t> Cluster_num,Cluster_size,module,Plane,Chip;
   //   TVector3 PlaneChipStrip;
   //   TVector2 XZ;
   std::vector<Double_t> Strip,DecaltX,DecaltZ,ADC,TDC,chiSquareSum,ndfSum;

   out_tree->Branch("Cluster_num", &Cluster_num);
   out_tree->Branch("Cluster_size", &Cluster_size);
   out_tree->Branch("Plane", &Plane);
   out_tree->Branch("Chip", &Chip);
   out_tree->Branch("Strip", &Strip);
   //out_tree->Branch("module", &module, "module/I");
   out_tree->Branch("DecaltX", &DecaltX);
   out_tree->Branch("DecaltZ", &DecaltZ);
   out_tree->Branch("ADC", &ADC);
   out_tree->Branch("TDC", &TDC);
   out_tree->Branch("chiSquareSum", &chiSquareSum);
   out_tree->Branch("ndfSum", &ndfSum);

   int CLUSTERinENTRY=0;
   int entry_max;
   int c_num_max;
   int c_counter=0;
   int c_num_entry=0;

   for(arientry=0; arientry<n_entry; arientry++){
      for(ariplane=0; ariplane<6; ariplane++){
         for(arichip=0; arichip<6; arichip++){
            c_num_max = v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.size();
            for(c_num=0; c_num<c_num_max; c_num++){
               Cluster_num.push_back(c_num_entry);
               Cluster_size.push_back(v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(c_num).w_num);
               Plane.push_back(ariplane);
               Chip.push_back(arichip);
               Strip.push_back(v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(c_num).strip_cluster);
               ADC.push_back(v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(c_num).adc_max_cluster);
               TDC.push_back(v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(c_num).tdc_cluster);
               cout << v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(c_num).tdc_cluster << endl;
               DecaltX.push_back(calcPosDecalt((double)ariplane,(double)arichip,v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(c_num).strip_cluster).X());
               DecaltZ.push_back(calcPosDecalt((double)ariplane,(double)arichip,v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(c_num).strip_cluster).Y());
               chiSquareSum.push_back(v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(c_num).chiSquareSum);
               ndfSum.push_back(v_entry.at(arientry).m_chip[ariplane][arichip].v_cluster.at(c_num).ndfSum);
               /*
                  std::cout << "entry=" << arientry;
                  std::cout << " c_num=" << Cluster_num.at(c_counter);
                  std::cout << " plane=" << Plane.at(c_counter);
                  std::cout << " chip=" << Chip.at(c_counter) << std::endl;
                  */
               //Strip.push_back(???);
               c_num_entry++;
               c_counter++;
            }//for c_num
         }//for arichip
      }//for ariplane
      out_tree->Fill();
      /*
         std::cout << "entry=" << arientry << std::endl;
         for(int uno=0; uno<ADC.size(); uno++){
         std::cout << Cluster_num.at(uno) << " " << Wave_num.at(uno) << " " << Plane.at(uno) << " " << Chip.at(uno) << " " << Strip.at(uno) << " " << ADC.at(uno) << " " << TDC.at(uno) << " " << DecaltX.at(uno) << " " << DecaltZ.at(uno) << std::endl;
         }
         */
      c_num_entry=0;//initialize
      Cluster_num.clear();
      Cluster_size.clear();
      Plane.clear();
      Chip.clear();
      Strip.clear();
      ADC.clear();
      TDC.clear();
      DecaltX.clear();
      DecaltZ.clear();
   }//for arientry

   out_tree->Write();
   fout->Close();




   // for(ariplane=0; ariplane<6; ariplane++){
   //     for(arichip=0; arichip<6; arichip++){
   // 	   c_num_max = m_chip[ariplane][arichip].v_cluster.size();
   // 	   if(c_num_max != 0){
   // 	       entry_max = m_chip[ariplane][arichip].v_cluster.at(c_num_max-1).entry_cluster;
   // 	       std::cout << "cnummax,entrymax=" << c_num_max << "," << entry_max << std::endl;
   // 	       for(int Entry=0; Entry<entry_max+1; Entry++){
   // 		   for(c_num=0; c_num<c_num_max; c_num++){
   // 		       if(m_chip[ariplane][arichip].v_cluster.at(c_num).entry_cluster == Entry){
   // 			   CLUSTERinENTRY++;
   // 		       }//if Entry
   // 		   }//for c_num

   // 		   Cluster_num.resize(CLUSTERinENTRY);
   // 		   //Wave_num.resize(CLUSTERinENTRY);
   // 		   Plane.resize(CLUSTERinENTRY);
   // 		   Chip.resize(CLUSTERinENTRY);
   // 		   Strip.resize(CLUSTERinENTRY);
   // 		   //module.resize(CLUSTERinENTRY);
   // 		   DecaltX.resize(CLUSTERinENTRY);
   // 		   DecaltZ.resize(CLUSTERinENTRY);
   // 		   ADC.resize(CLUSTERinENTRY);
   // 		   TDC.resize(CLUSTERinENTRY);
   // 		   CLUSTERinENTRY=0;//initialize


   // 		   std::cout << "size=" << Cluster_num.size() << std::endl;

   // 		   int AT=0;
   // 		   for(c_num=0; c_num<c_num_max; c_num++){
   // 		       if(m_chip[ariplane][arichip].v_cluster.at(c_num).entry_cluster == Entry){
   // 			   Cluster_num.at(AT)=m_chip[ariplane][arichip].v_cluster.at(c_num).c_cluster_num;
   // 			   //Wave_num.at(AT)=???;
   // 			   Plane.at(AT)=ariplane;
   // 			   Chip.at(AT)=arichip;
   // 			   Strip.at(AT)=m_chip[ariplane][arichip].v_cluster.at(c_num).strip_cluster;
   // 			   //module.at(AT)=???;
   // 			   DecaltX.at(AT)=calcPosDecalt((double)ariplane,(double)arichip,m_chip[ariplane][arichip].v_cluster.at(c_num).strip_cluster).X();
   // 			   DecaltZ.at(AT)=calcPosDecalt((double)ariplane,(double)arichip,m_chip[ariplane][arichip].v_cluster.at(c_num).strip_cluster).Y();
   // 			   AT++;
   // 		       }//if Entry(2kaime)
   // 		   }//for c_num(2kaime)

   // 		   out_tree->Fill();

   // 	       }//for Entry
   // 	   }//if cmax!=0
   //     }//for arichip
   // }//for ariplane

   // out_tree->Write();
   // fout->Close();

   //    for(ariplane=0; ariplane<6; ariplane++){
   //        for(arichip=0; arichip<6; arichip++){
   // 	   cluster_num_max = m_chip[ariplane][arichip].v_cluster.size();
   // 	   w_num_max = m_chip[ariplane][arichip].c_waveEvent.size();
   // 	   for(c_num=0; c_num<cluster_num_max; c_num++){
   // 	       for(int w_num=0; w_num<m_chip[ariplane][arichip].c_waveEvent.size(); w_num++){
   // 		   if(m_chip[ariplane][arichip].c_waveEvent.at(w_num).cluster_num==c_num){
   // 		       Wave_num++;
   // 		   }//if cluster_num == c_num
   // 		   for(int Entry=0; Entry<m_chip[ariplane][arichip].c_waveEvent.at(w_num_max).entry; Entry++){
   // 		       if(m_chip[ariplane][arichip].v_cluster.at(w_num).cluster_entry == Entry){
   // 		       }
   // 		   }//for Entry
   // 	       }//for w_num
   // 	       Cluster_num = m_chip[ariplane][arichip].v_cluster.at(c_num).c_cluster_num;
   // 	       Plane = ariplane;
   // 	       Chip = arichip;
   // 	       Strip = m_chip[ariplane][arichip].v_cluster.at(c_num).strip_cluster;
   // //	       PlaneChipStrip.SetX((double)ariplane);
   // //	       PlaneChipStrip.SetY((double)arichip);
   // //	       PlaneChipStrip.SetZ(m_chip[ariplane][arichip].v_cluster.at(c_num).strip_cluster);
   // 	       DecaltX = calcPosDecalt((double)ariplane,(double)arichip,m_chip[ariplane][arichip].v_cluster.at(c_num).strip_cluster).X();
   // 	       DecaltZ = calcPosDecalt((double)ariplane,(double)arichip,m_chip[ariplane][arichip].v_cluster.at(c_num).strip_cluster).Y();
   // //	       XZ.SetX(DecaltX);
   // //	       XZ.SetY(DecaltZ);
   // 	       ADC = m_chip[ariplane][arichip].v_cluster.at(c_num).adc_max_cluster;
   // 	       TDC = m_chip[ariplane][arichip].v_cluster.at(c_num).tdc_cluster;
   // 	       out_tree->Fill();

   // 	       Wave_num=0;//initialize

   // 	   }//for c_num
   //        }//for arichip
   //    }//for ariplane

   //   out_tree->Write();
   //   fout->Close();


   //--- end output file -------------------------------------------------




}//waveAnalyzer







bool isMaxAt3(vector<double> v_waveform){
   int max_sample=0;
   int max_adc=0;
   for(int sample=0; sample<v_waveform.size(); sample++){
      if(max_adc < v_waveform.at(sample)){
         max_sample = sample;
         max_adc = v_waveform.at(sample);
      }
   }
   if(max_sample == 3){
      return true;
   }
   else{
      return false;
   }

};


bool isMaxAt4(vector<double> v_waveform){
   int max_sample=0;
   int max_adc=0;
   for(int sample=0; sample<v_waveform.size(); sample++){
      if(max_adc < v_waveform.at(sample)){
         max_sample = sample;
         max_adc = v_waveform.at(sample);
      }
   }
   if(max_sample == 4){
      return true;
   }
   else{
      return false;
   }

};


bool isMaxAt5(vector<double> v_waveform){
   int max_sample=0;
   int max_adc=0;
   for(int sample=0; sample<v_waveform.size(); sample++){
      if(max_adc < v_waveform.at(sample)){
         max_sample = sample;
         max_adc = v_waveform.at(sample);
      }
   }
   if(max_sample == 5){
      return true;
   }
   else{
      return false;
   }

};

bool isBiggerAt8Than0(vector<double> v_waveform){
   if(v_waveform.at(0) < v_waveform.at(v_waveform.size()-1)){
      return true;
   }
   else{
      return false;
   }
};


bool isTailLow(vector<double> v_waveform){
   double min = 99999;
   double min_sample = 100;
   int n_sample = v_waveform.size();
   for(int sample=0; sample<n_sample; sample++){
      if(min > v_waveform.at(sample)){
         min = v_waveform.at(sample);
         min_sample = sample;
      }
   } 
   if(v_waveform.at(min_sample+1) > v_waveform.at(n_sample-1)){
      return true;
   }
   else{
      return false;
   }
};

bool isNotNoisyStrip(int plane, int chip,int strip){
   int stripNo = chip*n_strip + strip;
   if(plane == 4 && stripNo == 383){
      return false;
   }
   else{
      return true;
   }
};

double  getMax(vector<double> v_waveform){
   double max = -99999;
   for(int sample=0; sample<v_waveform.size(); sample++){
      if(v_waveform.at(sample) > max){
         max = v_waveform.at(sample);
      }
   }
   return max;
};

double getMin(vector<double> v_waveform){
   double min = 99999;
   for(int sample=0; sample<v_waveform.size(); sample++){
      if(min > v_waveform.at(sample)){
         min = v_waveform.at(sample);
      }
   }
   return min;
};

double  getMean(vector<double> v_data){
   int n_datasize = v_data.size();
   double mean = 0;
   double num;
   for(int data=0; data<n_datasize; data++){
      num = v_data.at(data);
      mean += num;
   }
   mean /= (double)n_datasize;
   return mean;
};

double  getStdDev(vector<double> v_data){
   int n_datasize = v_data.size();
   double mean = 0;
   double meanPow2 = 0;
   double stdDev = 0;
   double num;
   for(int data=0; data<n_datasize; data++){
      num = v_data.at(data);
      mean += num;
      meanPow2 += num*num;
   }
   mean /= (double)n_datasize;
   meanPow2 /= (double)n_datasize;
   stdDev = sqrt(meanPow2 - mean*mean);
   return stdDev;

};


int getMaxSampleNum(vector<double> v_waveform){
   int maxSample = -100;
   double maxValue = -99999;
   for(int sample=0; sample<n_sample; sample++){
      if(v_waveform.at(sample) > maxValue){
         maxSample = sample;
         maxValue = v_waveform.at(sample);
      }//if
   }//for sample
   return maxSample;
};


double  getWaveHeight(vector<double> v_waveform){
   double min = 99999;
   double max = -99999;
   double waveHeight;
   for(int sample=0; sample<v_waveform.size(); sample++){
      if(min > v_waveform.at(sample)){
         min = v_waveform.at(sample);
      }//if min
      if(max < v_waveform.at(sample)){
         max = v_waveform.at(sample);
      }//if max
   }//for sample
   waveHeight = max - min;
   return waveHeight;
};

double  getTimingWeightedMean(vector<double> v_waveform, vector<double> v_timing ,int maxSample){
   double timingWeighedMean = -9999;
   int    n_data = 3;
   double timing[n_data];
   double adc[n_data];
   double adc_sum = 0;
   double weightedTotal = 0;
   if(maxSample == 0){
      cerr << "peak sample num of wave is 0" << endl;
      return 0;
   }
   for(int data=0; data<n_data; data++){
      timing[data] = (double)v_timing.at(maxSample-1+data);
      adc[data]    = (double)v_waveform.at(maxSample-1+data);
      adc_sum += adc[data];
      weightedTotal += timing[data]*adc[data];
   }//for data
   timingWeighedMean = weightedTotal / adc_sum;
   return timingWeighedMean;

};


double  getMinSampleNum(vector<double> v_waveform){
   int minSample = -100;
   double minValue = 99999;
   for(int sample=0; sample<v_waveform.size(); sample++){
      if(minValue > v_waveform.at(sample)){
         minValue = v_waveform.at(sample);
         minSample = sample;
      }
   }//for sample
   return minSample;

};

bool waveThreshold(vector<double> v_waveform, _chip m_chip[n_plane][n_chip], int plane, int chip, int strip){
   double threshold_mtp = 5.0;
   vector<double> v_baseline;
   v_baseline.resize(n_sample);
   vector<double> v_stdDev;
   v_stdDev.resize(n_sample);
   vector<double> v_fluctuation;
   v_fluctuation.resize(n_sample);
   for(int sample=0; sample<n_sample; sample++){
      v_baseline.at(sample)
         = m_chip[plane][chip].m_calib_mean[strip][sample];
      v_stdDev.at(sample)
         = m_chip[plane][chip].m_calib_stdDev[strip][sample];
   }//for sample
   for(int sample=0; sample<n_sample; sample++){
      v_fluctuation.at(sample) = v_baseline.at(sample) + v_stdDev.at(sample)*threshold_mtp;
   }
   double maxSampleNum = getMaxSampleNum(v_waveform);
   double maxAdc = getMax(v_waveform);
   double threshold = getMax(v_fluctuation);

   if(maxAdc > threshold){
      return true;
   }
   else{
      return false;
   }

};

int  classifyWaveType(vector<double> v_waveform){
   int TYPE;
   double minSampleNum = getMinSampleNum(v_waveform);
   double maxSampleNum = getMaxSampleNum(v_waveform);
   double minAdcValue = getMin(v_waveform);
   for(int sample=0; sample<v_waveform.size(); sample++){
      v_waveform.at(sample) -= minAdcValue;
   }//for sample

   if(maxSampleNum == 4){
      if(maxSampleNum == minSampleNum + 3){
         TYPE = 41;
         if(v_waveform.at(minSampleNum+1) 
               > (v_waveform.at(minSampleNum)+v_waveform.at(minSampleNum+2))/2.0)
         {
            TYPE = 411;
         }
      }
      else if(maxSampleNum == minSampleNum + 2){
         if(v_waveform.at(maxSampleNum-2)<10){
            TYPE = 421;
         }
         else{
            TYPE = 42;
         }
      }
      else if(maxSampleNum == minSampleNum + 1){
         TYPE = 43;
      }
      else{
         TYPE = 41;
      }
   }//if maxSampleNum == 4
   else if(maxSampleNum == 3){
      if(maxSampleNum == minSampleNum + 2){
         TYPE = 32;
      }
      else if(maxSampleNum == minSampleNum + 3){
         if(v_waveform.at(minSampleNum)<15){
            TYPE = 331;
         }
         else{ 
            TYPE = 33;
         }
      }
      else{
         TYPE = 31;
      }
   }
   else if(maxSampleNum == 5){
      if(maxSampleNum == minSampleNum + 2){
         TYPE = 511;
      }
      else{
         TYPE = 51;
      }
   }
   else{
      TYPE = 41;
   }

   return TYPE;
};


TVector2 calcPosDecalt(double plane, double chip, double strip)
{

   int const moduleLigicalIDOrderTable[6]
      = {2, 1, 0, 3, 4, 5};

   int const TYPE[6]
      = {2, 1, 2, 1, 2, 1};

   //order in 102, 103, 104, 106, 107, 108
   double const pin1PosXZ[6][2]
      = { {158.8503, -8.0927}, {65.4744, 109.6316}, {108.3063, 116.4844},
         {-27.8569, 156.5644}, {-121.6067, 39.1093}, {-137.8413, 79.2980}
      };

   double const pin2PosXZ[6][2]
      = { {137.8126, 79.4139}, {121.5242, 39.2155}, {27.7495, 156.6169},
         {-108.4129, 116.4305}, {-65.5595, 109.5273}, {-158.8774, -8.2090}
      };
   double stripNo = chip*(double)n_strip + strip;
   int moduleLogicalIDOrder = moduleLigicalIDOrderTable[(int)plane];
   int type = TYPE[moduleLogicalIDOrder];
   //valid on analysis for SSD standalone
   //--------------------------------------
   if(type == 2){
      stripNo = n_stripOnPlane - stripNo;
   }
   //--------------------------------------


   TVector2 tv_pin1PosXZ;
   TVector2 tv_pin2PosXZ;
   tv_pin1PosXZ.SetX(pin1PosXZ[moduleLogicalIDOrder][0]);
   tv_pin1PosXZ.SetY(pin1PosXZ[moduleLogicalIDOrder][1]);
   tv_pin2PosXZ.SetX(pin2PosXZ[moduleLogicalIDOrder][0]);
   tv_pin2PosXZ.SetY(pin2PosXZ[moduleLogicalIDOrder][1]);

   TVector2 tv_pinCenter = (tv_pin1PosXZ + tv_pin2PosXZ) / 2.0;
   TVector2 tv_pin1To2 = tv_pin2PosXZ - tv_pin1PosXZ;
   TVector2 tv_pin1To2Unit = tv_pin1To2.Unit();
   TVector2 tv_vertical = tv_pin1To2Unit.Rotate(M_PI*90.0/180.0);

   TVector2 tv_sensorCenter = tv_pinCenter + 5.0*tv_vertical;

   double strip_center = ((double)n_stripOnPlane - 1.0) / 2.0;
   double strip_relative = stripNo - strip_center;

   TVector2 tv_localPosXZ = 0.080 * strip_relative * tv_pin1To2Unit;

   TVector2 tv_globalPosXZ = tv_sensorCenter + tv_localPosXZ;

   return tv_globalPosXZ;

};


