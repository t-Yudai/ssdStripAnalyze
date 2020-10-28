#include <iostream>
#include <vector>

const string input_filename = "../../tmpdatalink/tmprootfile/run00385.root.txt";
const string output_param_filename = "../paramTable/fitParam.dat";
TString input_rootfilename = "../../tmpdatalink/tmprootfile/run00385.root";
TString output_pdfname = "../pdf/hoge.pdf";


const int n_sample = 8;
const int n_sample_spline = 200;
const int n_plane = 6;
const int n_chip = 6;
const int PLANE = 2;
const int CHIP = 1;
const int n_parameter = 3;

bool isMaxAt3(vector<double>);
bool isMaxAt4(vector<double>);
bool isBiggerAt8Than0(vector<double>);
bool isTailLow(vector<double>);
double  getMax(vector<double>);
double  getMin(vector<double>);
double  getMean(vector<double>);
double  getStdDev(vector<double>);


struct _waveEvent{
   int entry;
   int tdc;
   int flag;
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
   
   _waveEvent(void) : flag(1) {}
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

struct _chip{
   vector<_waveEvent> c_waveEvent;
   vector<_loop>      c_loop;
};



void waveFitter(){
   _chip m_chip[6][6];
   ifstream ifs;
   ifs.open(input_filename);
   if(ifs.fail()){
      cerr << "cannot open .txt file" << endl;
      return 0;
   }
   vector< _wave > c_wave;
   int tdc;

   TTree *tree;
   TFile *file;
   file = new TFile(input_rootfilename, "read");
   tree = (TTree*)file->Get("ssd_rawdata");
   tree->SetBranchAddress("tdc", &tdc);
   int n_entry = tree->GetEntries();
   //int n_entry = 1000;
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

   for(int entry=0; entry<n_entry; entry++){
      int plane, chip, strip;
      tree->GetEntry(entry);
      _wave s_wave;

      while(true){
         ifs >> s;
         if(s[0] == '='){
            break;
         }
      }

      ifs >> entry_input;
      if(entry_input % 500 == 0){
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
         if(isMaxAt4(v_wave) && isBiggerAt8Than0(v_wave) && isTailLow(v_wave) ){
            _waveEvent s_waveEvent;
            s_waveEvent.v_waveform.resize(n_sample);
            s_waveEvent.v_timing.resize(n_sample);
            s_waveEvent.v_waveform_original.resize(n_sample);
            s_waveEvent.v_timing_original.resize(n_sample);
            for(int sample=0; sample<n_sample; sample++){
               s_waveEvent.v_waveform_original.at(sample) = v_wave.at(sample);
               s_waveEvent.v_timing_original.at(sample) = sample*25.0;
               s_waveEvent.v_waveform.at(sample)
                  =  (v_wave.at(sample) - min) * 50.0 / range;
               s_waveEvent.v_timing.at(sample) = sample*25.0 + (tdc-1000)*0.035;
            }//for sample
            s_waveEvent.entry = entry;
            s_waveEvent.tdc = tdc;
            m_chip[plane][chip].c_waveEvent.push_back(s_waveEvent);
         }//if v_wave
      }//for multi
   }//for entry




   vector< _waveEvent > c_waveEvent;
   vector< _loop > v_loop;
   v_loop.resize(100);

   for(int plane=0; plane<n_plane; plane++){
      vector<_slice> c_slice;
      _loop s_loop;
      for(int chip=0; chip<n_chip; chip++){
         cout << "plane " << plane << " chip " << chip << endl;

         for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
            m_chip[plane][chip].c_waveEvent.at(data).graph_wave
               = new TGraph(n_sample, &m_chip[plane][chip].c_waveEvent.at(data).v_timing[0], &m_chip[plane][chip].c_waveEvent.at(data).v_waveform[0]);
            m_chip[plane][chip].c_waveEvent.at(data).spline
               = new TSpline3(Form("spline3_%d", data), m_chip[plane][chip].c_waveEvent.at(data).graph_wave);
         }//for data

         for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
            m_chip[plane][chip].c_waveEvent.at(data).v_waveform_spline.resize(n_sample_spline);
            m_chip[plane][chip].c_waveEvent.at(data).v_timing_spline.resize(n_sample_spline);
            for(int sample=0; sample<n_sample_spline; sample++){
               m_chip[plane][chip].c_waveEvent.at(data).v_waveform_spline.at(sample)
                  = m_chip[plane][chip].c_waveEvent.at(data).spline->Eval(sample);
               m_chip[plane][chip].c_waveEvent.at(data).a_waveform_spline[sample]
                  = m_chip[plane][chip].c_waveEvent.at(data).spline->Eval(sample);
               m_chip[plane][chip].c_waveEvent.at(data).v_timing_spline.at(sample) = sample;
               m_chip[plane][chip].c_waveEvent.at(data).a_timing_spline[sample] = sample;
            }//for sample
            m_chip[plane][chip].c_waveEvent.at(data).graph_spline
               = new TGraph(n_sample_spline, &m_chip[plane][chip].c_waveEvent.at(data).v_timing_spline[0], &m_chip[plane][chip].c_waveEvent.at(data).v_waveform_spline[0]);
         }//for data

         for(int sample_spline=0; sample_spline<n_sample_spline; sample_spline++){
            _slice s_slice;
            s_slice.h_projection
               = new TH1D(Form("h_projection_%d_%d_%d", plane, chip, sample_spline),Form("at_%d_%d_%d#musec;adc;count",plane, chip, sample_spline),50,0,50);
            for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
               s_slice.v_histdata.push_back(m_chip[plane][chip].c_waveEvent.at(data).spline->Eval(sample_spline));
               s_slice.h_projection->Fill(m_chip[plane][chip].c_waveEvent.at(data).spline->Eval(sample_spline));
            }//for data
            s_slice.mean = getMean(s_slice.v_histdata);
            s_slice.stdDev = getStdDev(s_slice.v_histdata);

            s_loop.c_slice.push_back(s_slice);

         }//for sample_spline
         m_chip[plane][chip].c_loop.push_back(s_loop);
      }//for chip
   }//for plane

//----------------------------------------------
//EVENT SELECTION
//
   int n_loop = 2;
   vector<double> v_timing_eval{30, 130, 190};

   for(int plane=0; plane<n_plane; plane++){
      for(int chip=0; chip<n_chip; chip++){
         m_chip[plane][chip].c_loop.resize(n_loop+1);
         for(int loop=0; loop<n_loop; loop++){
            for(int t_data=0; t_data<v_timing_eval.size(); t_data++){
               double timing_eval = v_timing_eval.at(t_data);
               double th_mean
                  = m_chip[plane][chip].c_loop.at(loop).c_slice.at(timing_eval).mean;
               double th_stdDev
                  = m_chip[plane][chip].c_loop.at(loop).c_slice.at(timing_eval).stdDev;
               double th_min = th_mean - th_stdDev;
               double th_max = th_mean + th_stdDev;

               for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
                  double value
                     = m_chip[plane][chip].c_waveEvent.at(data).spline->Eval((double)timing_eval);
                  if(value < th_min || th_max < value){
                     m_chip[plane][chip].c_waveEvent.at(data).flag = 0;
                  }//if value
               }//for data
            }//for t_data
            _slice s_slice;
            for(int sample_spline=0; sample_spline<n_sample_spline; sample_spline++){
               s_slice.h_projection
                  = new TH1D(Form("h_projection l%d p%d c%d time%d", loop, plane, chip, sample_spline),
                             Form("p%d c%d %d#musec;adc;count", plane, chip, sample_spline),
                             50, 0, 50
                        );
               for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
                  if(m_chip[plane][chip].c_waveEvent.at(data).flag == 1){
                     s_slice.v_histdata.push_back(m_chip[plane][chip].c_waveEvent.at(data).spline->Eval(sample_spline));
                     s_slice.h_projection->Fill(m_chip[plane][chip].c_waveEvent.at(data).spline->Eval(sample_spline));
                  }//if chip
               }//for data
               s_slice.mean = getMean(s_slice.v_histdata);
               s_slice.stdDev = getStdDev(s_slice.v_histdata);
               m_chip[plane][chip].c_loop.at(loop+1).c_slice.push_back(s_slice);
            }//for sample_spline
         }//for loop
      }//for chip
   }//for plane

cerr << "successed reading" << endl;


   

  TCanvas *c1 = new TCanvas("c1", "My Canvas", 0, 0, 596, 804);
  c1->Divide(6,6);

  for(int plane=0; plane<n_plane; plane++){
     for(int chip=0; chip<n_chip; chip++){
        c1->cd(plane*6 + chip + 1);
        m_chip[plane][chip].c_loop.at(n_loop).h_waves_piled
           = new TH2D(Form("hist_piled_p%d_c%d", plane, chip),
                      Form("waves_piled_p%d_c%d", plane, chip),
                      200, 0, 200, 200, 0, 55);
        for(int data=0; data<m_chip[plane][chip].c_waveEvent.size(); data++){
           for(int sample=0; sample<n_sample; sample++){
              if(m_chip[plane][chip].c_waveEvent.at(data).flag==1){
                 m_chip[plane][chip].c_loop.at(n_loop).h_waves_piled
                    ->Fill(m_chip[plane][chip].c_waveEvent.at(data).v_timing.at(sample),
                          m_chip[plane][chip].c_waveEvent.at(data).v_waveform.at(sample));
              }//for sample
           }//if flag
        }//for data
        m_chip[plane][chip].c_loop.at(n_loop).fit_func
           = new TF1(Form("fit_l%d_p%d_c%d", n_loop, plane, chip),
                     "[0]*(x-[1])/[2]*exp(-(x-[1])/[2])", 60, 160
                 );
        m_chip[plane][chip].c_loop.at(n_loop).fit_func -> SetParameter(0,50);
        m_chip[plane][chip].c_loop.at(n_loop).fit_func -> SetParameter(1,40);
        m_chip[plane][chip].c_loop.at(n_loop).fit_func -> SetParameter(2, 100);
        //m_chip[plane][chip].c_loop.at(n_loop).h_waves_piled -> Draw("colz");
        m_chip[plane][chip].c_loop.at(n_loop).h_waves_piled
           -> Fit(Form("fit_l%d_p%d_c%d", n_loop, plane, chip), "", "", 60, 160);
        m_chip[plane][chip].c_loop.at(n_loop).fit_func
           -> GetParameters(m_chip[plane][chip].c_loop.at(n_loop).fit_parameter);
     }//for chip
  }//for plane


  ofstream ofs(output_param_filename);
  ofs << "plane chip p[0] p[1] p[2]" << endl;
  for(int plane=0; plane<n_plane; plane++){
     for(int chip=0; chip<n_chip; chip++){
        ofs << plane << " " << chip << " " ;
        cout << plane << " " << chip << " " ;
        double param[n_parameter];
        for(int data=0; data<n_parameter; data++){
           param[data] = m_chip[plane][chip].c_loop.at(n_loop).fit_parameter[data];
        }
        for(int data=0; data<n_parameter; data++){
           cout << param[data] << " " ;
           ofs << param[data] << " " ;
        }//for data
        cout << endl;
        ofs << endl;
     }//for chip
  }//for plane

/* 
  c1->Divide(3,7);
  for(int slice=0; slice<20; slice++){
     c1->cd(slice+1);
     v_loop.at(n_loop).c_slice.at(10*slice).h_projection->Draw();
  }
*/

 /* 
  c1->Divide(8,8);
  c1->cd(1);
  int n_canvas = 64;
  int canvas = 0;
  int waveNum = 0;
  while(canvas<n_canvas && waveNum < c_waveEvent.size()){
     if(c_waveEvent.at(waveNum).flag == 1){
        c_waveEvent.at(waveNum).graph_wave->SetMarkerStyle(21);
        c_waveEvent.at(waveNum).graph_wave->SetMarkerSize(0.5);
        c_waveEvent.at(waveNum).graph_wave->SetMarkerColor(2);
        c_waveEvent.at(waveNum).graph_wave->Draw("AP");
        c_waveEvent.at(waveNum).graph_spline->SetMarkerStyle(21);
        c_waveEvent.at(waveNum).graph_spline->SetMarkerSize(0.5);
        c_waveEvent.at(waveNum).graph_spline->SetMarkerColor(3);
        c_waveEvent.at(waveNum).graph_spline->Draw("same");
        canvas++;
        c1->cd(canvas+1);
     }
     waveNum++;
  }//for canvas
*/
/*
  v_loop.at(n_loop).h_waves_piled
     = new TH2D(Form("hist_piled_%d",1), "waves_piled_loop2", 200, 0, 200, 200, 0, 55);

  for(int data=0; data<c_waveEvent.size(); data++){
     if(c_waveEvent.at(data).flag==1){
        for(int sample=0; sample<n_sample; sample++){
           v_loop.at(n_loop).h_waves_piled
              ->Fill(c_waveEvent.at(data).v_timing.at(sample), 
                    c_waveEvent.at(data).v_waveform.at(sample));
        }//for sample
     }//if flag
  }

  v_loop.at(n_loop).h_waves_piled->Draw();

  TF1 *func = new TF1("fit", "pol3", 40, 200);
  
  TF1 *func1 = new TF1("fit2", "[0]*(x-[1])/[2]*exp(-(x-[1])/[2])", 60, 160);
  func1 -> SetParameter(0,50);
  func1 -> SetParameter(1,40);
  func1 -> SetParameter(2,100);
 
  v_loop.at(n_loop).h_waves_piled->Fit("fit2", "", "", 60, 160);



 
  c1->SaveAs(output_pdfname+"[", "pdf");


  v_loop.at(n_loop).h_waves_piled->Draw();
  v_loop.at(n_loop).h_waves_piled->Fit("fit2", "", "", 60, 160);
  c1->SaveAs(output_pdfname, "pdf");

  c1->Clear();
  c1->Divide(3,7);
  for(int slice=0; slice<20; slice++){
     c1->cd(slice+1);
     v_loop.at(n_loop).c_slice.at(10*slice).h_projection->Draw();
  }
  c1->SaveAs(output_pdfname, "pdf");


  c1->Clear();
  c1->Divide(8,8);
  c1->cd(1);
  int n_canvas = 64;
  int canvas = 0;
  int waveNum = 0;
  while(canvas<n_canvas && waveNum < c_waveEvent.size()){
     if(c_waveEvent.at(waveNum).flag == 1){
        c_waveEvent.at(waveNum).graph_wave->SetMarkerStyle(21);
        c_waveEvent.at(waveNum).graph_wave->SetMarkerSize(0.5);
        c_waveEvent.at(waveNum).graph_wave->SetMarkerColor(2);
        c_waveEvent.at(waveNum).graph_wave->Draw("AP");
        c_waveEvent.at(waveNum).graph_spline->SetMarkerStyle(21);
        c_waveEvent.at(waveNum).graph_spline->SetMarkerSize(0.5);
        c_waveEvent.at(waveNum).graph_spline->SetMarkerColor(3);
        c_waveEvent.at(waveNum).graph_spline->Draw("same");
        canvas++;
        c1->cd(canvas+1);
     }
     waveNum++;
  }//for canvas
  c1->SaveAs(output_pdfname, "pdf");

  c1->SaveAs(output_pdfname+"]", "pdf");
*/
}//waveDrawer2





bool isMaxAt3(vector<double> v_waveform){
   int max_sample=0;
   int max_adc=0;
   for(int sample=0; sample<v_waveform.size(); sample++){
      if(max_adc < v_waveform.at(sample)){
         max_sample = sample;
         max_adc = v_waveform.at(sample);
      }
   }
   if(max_sample == 2){
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
   if(max_sample == 3){
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








