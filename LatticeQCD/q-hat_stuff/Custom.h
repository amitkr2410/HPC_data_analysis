void DrawJetscapeLogo(double x1, double y1, double x2, double y2)
{

  // TPad * pad1 = new TPad("pad1","pad1",x1,y1,x2,y2);pad1->Draw();pad1->cd();
  //TImage *img = TImage::Open("/home/amit/Dropbox/JetscapeLogo.jpg");
  //img->Draw();
}

void Custom(TH1D *h1)
{
gPad->SetTickx();
gPad->SetTicky();
// h1->SetLineColor(color);
//h1->GetYaxis()->SetTitle("dN/dE per 1GeV");h1->GetXaxis()->SetTitle("E (GeV)");
h1->GetYaxis()->CenterTitle();
h1->GetXaxis()->CenterTitle();
//h1->GetYaxis()->SetRangeUser(0.001,10.0);
//h1->GetXaxis()->SetRangeUser(0.0,60);
h1->GetXaxis()->SetNdivisions(307);
h1->GetYaxis()->SetNdivisions(307);
h1->GetYaxis()->SetTitleSize(0.05);h1->GetYaxis()->SetTitleFont(22);
h1->GetXaxis()->SetTitleSize(0.05);h1->GetXaxis()->SetTitleFont(22);
h1->GetYaxis()->SetLabelSize(0.05);h1->GetYaxis()->SetLabelFont(22);
h1->GetXaxis()->SetLabelSize(0.05);h1->GetXaxis()->SetLabelFont(22);
h1->GetYaxis()->SetTitleOffset(1.2);
h1->GetXaxis()->SetTitleOffset(1.2);
h1->GetYaxis()->SetTickSize(0.04);
h1->GetXaxis()->SetTickSize(0.04);
h1->SetMarkerSize(3);h1->SetLineWidth(3);

}

void Custom(TH1D *h1, int color)
{
gPad->SetTickx();
gPad->SetTicky();
h1->SetLineColor(color);
h1->SetMarkerColor(color);
//h1->GetYaxis()->SetTitle("dN/dE per 1GeV");h1->GetXaxis()->SetTitle("E (GeV)");
h1->GetYaxis()->CenterTitle();
h1->GetXaxis()->CenterTitle();
//h1->GetYaxis()->SetRangeUser(0.001,10.0);
//h1->GetXaxis()->SetRangeUser(0.0,60);
h1->GetXaxis()->SetNdivisions(307);
h1->GetYaxis()->SetNdivisions(307);
h1->GetYaxis()->SetTitleSize(0.07);h1->GetYaxis()->SetTitleFont(22);
h1->GetXaxis()->SetTitleSize(0.07);h1->GetXaxis()->SetTitleFont(22);
h1->GetYaxis()->SetLabelSize(0.07);h1->GetYaxis()->SetLabelFont(22);
h1->GetXaxis()->SetLabelSize(0.07);h1->GetXaxis()->SetLabelFont(22);
h1->GetYaxis()->SetTitleOffset(0.8);
h1->GetXaxis()->SetTitleOffset(1.);
h1->GetYaxis()->SetTickSize(0.04);
h1->GetXaxis()->SetTickSize(0.04);
h1->SetMarkerSize(3);h1->SetLineWidth(3);

}
void Custom(TH2D *h1)
{
 gPad->SetTickx();
 gPad->SetTicky();
  //  h1->SetLineColor(color);
  //h1->SetMarkerColor(color);
h1->GetYaxis()->CenterTitle();
h1->GetXaxis()->CenterTitle();
h1->GetXaxis()->SetNdivisions(307);
h1->GetYaxis()->SetNdivisions(307);
h1->GetZaxis()->SetNdivisions(307);
h1->GetZaxis()->SetTitleSize(0.07);h1->GetZaxis()->SetTitleFont(22);
h1->GetYaxis()->SetTitleSize(0.07);h1->GetYaxis()->SetTitleFont(22);
h1->GetXaxis()->SetTitleSize(0.07);h1->GetXaxis()->SetTitleFont(22);
h1->GetZaxis()->SetLabelSize(0.07);h1->GetZaxis()->SetLabelFont(22);
h1->GetYaxis()->SetLabelSize(0.07);h1->GetYaxis()->SetLabelFont(22);
h1->GetXaxis()->SetLabelSize(0.07);h1->GetXaxis()->SetLabelFont(22);
h1->GetZaxis()->SetTitleOffset(1.5);
h1->GetYaxis()->SetTitleOffset(1.5);
h1->GetXaxis()->SetTitleOffset(1.5);
h1->GetZaxis()->SetTickSize(0.04);
h1->GetYaxis()->SetTickSize(0.04);
h1->GetXaxis()->SetTickSize(0.04);
h1->SetMarkerSize(3);h1->SetLineWidth(3);

}

void CustomGlobal()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLineWidth(3);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetEndErrorSize(2);
}

void Custom(TH1D & h, int color)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h.SetMarkerColor(color);
  h.SetLineColor(color);
  //h.GetYaxis()->SetTitle("dN/dE per 1GeV");h.GetXaxis()->SetTitle("E (GeV)");
  h.GetYaxis()->CenterTitle();
  h.GetXaxis()->CenterTitle();
  h.GetYaxis()->SetRangeUser(0.001,10.0);
  h.GetXaxis()->SetRangeUser(0.0,60);
  h.GetXaxis()->SetNdivisions(307);
  h.GetYaxis()->SetNdivisions(307);
  h.GetYaxis()->SetTitleSize(0.05);h.GetYaxis()->SetTitleFont(22);
  h.GetXaxis()->SetTitleSize(0.05);h.GetXaxis()->SetTitleFont(22);
  h.GetYaxis()->SetLabelSize(0.05);h.GetYaxis()->SetLabelFont(22);
  h.GetXaxis()->SetLabelSize(0.05);h.GetXaxis()->SetLabelFont(22);
  h.GetYaxis()->SetTitleOffset(1.5);
  h.GetXaxis()->SetTitleOffset(1.2);
  h.GetYaxis()->SetTickSize(0.04);
  h.GetXaxis()->SetTickSize(0.04);
  h.SetMarkerSize(3);h.SetLineWidth(3);

}

void Custom(TGraph *h1)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h1->SetLineColor(kBlack);  
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetTitleSize(0.05);h1->GetYaxis()->SetTitleFont(22);
  h1->GetXaxis()->SetTitleSize(0.05);h1->GetXaxis()->SetTitleFont(22);
  h1->GetYaxis()->SetLabelSize(0.05);h1->GetYaxis()->SetLabelFont(22);
  h1->GetXaxis()->SetLabelSize(0.05);h1->GetXaxis()->SetLabelFont(22);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTickSize(0.04);
  h1->GetXaxis()->SetTickSize(0.04);
  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(1);h1->SetLineWidth(3);
  h1->SetFillColor(0);
}
void Custom(TGraph *h1, int colornum)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h1->SetLineColor(colornum);
  h1->SetMarkerColor(colornum);
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetTitleSize(0.05);h1->GetYaxis()->SetTitleFont(22);
  h1->GetXaxis()->SetTitleSize(0.05);h1->GetXaxis()->SetTitleFont(22);
  h1->GetYaxis()->SetLabelSize(0.05);h1->GetYaxis()->SetLabelFont(22);
  h1->GetXaxis()->SetLabelSize(0.05);h1->GetXaxis()->SetLabelFont(22);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTickSize(0.04);
  h1->GetXaxis()->SetTickSize(0.04);
  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(1);h1->SetLineWidth(3);
  h1->SetFillColor(0);
}


void Custom(TGraphErrors *h1)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h1->SetLineColor(kRed);
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetTitleSize(0.05);h1->GetYaxis()->SetTitleFont(22);
  h1->GetXaxis()->SetTitleSize(0.05);h1->GetXaxis()->SetTitleFont(22);
  h1->GetYaxis()->SetLabelSize(0.05);h1->GetYaxis()->SetLabelFont(22);
  h1->GetXaxis()->SetLabelSize(0.05);h1->GetXaxis()->SetLabelFont(22);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTickSize(0.04);
  h1->GetXaxis()->SetTickSize(0.04);
  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(1);h1->SetLineWidth(3);

}

void Custom(TGraphErrors *h1, int colornum)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h1->SetLineColor(colornum);
  h1->SetMarkerColor(colornum);
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetTitleSize(0.07);h1->GetYaxis()->SetTitleFont(22);
  h1->GetXaxis()->SetTitleSize(0.07);h1->GetXaxis()->SetTitleFont(22);
  h1->GetYaxis()->SetLabelSize(0.07);h1->GetYaxis()->SetLabelFont(22);
  h1->GetXaxis()->SetLabelSize(0.07);h1->GetXaxis()->SetLabelFont(22);
  h1->GetYaxis()->SetTitleOffset(1.0);
  h1->GetXaxis()->SetTitleOffset(1.0);
  h1->GetYaxis()->SetTickSize(0.05);
  h1->GetXaxis()->SetTickSize(0.05);
  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(1);h1->SetLineWidth(4);
  h1->SetFillColor(0);
}

void TLEGEND8(double x, double y,double x2, double y2, TObject *g,string s,string st, TObject *g2,string s2,string st2, TObject *g3,string s3, string st3, TObject *g4,string s4,string st4, TObject *g5,string s5,string st5, TObject *g6,string s6,string st6, TObject *g7,string s7,string st7, TObject *g8,string s8,string st8,double TEXTSIZE)
{
  TLegend *L = new TLegend(x,y,x2,y2);
  L->AddEntry(g,s.c_str(),st.c_str());
  L->AddEntry(g2,s2.c_str(),st2.c_str());
  L->AddEntry(g3,s3.c_str(),st3.c_str());
  L->AddEntry(g4,s4.c_str(),st4.c_str());
  L->AddEntry(g5,s5.c_str(),st5.c_str());
  L->AddEntry(g6,s6.c_str(),st6.c_str());
  L->AddEntry(g7,s7.c_str(),st7.c_str());
  L->AddEntry(g8,s8.c_str(),st8.c_str());
  L->SetFillStyle(0);L->SetTextFont(22);L->SetTextSize(TEXTSIZE);
  L->Draw();
}

void TLEGEND6(double x, double y,double x2, double y2, TObject *g,string s,string st, TObject *g2,string s2,string st2, TObject *g3,string s3, string st3, TObject *g4,string s4,string st4, TObject *g5,string s5,string st5, TObject *g6,string s6,string st6, double TEXTSIZE)
{
  TLegend *L = new TLegend(x,y,x2,y2);
  L->AddEntry(g,s.c_str(),st.c_str());
  L->AddEntry(g2,s2.c_str(),st2.c_str());
  L->AddEntry(g3,s3.c_str(),st3.c_str());
  L->AddEntry(g4,s4.c_str(),st4.c_str());
  L->AddEntry(g5,s5.c_str(),st5.c_str());
  L->AddEntry(g6,s6.c_str(),st6.c_str());
  L->SetFillStyle(0);L->SetTextFont(22);L->SetTextSize(TEXTSIZE);
  L->Draw();
}

void TLEGEND5(double x, double y,double x2, double y2, TObject *g,string s,string st, TObject *g2,string s2,string st2, TObject *g3,string s3, string st3, TObject *g4,string s4,string st4, TObject *g5,string s5,string st5, double TEXTSIZE)
{
  TLegend *L = new TLegend(x,y,x2,y2);
  L->AddEntry(g,s.c_str(),st.c_str());
  L->AddEntry(g2,s2.c_str(),st2.c_str());
  L->AddEntry(g3,s3.c_str(),st3.c_str());
  L->AddEntry(g4,s4.c_str(),st4.c_str());
  L->AddEntry(g5,s5.c_str(),st5.c_str());
  L->SetFillStyle(0);L->SetTextFont(22);L->SetTextSize(TEXTSIZE);
  L->Draw();
}

void TLEGEND4(double x, double y,double x2, double y2, TObject *g,string s,string st, TObject *g2,string s2,string st2, TObject *g3,string s3,string st3, TObject *g4,string s4,string st4,double TEXTSIZE)
{
  TLegend *L = new TLegend(x,y,x2,y2);
 L->AddEntry(g,s.c_str(),st.c_str());
 L->AddEntry(g2,s2.c_str(),st2.c_str());
 L->AddEntry(g3,s3.c_str(),st3.c_str());
 L->AddEntry(g4,s4.c_str(),st4.c_str());
L->SetFillStyle(0);L->SetTextFont(22);L->SetTextSize(TEXTSIZE);
L->Draw();
}
void TLEGEND3(double x, double y,double x2, double y2, TObject *g,string s,string st, TObject *g2,string s2,string st2, TObject *g3,string s3,string st3,double TEXTSIZE)
{
  TLegend *L = new TLegend(x,y,x2,y2);
  L->AddEntry(g,s.c_str(),st.c_str());
  L->AddEntry(g2,s2.c_str(),st2.c_str());
  L->AddEntry(g3,s3.c_str(),st3.c_str());
  L->SetFillStyle(0);L->SetTextFont(22);L->SetTextSize(TEXTSIZE);  if(s.length()==0){L->SetTextColor(0);}
  L->Draw();
}


void TLEGEND2(double x, double y,double x2, double y2, TObject *g,string s,string st, TObject *g2,string s2,string st2,double TEXTSIZE)
{
  TLegend *L = new TLegend(x,y,x2,y2);
  L->AddEntry(g,s.c_str(),st.c_str());
  L->AddEntry(g2,s2.c_str(),st2.c_str());
  L->SetFillStyle(0);L->SetTextFont(22);L->SetTextSize(TEXTSIZE); if(s.length()==0){L->SetTextColor(0);}
  L->Draw();
}

void TLEGEND1(double x, double y,double x2, double y2, TObject *g,string s,string st, double TEXTSIZE)
{
  TLegend *L = new TLegend(x,y,x2,y2);
  if(st.length() == 0 ){ L->AddEntry(g,s.c_str());}
  else { L->AddEntry(g,s.c_str(),st.c_str());}
  L->SetFillStyle(0);L->SetTextFont(22);L->SetTextSize(TEXTSIZE); if(s.length()==0){L->SetTextColor(0);}
  L->Draw();
}
void TPAVETEXT(double x, double y,double x2, double y2, string s, double TEXTSIZE)
{
  TPaveText *Pt1 = new TPaveText(x,y,x2,y2,"NB");
  Pt1->AddText(s.c_str());
  Pt1->SetTextAlign(12);Pt1->SetFillColor(0); Pt1->SetTextFont(22);Pt1->SetTextSize(TEXTSIZE);  Pt1->Draw();
}


TGraphErrors * ScaleTGraphErrors(TGraphErrors *g, double PreFactor)
{
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  double x[N], y[N], ex[N], ey[N];
  for(int i=0;i<N;i++)
    {
      x[i]=X[i]; ex[i]=EX[i];
      y[i]  = Y[i]*PreFactor;
      ey[i] = EY[i]*fabs(PreFactor);
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}

TGraph * ScaleTGraph(TGraph *g, double PreFactor)
{
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();
  double x[N], y[N];
  for(int i=0;i<N;i++)
    {
      x[i]  = X[i]; 
      y[i]  = Y[i]*PreFactor;
    }
  TGraph *gr = new TGraph(N,x,y);
  return gr;
}
