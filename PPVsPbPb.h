// define a function of calculating the ratio between two root graphs

void PPVsPbPb(TGraphErrors *JetCrossSectionPP, TGraphErrors *JetCrossSectionPbPb, TGraphErrors *SingleHadronYieldPP, TGraphErrors *SingleHadronYieldPbPb, TGraphErrors *JetShapePP, TGraphErrors *JetShapePbPb, TGraphErrors *JetFFPP, TGraphErrors *JetFFPbPb, TGraphErrors *JetMassPP, TGraphErrors *JetMassPbPb, TH2D *CountVsDeltaEtaDeltaPhiSpectrumPP, TH2D *CountVsDeltaEtaDeltaPhiSpectrumPbPb)
{
  CustomGlobal();
  char FileName[5000]; double RMax=pow(10,0), RMin=pow(10,-7);
  TCanvas *canvasN = new TCanvas("canvasN","canvasN",800,700);
  canvasN->Divide(2,3);
  canvasN->cd(1);   
  JetCrossSectionPP->Draw("AP");Custom(JetCrossSectionPP,1); Custom(JetCrossSectionPbPb,2); 
  //Sys->Draw("p2 same");Custom(Sys,kBlack);Sys->SetFillColorAlpha(kGray,0.01);Sys->SetFillStyle(1001);
  //Stat->Draw("P same");Custom(Stat,kBlack); Stat->SetLineWidth(2);
  JetCrossSectionPbPb->Draw("same");
  JetCrossSectionPP->GetXaxis()->SetRangeUser(10,300);JetCrossSectionPP->GetYaxis()->SetRangeUser(RMin, RMax);
  JetCrossSectionPP->GetXaxis()->SetTitle("Jet #it{p}_{T}  (GeV)");JetCrossSectionPP->GetYaxis()->SetTitle("diff jet cross section");
  gPad->SetMargin(0.14,0.03,0.18,0.15);  
  gPad->SetLogy();JetCrossSectionPP->GetYaxis()->SetTitleOffset(1.0);
  
  TLegend *leg1 = new TLegend(0.18,0.2,0.52 ,0.37);
  leg1->SetTextFont(22);leg1->SetTextSize(0.06);
  //leg1->AddEntry(Sys,"CMS 2.76 TeV","lepf"); 
  leg1->AddEntry(JetCrossSectionPP, "pp","l");
  leg1->AddEntry(JetCrossSectionPbPb, "PbPb","l");
  leg1->SetFillStyle(0);  leg1->Draw();
  TPaveText *Pt1 = new TPaveText(50,RMax*1.02,300,RMax*8,"NB");
  Pt1->AddText("diff jet crossection, R=0.3, |#eta_{jet}|<2.0");
  Pt1->SetTextFont(22);Pt1->SetTextSize(0.06);Pt1->Draw();

  canvasN->cd(2);
  RMin=pow(10,-11);
  SingleHadronYieldPP->Draw("AP"); Custom(SingleHadronYieldPP,1); Custom(SingleHadronYieldPbPb,2);
  //Sys->Draw("p2 same");Custom(Sys,kBlack);Sys->SetFillColorAlpha(kGray,0.01);Sys->SetFillStyle(1001);
  //Stat->Draw("P same");Custom(Stat,kBlack); Stat->SetLineWidth(2);
  SingleHadronYieldPbPb->Draw("P same");
  SingleHadronYieldPP->GetXaxis()->SetRangeUser(0,100);SingleHadronYieldPP->GetYaxis()->SetRangeUser(RMin, RMax);
  SingleHadronYieldPP->GetXaxis()->SetTitle("#it{p}_{T}  (GeV)");SingleHadronYieldPP->GetYaxis()->SetTitle(" single hadron Yield");
  gPad->SetMargin(0.14,0.03,0.18,0.15);
  gPad->SetLogy();SingleHadronYieldPP->GetYaxis()->SetTitleOffset(1.0);

  TLegend *leg2 = new TLegend(0.56,0.5,0.82 ,0.67);
  leg2->SetTextFont(22);leg2->SetTextSize(0.06);
  //leg2->AddEntry(G5,"CMS 2.76 TeV","lepf"); 
  leg2->AddEntry(SingleHadronYieldPP, "pp","lp");
  leg2->AddEntry(SingleHadronYieldPbPb, "PbPb","lp");
  leg2->SetFillStyle(0);  leg2->Draw();
  TPaveText *Pt2 = new TPaveText(5,RMax*1.02,100,RMax*14.2,"NB");
  Pt2->AddText("Single hadron spectrum,  |#eta|<3.0");
  Pt2->SetTextFont(22);Pt2->SetTextSize(0.06);Pt2->Draw();
  
  canvasN->cd(3);
  JetShapePP->Draw("AP"); Custom(JetShapePP,1); Custom(JetShapePbPb,2);
  //Sys->Draw("p2 same");Custom(Sys,kBlack);Sys->SetFillColorAlpha(kGray,0.01);Sys->SetFillStyle(1001);
  //Stat->Draw("P same");Custom(Stat,kBlack); Stat->SetLineWidth(2);
  JetShapePbPb->Draw("P same");
  JetShapePP->GetXaxis()->SetRangeUser(0,0.35);JetShapePP->GetYaxis()->SetRangeUser(0.01,20);
  JetShapePP->GetXaxis()->SetTitle("r");JetShapePP->GetYaxis()->SetTitle("Jet shape #rho(r)");
  gPad->SetMargin(0.12,0.03,0.18,0.15);
  gPad->SetLogy();

  TLegend *leg3 = new TLegend(0.16,0.16,0.52 ,0.37);
  leg3->SetTextFont(22);leg3->SetTextSize(0.06);
  //leg3->AddEntry(G5,"CMS 2.76 TeV","lepf"); 
  leg3->AddEntry(JetShapePP, "pp","l");
  leg3->AddEntry(JetShapePbPb, "PbPb","l");
  leg3->SetFillStyle(0);  leg3->Draw();
  TPaveText *Pt3 = new TPaveText(0, 20*1.05, 0.35, 20*10.2,"NB");
  Pt3->AddText("JetShape #rho(r), p^{jet}_{T}>10 GeV, |#eta|<2.0");
  Pt3->SetTextFont(22);Pt3->SetTextSize(0.06);Pt3->Draw();

  canvasN->cd(4);
  JetFFPP->Draw("AP"); Custom(JetFFPP,1); Custom(JetFFPbPb,2);
  //Sys->Draw("p2 same");Custom(Sys,kBlack);Sys->SetFillColorAlpha(kGray,0.01);Sys->SetFillStyle(1001);
  //Stat->Draw("P same");Custom(Stat,kBlack); Stat->SetLineWidth(2);
  JetFFPbPb->Draw("Hist ][ same");
  JetFFPP->GetXaxis()->SetRangeUser(0.0,1.0);JetFFPP->GetYaxis()->SetRangeUser(0.01, 1000);
  JetFFPP->GetXaxis()->SetTitle("Z");JetFFPP->GetYaxis()->SetTitle(" Jet FF D_{z} ");
  gPad->SetMargin(0.12,0.03,0.18,0.15);
  gPad->SetLogy();

  TLegend *leg4 = new TLegend(0.16,0.16,0.52 ,0.37);
  leg4->SetTextFont(22);leg4->SetTextSize(0.06);
  //leg4->AddEntry(G5,"CMS 2.76 TeV","lepf"); 
  leg4->AddEntry(JetFFPP, "PP","l");
  leg4->AddEntry(JetFFPbPb, "PbPb","l");
  leg4->SetFillStyle(0);  leg4->Draw();
  TPaveText *Pt4 = new TPaveText(0,1000*2.01,1,1000*8,"NB");
  Pt4->AddText("Jet Fragmentation function, |#eta|<2.0");
  Pt4->SetTextFont(22);Pt4->SetTextSize(0.06);Pt4->Draw();

  canvasN->cd(5);
  JetMassPbPb->Draw("AP"); Custom(JetMassPbPb,2);
  JetMassPP->Draw("P same"); Custom(JetMassPP,1);
  JetMassPbPb->GetXaxis()->SetRangeUser(0,30);JetMassPbPb->GetYaxis()->SetRangeUser(-0.05, 0.15);
  JetMassPbPb->GetXaxis()->SetTitle("Jet Mass (GeV)");JetMassPbPb->GetYaxis()->SetTitle("1/N dN/dM_{jet} (GeV^{-1})");
 //MCJetMassRatio->Draw("AXIS"); Custom(MCJetMassRatio,2);
  //Sys->Draw("p2 same");Custom(Sys,kBlack);Sys->SetFillColorAlpha(kGray,0.01);Sys->SetFillStyle(1001);
  //Stat->Draw("P same");Custom(Stat,kBlack); Stat->SetLineWidth(2);
  //MCJetMassRatio->Draw("Hist E X0 ][  same");MCJetMassRatio->Draw("Hist ][ same");
  //MCJetMassRatio->GetXaxis()->SetRangeUser(50,300);MCJetMassRatio->GetYaxis()->SetRangeUser(0.0, 1.5);
  //MCJetMassRatio->GetXaxis()->SetTitle("Jet Mass");MCJetMassRatio->GetYaxis()->SetTitle("Jet R_{AA}");
  gPad->SetMargin(0.12,0.03,0.15,0.15);

  TLegend *leg5 = new TLegend(0.36,0.16,0.72 ,0.37);
  leg5->SetTextFont(22);leg5->SetTextSize(0.06);
  //leg5->AddEntry(G5,"CMS 2.76 TeV","lepf"); 
  leg5->AddEntry(JetMassPP, "pp","l");
  leg5->AddEntry(JetMassPbPb, "PbPb","l");
  leg5->SetFillStyle(0);  leg5->Draw();
  TPaveText *Pt5 = new TPaveText(0,0.155,30,0.19,"NB");
  Pt5->AddText("JetMass, |#eta|<2.0");
  Pt5->SetTextFont(22);Pt5->SetTextSize(0.06);Pt5->Draw();
 

  canvasN->cd(6);
  CountVsDeltaEtaDeltaPhiSpectrumPbPb->Draw("surf2"); Custom(CountVsDeltaEtaDeltaPhiSpectrumPbPb);
  //  CountVsDeltaEtaDeltaPhiSpectrumPbPb->Draw("colz"); Custom(CountVsDeltaEtaDeltaPhiSpectrumPbPb);
  CountVsDeltaEtaDeltaPhiSpectrumPbPb->GetXaxis()->SetRangeUser(-1,1);CountVsDeltaEtaDeltaPhiSpectrumPbPb->GetYaxis()->SetRangeUser(-1.0, 2.0);
  CountVsDeltaEtaDeltaPhiSpectrumPbPb->GetXaxis()->SetTitle("#Delta#eta");
  CountVsDeltaEtaDeltaPhiSpectrumPbPb->GetYaxis()->SetTitle("#Delta#phi");
  CountVsDeltaEtaDeltaPhiSpectrumPbPb->GetZaxis()->SetTitle("#frac{#DeltaN(1,2)}{#Delta#eta_{12}#Delta#phi_{12}}");
  gPad->SetMargin(0.25,0.03,0.15,0.2);
  TLegend *leg6 = new TLegend(0.36,0.6,0.72 ,0.8);
  leg6->SetTextFont(22);leg6->SetTextSize(0.06);
  //  leg6->AddEntry(JetMassPP, "pp","l");
  leg6->AddEntry(CountVsDeltaEtaDeltaPhiSpectrumPbPb, "PbPb","");
  leg6->SetFillStyle(0);  leg6->Draw();
  TView3D *view = (TView3D*) TView::CreateView(1);
  //view->Side();
  /*  TPaveText *Pt6 = new TPaveText(-10, 7000,10,11000,"NB");
  Pt6->AddText("Two-particle correlation,1.0 < p_{T}< 10.0, p_{T1}>p_{T2}, |#eta|<2.0");
  Pt6->SetTextFont(22);Pt6->SetTextSize(0.06);Pt6->Draw();
  */
  TPaveLabel *Pt6 = new TPaveLabel(-1.2,1.0,0.8,1.3,"Two-particle correlation, p_{T}#in [1,10]GeV, p_{T1}>p_{T2}, |#eta|<2.0","NB");
  Pt6->SetTextFont(22);Pt6->SetTextSize(0.6);Pt6->Draw();
  Pt6->Draw(); 
  sprintf(FileName,"PlotTwoStageHydroTest_HadronsFromMediumAndJet.pdf");
  //canvasN->SaveAs(FileName);

}
