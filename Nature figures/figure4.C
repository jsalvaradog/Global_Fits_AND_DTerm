
double tvals[5] = {0.11,0.15,0.2,0.26,0.34};
double dvals[5] = {-1.57,-1.36,-1.28,-1.10,-0.92};
double errors[5] = {0.1,0.08,0.08,0.07,0.07};

void figure4(void){
 gStyle->SetOptTitle(0);
 gStyle->SetPadTopMargin(0.01);
 gStyle->SetPadRightMargin(0.01);
 TLatex *myLatex = new TLatex();
 TCanvas *CAN = new TCanvas("CAN","CAN",800,800);
 CAN->SetGrid();
 TGraphErrors *myG = new TGraphErrors(5,tvals,dvals,0,errors);
 myG->SetMarkerStyle(21);
 myG->SetMarkerSize(2.5);
 myG->SetLineWidth(4);
 myG->GetXaxis()->SetLabelSize(0.05);
 myG->GetYaxis()->SetLabelSize(0.05);
 
 TF1 *powerlawfit = new TF1("powerlawfit","[0]*TMath::Power(1+x/([1]*[1]),-3)",0.0,1);
 powerlawfit->SetNpx(5000);
 powerlawfit->SetParNames("Constant","M");
 powerlawfit->SetLineColor(1);
 powerlawfit->SetParameters(-2.0,1.0);
 myG->Fit("powerlawfit","+");
 powerlawfit->Draw("same");

 Double_t C = powerlawfit->GetParameter(0);
 Double_t M = powerlawfit->GetParameter(1);
 Double_t alpha = -3;
 TF1 *fitErrBand0 = new TF1("fitErrBand1",Form("%f*TMath::Power(1+x/(%f*%f),%f)",C,M,M,alpha),0,1);
 fitErrBand0->SetLineColor(2);
 fitErrBand0->SetLineStyle(2);
 TF1 *fitErrBand1 = new TF1("fitErrBand1",Form("%f*TMath::Power(1+x/(%f*%f),%f)",C+0.5,M,M,alpha),0,1);
 TF1 *fitErrBandL1 = new TF1("fitErrBandL1",Form("%f*TMath::Power(1+x/(%f*%f),%f)",C+0.2,M,M,alpha),0,1);
 fitErrBand1->SetLineColor(4);
 TF1 *fitErrBand2 = new TF1("fitErrBand2",Form("%f*TMath::Power(1+x/(%f*%f),%f)",C,M+0.15,M+0.15,alpha),0,1);
 TF1 *fitErrBandL2 = new TF1("fitErrBandL2",Form("%f*TMath::Power(1+x/(%f*%f),%f)",C,M+0.1,M+0.1,alpha),0,1);
 fitErrBand2->SetLineColor(6);
 TF1 *fitErrBand3 = new TF1("fitErrBand3",Form("%f*TMath::Power(1+x/(%f*%f),%f)",C,M,M,alpha-0.5),0,1);
 fitErrBand3->SetLineColor(8);

 TF1 *fitErr12Band1 = new TF1("fitErr12Band1",Form("%f*TMath::Power(1+x/(%f*%f),%f)",C+0.25,M,M,alpha),0,1);
 TF1 *fitErr12Band2 = new TF1("fitErr12Band2",Form("%f*TMath::Power(1+x/(%f*%f),%f)",C,M+0.15,M+0.07,alpha),0,1);
 TF1 *fitErr12Band3 = new TF1("fitErr12Band3",Form("%f*TMath::Power(1+x/(%f*%f),%f)",C,M,M,alpha-0.25),0,1);

 Double_t XerrorBand_disp[1000], YerrorBand_disp[1000],  dYerrorBand_disp[1000];
 Double_t YerrorBandB_disp[1000],  dYerrorBandB_disp[1000];
 Double_t YerrorBand12_disp[1000],  dYerrorBand12_disp[1000];
 Double_t dYerrorBand1_disp[1000];
 Double_t dYerrorBand2_disp[1000];
 Double_t dYerrorBand3_disp[1000];
 for(int ip=0;ip<1000;ip++){
 	Double_t X = 0.0001 + 0.5*ip/(1.*1000);
	Double_t Y = powerlawfit->Eval(X);
	Double_t dY2 = 0;
	dY2 += TMath::Power(fitErrBandL1->Eval(X)-Y,2);
	dY2 += TMath::Power(fitErrBandL2->Eval(X)-Y,2);
	dY2 *= 0.5;
	Double_t dY1 = 0;
	dY1 +=  TMath::Power(fitErrBand1->Eval(X)-Y,2);
	dY1 +=  TMath::Power(fitErrBand2->Eval(X)-Y,2);
	dY1 +=  TMath::Power(fitErrBand3->Eval(X)-Y,2);
	dY1 *= 0.5;
	dY1 -= dY2;
	Double_t dY12 = 0;
	dY12 +=  TMath::Power(fitErr12Band1->Eval(X)-Y,2);
	dY12 +=  TMath::Power(fitErr12Band2->Eval(X)-Y,2);
	dY12 +=  TMath::Power(fitErr12Band3->Eval(X)-Y,2);
	dY12 *= 0.5;
	dY12 -= dY2;
	dY1 = TMath::Sqrt(dY1);
	dY2 = TMath::Sqrt(dY2);
	dY12 = TMath::Sqrt(dY12);
	dYerrorBand1_disp[ip] = TMath::Abs(fitErrBand1->Eval(X)-Y);
	dYerrorBand2_disp[ip] = TMath::Abs(fitErrBand2->Eval(X)-Y);
	dYerrorBand3_disp[ip] = TMath::Abs(fitErrBand3->Eval(X)-Y);
	XerrorBand_disp[ip] = X;
	YerrorBand_disp[ip] = Y;
	dYerrorBand_disp[ip] = dY2;
	YerrorBandB_disp[ip] = -2.2;
	dYerrorBandB_disp[ip] = dY1;
	YerrorBand12_disp[ip] = -2.2;
	dYerrorBand12_disp[ip] = dY12;
 }
 TGraph *gErrorBand   = new TGraph(2000);
 TGraph *gErrorBandB  = new TGraph(2000);
 TGraph *gErrorBand12 = new TGraph(2000);
 for(int ip=0;ip<1000;ip++){
 	Double_t dY = YerrorBand_disp[ip]+dYerrorBand_disp[ip];
	gErrorBand->SetPoint(ip,XerrorBand_disp[ip],dY);
	dY = YerrorBand_disp[999-ip]-dYerrorBand_disp[999-ip];
	gErrorBand->SetPoint(1000+ip,XerrorBand_disp[999-ip],dY);

	dY = YerrorBandB_disp[ip]+dYerrorBandB_disp[ip];
	gErrorBandB->SetPoint(ip,XerrorBand_disp[ip],dY);
	dY = YerrorBandB_disp[999-ip];
	gErrorBandB->SetPoint(1000+ip,XerrorBand_disp[999-ip],dY);

	dY = YerrorBand12_disp[ip]+dYerrorBand12_disp[ip];
	gErrorBand12->SetPoint(ip,XerrorBand_disp[ip],dY);
	dY = YerrorBand12_disp[999-ip];
	gErrorBand12->SetPoint(1000+ip,XerrorBand_disp[999-ip],dY);
 }

 gErrorBand->SetFillColor(1);
 gErrorBand->SetFillStyle(3004);
 gErrorBandB->SetFillColor(1);
 gErrorBandB->SetFillStyle(3005);
 gErrorBand12->SetFillColor(1);
 gErrorBand12->SetFillStyle(3001);

 powerlawfit->Draw();
 powerlawfit->GetXaxis()->SetRangeUser(0,0.5);
 powerlawfit->GetYaxis()->SetRangeUser(-2.25,-0.61);
 
 myG->GetXaxis()->SetRangeUser(0,0.5);
 myG->GetYaxis()->SetRangeUser(-1.99,-0.61);
 gErrorBandB->Draw("fsame");
 gErrorBand12->Draw("fsame");
 myG->Draw("Psame");
 myLatex->DrawLatex(0.165,-0.75,"d_{1}(t) = d_{1}#left(1 - t/M^{2}#right)^{-#alpha}");
 TPave *CoverPave1 = new TPave(0.415,-2.4,0.525,-2.26);
 CoverPave1->SetBorderSize(0);
 CoverPave1->SetFillColor(0);
 CoverPave1->SetLineWidth(0);
 CoverPave1->Draw();
 myLatex->DrawLatex(0.405,-2.37,"-t (GeV^{2})");
 TPave *CoverPave2 = new TPave(-0.04,-0.75,-0.01,-0.6);
 CoverPave2->SetBorderSize(0);
 CoverPave2->SetFillColor(0);
 CoverPave2->SetLineWidth(0);
 CoverPave2->Draw();
 myLatex->SetTextAngle(90);
 myLatex->DrawLatex(-0.02,-0.76,"d_{1}(t)");
 CAN->Print("figure4.png");
 gSystem->Exit(0);
}
