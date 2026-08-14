void figure2()
{
//=========Macro generated from canvas: CANb/CANb
//=========  (Tue Feb 13 18:53:13 2018) by ROOT version 6.12/04
   TCanvas *CANb = new TCanvas("CANb", "CANb",0,0,1600,800);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   CANb->Range(0,0,1,1);
   CANb->SetFillColor(0);
   CANb->SetBorderMode(0);
   CANb->SetBorderSize(2);
   CANb->SetLeftMargin(0.05);
   CANb->SetRightMargin(0);
   CANb->SetTopMargin(0);
   CANb->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: CANb_1
   TPad *CANb_1 = new TPad("CANb_1", "CANb_1",0,0,0.5011905,0.9);
   CANb_1->Draw();
   CANb_1->cd();
   CANb_1->Range(-40,-0.3055556,360,0.25);
   CANb_1->SetFillColor(0);
   CANb_1->SetBorderMode(0);
   CANb_1->SetBorderSize(2);
   CANb_1->SetGridx();
   CANb_1->SetGridy();
   CANb_1->SetRightMargin(0);
   CANb_1->SetTopMargin(0);
   CANb_1->SetFrameBorderMode(0);
   CANb_1->SetFrameBorderMode(0);
   
   TH2F *emb2__1 = new TH2F("emb2__1","emb2",200,0,360,200,-0.25,0.25);
   emb2__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   emb2__1->SetLineColor(ci);
   emb2__1->GetXaxis()->SetNdivisions(-606);
   emb2__1->GetXaxis()->SetLabelFont(42);
   emb2__1->GetXaxis()->SetLabelSize(0.06);
   emb2__1->GetXaxis()->SetTitleSize(0.035);
   emb2__1->GetXaxis()->SetTitleFont(42);
   emb2__1->GetYaxis()->SetNdivisions(505);
   emb2__1->GetYaxis()->SetLabelFont(42);
   emb2__1->GetYaxis()->SetLabelSize(0.06);
   emb2__1->GetYaxis()->SetTitleSize(0.035);
   emb2__1->GetYaxis()->SetTitleOffset(0);
   emb2__1->GetYaxis()->SetTitleFont(42);
   emb2__1->GetZaxis()->SetLabelFont(42);
   emb2__1->GetZaxis()->SetLabelSize(0.035);
   emb2__1->GetZaxis()->SetTitleSize(0.035);
   emb2__1->GetZaxis()->SetTitleFont(42);
   emb2__1->Draw("");
   
   Double_t f_g_3_1_fx1001[12] = {
   19.6015,
   48.3572,
   75.0192,
   102.6721,
   135.6295,
   163.5866,
   197.0813,
   224.426,
   256.8322,
   285.7622,
   311.4636,
   339.83};
   Double_t f_g_3_1_fy1001[12] = {
   0.1643,
   0.1782,
   0.189,
   0.1346,
   0.1834,
   0.0389,
   -0.0045,
   -0.1473,
   -0.2061,
   -0.1975,
   -0.0978,
   -0.1831};
   Double_t f_g_3_1_fex1001[12] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t f_g_3_1_fey1001[12] = {
   0.0582,
   0.0233,
   0.0194,
   0.024,
   0.0307,
   0.0327,
   0.0328,
   0.0312,
   0.0228,
   0.0195,
   0.0234,
   0.067};
   TGraphErrors *gre = new TGraphErrors(12,f_g_3_1_fx1001,f_g_3_1_fy1001,f_g_3_1_fex1001,f_g_3_1_fey1001);
   gre->SetName("f_g_3_1");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);
   gre->SetLineWidth(3);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(2.5);
   
   TH1F *Graph_f_g_3_11001 = new TH1F("Graph_f_g_3_11001","Graph",100,0,371.8528);
   Graph_f_g_3_11001->SetMinimum(-0.29736);
   Graph_f_g_3_11001->SetMaximum(0.26976);
   Graph_f_g_3_11001->SetDirectory(0);
   Graph_f_g_3_11001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_f_g_3_11001->SetLineColor(ci);
   Graph_f_g_3_11001->GetXaxis()->SetLabelFont(42);
   Graph_f_g_3_11001->GetXaxis()->SetLabelSize(0.035);
   Graph_f_g_3_11001->GetXaxis()->SetTitleSize(0.035);
   Graph_f_g_3_11001->GetXaxis()->SetTitleFont(42);
   Graph_f_g_3_11001->GetYaxis()->SetLabelFont(42);
   Graph_f_g_3_11001->GetYaxis()->SetLabelSize(0.035);
   Graph_f_g_3_11001->GetYaxis()->SetTitleSize(0.035);
   Graph_f_g_3_11001->GetYaxis()->SetTitleOffset(0);
   Graph_f_g_3_11001->GetYaxis()->SetTitleFont(42);
   Graph_f_g_3_11001->GetZaxis()->SetLabelFont(42);
   Graph_f_g_3_11001->GetZaxis()->SetLabelSize(0.035);
   Graph_f_g_3_11001->GetZaxis()->SetTitleSize(0.035);
   Graph_f_g_3_11001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_f_g_3_11001);
   
   
   TF1 *fun1002 = new TF1("fun","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,371.8528, TF1::EAddToList::kNo);
   fun1002->SetFillColor(19);
   fun1002->SetFillStyle(0);
   fun1002->SetLineColor(4);
   fun1002->SetLineWidth(2);
   fun1002->SetChisquare(23.41053);
   fun1002->SetNDF(10);
   fun1002->GetXaxis()->SetLabelFont(42);
   fun1002->GetXaxis()->SetLabelSize(0.035);
   fun1002->GetXaxis()->SetTitleSize(0.035);
   fun1002->GetXaxis()->SetTitleFont(42);
   fun1002->GetYaxis()->SetLabelFont(42);
   fun1002->GetYaxis()->SetLabelSize(0.035);
   fun1002->GetYaxis()->SetTitleSize(0.035);
   fun1002->GetYaxis()->SetTitleOffset(0);
   fun1002->GetYaxis()->SetTitleFont(42);
   fun1002->SetParameter(0,0.1929328);
   fun1002->SetParError(0,0.009561767);
   fun1002->SetParLimits(0,0,0);
   fun1002->SetParameter(1,-0.03626647);
   fun1002->SetParError(1,0.1193446);
   fun1002->SetParLimits(1,0,0);
   fun1002->SetParent(gre);
   gre->GetListOfFunctions()->Add(fun1002);
   gre->Draw("p");
   
   TF1 *funmm1003 = new TF1("funmm","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,360, TF1::EAddToList::kDefault);
   funmm1003->SetFillColor(19);
   funmm1003->SetFillStyle(0);
   funmm1003->SetLineColor(4);
   funmm1003->SetLineWidth(2);
   funmm1003->SetLineStyle(2);
   funmm1003->GetXaxis()->SetLabelFont(42);
   funmm1003->GetXaxis()->SetLabelSize(0.035);
   funmm1003->GetXaxis()->SetTitleSize(0.035);
   funmm1003->GetXaxis()->SetTitleFont(42);
   funmm1003->GetYaxis()->SetLabelFont(42);
   funmm1003->GetYaxis()->SetLabelSize(0.035);
   funmm1003->GetYaxis()->SetTitleSize(0.035);
   funmm1003->GetYaxis()->SetTitleOffset(0);
   funmm1003->GetYaxis()->SetTitleFont(42);
   funmm1003->SetParameter(0,0.183371);
   funmm1003->SetParError(0,0);
   funmm1003->SetParLimits(0,0,0);
   funmm1003->SetParameter(1,-0.1556111);
   funmm1003->SetParError(1,0);
   funmm1003->SetParLimits(1,0,0);
   funmm1003->Draw("same");
   
   TF1 *funpm1004 = new TF1("funpm","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,360, TF1::EAddToList::kDefault);
   funpm1004->SetFillColor(19);
   funpm1004->SetFillStyle(0);
   funpm1004->SetLineColor(4);
   funpm1004->SetLineWidth(2);
   funpm1004->SetLineStyle(2);
   funpm1004->GetXaxis()->SetLabelFont(42);
   funpm1004->GetXaxis()->SetLabelSize(0.035);
   funpm1004->GetXaxis()->SetTitleSize(0.035);
   funpm1004->GetXaxis()->SetTitleFont(42);
   funpm1004->GetYaxis()->SetLabelFont(42);
   funpm1004->GetYaxis()->SetLabelSize(0.035);
   funpm1004->GetYaxis()->SetTitleSize(0.035);
   funpm1004->GetYaxis()->SetTitleOffset(0);
   funpm1004->GetYaxis()->SetTitleFont(42);
   funpm1004->SetParameter(0,0.2024945);
   funpm1004->SetParError(0,0);
   funpm1004->SetParLimits(0,0,0);
   funpm1004->SetParameter(1,-0.1556111);
   funpm1004->SetParError(1,0);
   funpm1004->SetParLimits(1,0,0);
   funpm1004->Draw("same");
   
   TF1 *funmp1005 = new TF1("funmp","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,360, TF1::EAddToList::kDefault);
   funmp1005->SetFillColor(19);
   funmp1005->SetFillStyle(0);
   funmp1005->SetLineColor(4);
   funmp1005->SetLineWidth(2);
   funmp1005->SetLineStyle(2);
   funmp1005->GetXaxis()->SetLabelFont(42);
   funmp1005->GetXaxis()->SetLabelSize(0.035);
   funmp1005->GetXaxis()->SetTitleSize(0.035);
   funmp1005->GetXaxis()->SetTitleFont(42);
   funmp1005->GetYaxis()->SetLabelFont(42);
   funmp1005->GetYaxis()->SetLabelSize(0.035);
   funmp1005->GetYaxis()->SetTitleSize(0.035);
   funmp1005->GetYaxis()->SetTitleOffset(0);
   funmp1005->GetYaxis()->SetTitleFont(42);
   funmp1005->SetParameter(0,0.183371);
   funmp1005->SetParError(0,0);
   funmp1005->SetParLimits(0,0,0);
   funmp1005->SetParameter(1,0.08307814);
   funmp1005->SetParError(1,0);
   funmp1005->SetParLimits(1,0,0);
   funmp1005->Draw("same");
   
   TF1 *funpp1006 = new TF1("funpp","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,360, TF1::EAddToList::kDefault);
   funpp1006->SetFillColor(19);
   funpp1006->SetFillStyle(0);
   funpp1006->SetLineColor(4);
   funpp1006->SetLineWidth(2);
   funpp1006->SetLineStyle(2);
   funpp1006->GetXaxis()->SetLabelFont(42);
   funpp1006->GetXaxis()->SetLabelSize(0.035);
   funpp1006->GetXaxis()->SetTitleSize(0.035);
   funpp1006->GetXaxis()->SetTitleFont(42);
   funpp1006->GetYaxis()->SetLabelFont(42);
   funpp1006->GetYaxis()->SetLabelSize(0.035);
   funpp1006->GetYaxis()->SetTitleSize(0.035);
   funpp1006->GetYaxis()->SetTitleOffset(0);
   funpp1006->GetYaxis()->SetTitleFont(42);
   funpp1006->SetParameter(0,0.2024945);
   funpp1006->SetParError(0,0);
   funpp1006->SetParLimits(0,0,0);
   funpp1006->SetParameter(1,0.08307814);
   funpp1006->SetParError(1,0);
   funpp1006->SetParLimits(1,0,0);
   funpp1006->Draw("same");
  
   TH1F *funBMK_3_11007 = new TH1F("funBMK_3_1","funBMK_3_1",100,0,360);
   funBMK_3_11007->SetName("funBMK_3_1");
   funBMK_3_11007->SetTitle("funBMK_3_1");
   funBMK_3_11007->SetBinContent(0,2.163896e-17);
   funBMK_3_11007->SetBinContent(1,0.01109651);
   funBMK_3_11007->SetBinContent(2,0.0221594);
   funBMK_3_11007->SetBinContent(3,0.03315511);
   funBMK_3_11007->SetBinContent(4,0.04405017);
   funBMK_3_11007->SetBinContent(5,0.05481123);
   funBMK_3_11007->SetBinContent(6,0.06540513);
   funBMK_3_11007->SetBinContent(7,0.07579891);
   funBMK_3_11007->SetBinContent(8,0.08595989);
   funBMK_3_11007->SetBinContent(9,0.09585564);
   funBMK_3_11007->SetBinContent(10,0.1054541);
   funBMK_3_11007->SetBinContent(11,0.1147235);
   funBMK_3_11007->SetBinContent(12,0.1236325);
   funBMK_3_11007->SetBinContent(13,0.1321501);
   funBMK_3_11007->SetBinContent(14,0.1402458);
   funBMK_3_11007->SetBinContent(15,0.1478896);
   funBMK_3_11007->SetBinContent(16,0.155052);
   funBMK_3_11007->SetBinContent(17,0.1617039);
   funBMK_3_11007->SetBinContent(18,0.1677888);
   funBMK_3_11007->SetBinContent(19,0.1721902);
   funBMK_3_11007->SetBinContent(20,0.1759028);
   funBMK_3_11007->SetBinContent(21,0.1789151);
   funBMK_3_11007->SetBinContent(22,0.1812186);
   funBMK_3_11007->SetBinContent(23,0.1828075);
   funBMK_3_11007->SetBinContent(24,0.1836783);
   funBMK_3_11007->SetBinContent(25,0.1838306);
   funBMK_3_11007->SetBinContent(26,0.1832662);
   funBMK_3_11007->SetBinContent(27,0.1819897);
   funBMK_3_11007->SetBinContent(28,0.1800078);
   funBMK_3_11007->SetBinContent(29,0.17733);
   funBMK_3_11007->SetBinContent(30,0.173968);
   funBMK_3_11007->SetBinContent(31,0.1699359);
   funBMK_3_11007->SetBinContent(32,0.16525);
   funBMK_3_11007->SetBinContent(33,0.159929);
   funBMK_3_11007->SetBinContent(34,0.1539936);
   funBMK_3_11007->SetBinContent(35,0.1474669);
   funBMK_3_11007->SetBinContent(36,0.1403739);
   funBMK_3_11007->SetBinContent(37,0.1327418);
   funBMK_3_11007->SetBinContent(38,0.1245997);
   funBMK_3_11007->SetBinContent(39,0.1159786);
   funBMK_3_11007->SetBinContent(40,0.1069113);
   funBMK_3_11007->SetBinContent(41,0.09743235);
   funBMK_3_11007->SetBinContent(42,0.087578);
   funBMK_3_11007->SetBinContent(43,0.07738589);
   funBMK_3_11007->SetBinContent(44,0.0668951);
   funBMK_3_11007->SetBinContent(45,0.0561459);
   funBMK_3_11007->SetBinContent(46,0.04517965);
   funBMK_3_11007->SetBinContent(47,0.03403865);
   funBMK_3_11007->SetBinContent(48,0.02276591);
   funBMK_3_11007->SetBinContent(49,0.01140504);
   funBMK_3_11007->SetBinContent(50,0);
   funBMK_3_11007->SetBinContent(51,-0.01140504);
   funBMK_3_11007->SetBinContent(52,-0.02276591);
   funBMK_3_11007->SetBinContent(53,-0.03403865);
   funBMK_3_11007->SetBinContent(54,-0.04517965);
   funBMK_3_11007->SetBinContent(55,-0.0561459);
   funBMK_3_11007->SetBinContent(56,-0.0668951);
   funBMK_3_11007->SetBinContent(57,-0.07738589);
   funBMK_3_11007->SetBinContent(58,-0.087578);
   funBMK_3_11007->SetBinContent(59,-0.09743235);
   funBMK_3_11007->SetBinContent(60,-0.1069113);
   funBMK_3_11007->SetBinContent(61,-0.1159786);
   funBMK_3_11007->SetBinContent(62,-0.1245997);
   funBMK_3_11007->SetBinContent(63,-0.1327418);
   funBMK_3_11007->SetBinContent(64,-0.1403739);
   funBMK_3_11007->SetBinContent(65,-0.1474669);
   funBMK_3_11007->SetBinContent(66,-0.1539936);
   funBMK_3_11007->SetBinContent(67,-0.159929);
   funBMK_3_11007->SetBinContent(68,-0.16525);
   funBMK_3_11007->SetBinContent(69,-0.1699359);
   funBMK_3_11007->SetBinContent(70,-0.173968);
   funBMK_3_11007->SetBinContent(71,-0.17733);
   funBMK_3_11007->SetBinContent(72,-0.1800078);
   funBMK_3_11007->SetBinContent(73,-0.1819897);
   funBMK_3_11007->SetBinContent(74,-0.1832662);
   funBMK_3_11007->SetBinContent(75,-0.1838306);
   funBMK_3_11007->SetBinContent(76,-0.1836783);
   funBMK_3_11007->SetBinContent(77,-0.1828075);
   funBMK_3_11007->SetBinContent(78,-0.1812186);
   funBMK_3_11007->SetBinContent(79,-0.1789151);
   funBMK_3_11007->SetBinContent(80,-0.1759028);
   funBMK_3_11007->SetBinContent(81,-0.1721902);
   funBMK_3_11007->SetBinContent(82,-0.1677888);
   funBMK_3_11007->SetBinContent(83,-0.1617039);
   funBMK_3_11007->SetBinContent(84,-0.155052);
   funBMK_3_11007->SetBinContent(85,-0.1478896);
   funBMK_3_11007->SetBinContent(86,-0.1402458);
   funBMK_3_11007->SetBinContent(87,-0.1321501);
   funBMK_3_11007->SetBinContent(88,-0.1236325);
   funBMK_3_11007->SetBinContent(89,-0.1147235);
   funBMK_3_11007->SetBinContent(90,-0.1054541);
   funBMK_3_11007->SetBinContent(91,-0.09585564);
   funBMK_3_11007->SetBinContent(92,-0.08595989);
   funBMK_3_11007->SetBinContent(93,-0.07579891);
   funBMK_3_11007->SetBinContent(94,-0.06540513);
   funBMK_3_11007->SetBinContent(95,-0.05481123);
   funBMK_3_11007->SetBinContent(96,-0.04405017);
   funBMK_3_11007->SetBinContent(97,-0.03315511);
   funBMK_3_11007->SetBinContent(98,-0.0221594);
   funBMK_3_11007->SetBinContent(99,-0.01109651);
   funBMK_3_11007->SetBinContent(100,-2.163896e-17);
   funBMK_3_11007->SetBinContent(101,0);
   funBMK_3_11007->SetBinContent(102,360);
   funBMK_3_11007->SetFillColor(19);
   funBMK_3_11007->SetFillStyle(0);
   funBMK_3_11007->SetLineColor(2);
   funBMK_3_11007->SetLineWidth(2);
   funBMK_3_11007->GetXaxis()->SetLabelFont(42);
   funBMK_3_11007->GetXaxis()->SetLabelSize(0.035);
   funBMK_3_11007->GetXaxis()->SetTitleSize(0.035);
   funBMK_3_11007->GetXaxis()->SetTitleFont(42);
   funBMK_3_11007->GetYaxis()->SetLabelFont(42);
   funBMK_3_11007->GetYaxis()->SetLabelSize(0.035);
   funBMK_3_11007->GetYaxis()->SetTitleSize(0.035);
   funBMK_3_11007->GetYaxis()->SetTitleOffset(0);
   funBMK_3_11007->GetYaxis()->SetTitleFont(42);
   funBMK_3_11007->Draw("Lsame");

   Double_t f_g_3_1_fx1003[12] = {
   19.6015,
   48.3572,
   75.0192,
   102.6721,
   135.6295,
   163.5866,
   197.0813,
   224.426,
   256.8322,
   285.7622,
   311.4636,
   339.83};
   Double_t f_g_3_1_fy1003[12] = {
   0.1643,
   0.1782,
   0.189,
   0.1346,
   0.1834,
   0.0389,
   -0.0045,
   -0.1473,
   -0.2061,
   -0.1975,
   -0.0978,
   -0.1831};
   Double_t f_g_3_1_fex1003[12] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t f_g_3_1_fey1003[12] = {
   0.0582,
   0.0233,
   0.0194,
   0.024,
   0.0307,
   0.0327,
   0.0328,
   0.0312,
   0.0228,
   0.0195,
   0.0234,
   0.067};
   gre = new TGraphErrors(12,f_g_3_1_fx1003,f_g_3_1_fy1003,f_g_3_1_fex1003,f_g_3_1_fey1003);
   gre->SetName("f_g_3_1");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);
   gre->SetLineWidth(3);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(2.5);
   
   TH1F *Graph_Graph_f_g_3_110011003 = new TH1F("Graph_Graph_f_g_3_110011003","Graph",100,0,371.8528);
   Graph_Graph_f_g_3_110011003->SetMinimum(-0.29736);
   Graph_Graph_f_g_3_110011003->SetMaximum(0.26976);
   Graph_Graph_f_g_3_110011003->SetDirectory(0);
   Graph_Graph_f_g_3_110011003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_f_g_3_110011003->SetLineColor(ci);
   Graph_Graph_f_g_3_110011003->GetXaxis()->SetLabelFont(42);
   Graph_Graph_f_g_3_110011003->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_f_g_3_110011003->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_f_g_3_110011003->GetXaxis()->SetTitleFont(42);
   Graph_Graph_f_g_3_110011003->GetYaxis()->SetLabelFont(42);
   Graph_Graph_f_g_3_110011003->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_f_g_3_110011003->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_f_g_3_110011003->GetYaxis()->SetTitleOffset(0);
   Graph_Graph_f_g_3_110011003->GetYaxis()->SetTitleFont(42);
   Graph_Graph_f_g_3_110011003->GetZaxis()->SetLabelFont(42);
   Graph_Graph_f_g_3_110011003->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_f_g_3_110011003->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_f_g_3_110011003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_f_g_3_110011003);
   
   
   TF1 *fun1004 = new TF1("fun","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,371.8528, TF1::EAddToList::kNo);
   fun1004->SetFillColor(19);
   fun1004->SetFillStyle(0);
   fun1004->SetLineColor(4);
   fun1004->SetLineWidth(2);
   fun1004->SetChisquare(23.41053);
   fun1004->SetNDF(10);
   fun1004->GetXaxis()->SetLabelFont(42);
   fun1004->GetXaxis()->SetLabelSize(0.035);
   fun1004->GetXaxis()->SetTitleSize(0.035);
   fun1004->GetXaxis()->SetTitleFont(42);
   fun1004->GetYaxis()->SetLabelFont(42);
   fun1004->GetYaxis()->SetLabelSize(0.035);
   fun1004->GetYaxis()->SetTitleSize(0.035);
   fun1004->GetYaxis()->SetTitleOffset(0);
   fun1004->GetYaxis()->SetTitleFont(42);
   fun1004->SetParameter(0,0.1929328);
   fun1004->SetParError(0,0.009561767);
   fun1004->SetParLimits(0,0,0);
   fun1004->SetParameter(1,-0.03626647);
   fun1004->SetParError(1,0.1193446);
   fun1004->SetParLimits(1,0,0);
   fun1004->SetParent(gre);
   gre->GetListOfFunctions()->Add(fun1004);
   gre->Draw("p");
   funBMK_3_11007->Draw("Lsame");
   CANb_1->Modified();
   CANb->cd();
  
// ------------>Primitives in pad: CANb_2
   TPad *CANb_2 = new TPad("CANb_2", "CANb_2",0.5011905,0,0.952381,0.9);
   CANb_2->Draw();
   CANb_2->cd();
   CANb_2->Range(0,-0.3055556,360,0.25);
   CANb_2->SetFillColor(0);
   CANb_2->SetBorderMode(0);
   CANb_2->SetBorderSize(2);
   CANb_2->SetGridx();
   CANb_2->SetGridy();
   CANb_2->SetLeftMargin(0);
   CANb_2->SetRightMargin(0);
   CANb_2->SetTopMargin(0);
   CANb_2->SetFrameBorderMode(0);
   CANb_2->SetFrameBorderMode(0);
   
   TH2F *emb2__2 = new TH2F("emb2__2","emb2",200,0,360,200,-0.25,0.25);
   emb2__2->SetStats(0);

   ci = TColor::GetColor("#000099");
   emb2__2->SetLineColor(ci);
   emb2__2->GetXaxis()->SetNdivisions(-606);
   emb2__2->GetXaxis()->SetLabelFont(42);
   emb2__2->GetXaxis()->SetLabelSize(0.06);
   emb2__2->GetXaxis()->SetTitleSize(0.035);
   emb2__2->GetXaxis()->SetTitleFont(42);
   emb2__2->GetYaxis()->SetNdivisions(505);
   emb2__2->GetYaxis()->SetLabelFont(42);
   emb2__2->GetYaxis()->SetLabelSize(0.06);
   emb2__2->GetYaxis()->SetTitleSize(0.035);
   emb2__2->GetYaxis()->SetTitleOffset(0);
   emb2__2->GetYaxis()->SetTitleFont(42);
   emb2__2->GetZaxis()->SetLabelFont(42);
   emb2__2->GetZaxis()->SetLabelSize(0.035);
   emb2__2->GetZaxis()->SetTitleSize(0.035);
   emb2__2->GetZaxis()->SetTitleFont(42);
   emb2__2->Draw("");
   
   Double_t f_g_3_2_fx1005[12] = {
   28.119,
   50.2425,
   73.7667,
   103.4781,
   134.5666,
   164.3223,
   197.3449,
   224.7145,
   257.1572,
   285.6241,
   309.35,
   336.3834};
   Double_t f_g_3_2_fy1005[12] = {
   -0.3824,
   0.162,
   0.1763,
   0.1786,
   0.1254,
   0.0811,
   -0.1069,
   -0.1401,
   -0.2071,
   -0.1902,
   -0.1448,
   -0.6382};
   Double_t f_g_3_2_fex1005[12] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t f_g_3_2_fey1005[12] = {
   0.6852,
   0.0274,
   0.0258,
   0.0317,
   0.0424,
   0.0402,
   0.0427,
   0.0377,
   0.0324,
   0.0246,
   0.0322,
   0.3839};
   gre = new TGraphErrors(12,f_g_3_2_fx1005,f_g_3_2_fy1005,f_g_3_2_fex1005,f_g_3_2_fey1005);
   gre->SetName("f_g_3_2");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);
   gre->SetLineWidth(3);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(2.5);
   
   TH1F *Graph_f_g_3_21005 = new TH1F("Graph_f_g_3_21005","Graph",100,0,367.2098);
   Graph_f_g_3_21005->SetMinimum(-1.20464);
   Graph_f_g_3_21005->SetMaximum(0.43984);
   Graph_f_g_3_21005->SetDirectory(0);
   Graph_f_g_3_21005->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_f_g_3_21005->SetLineColor(ci);
   Graph_f_g_3_21005->GetXaxis()->SetLabelFont(42);
   Graph_f_g_3_21005->GetXaxis()->SetLabelSize(0.035);
   Graph_f_g_3_21005->GetXaxis()->SetTitleSize(0.035);
   Graph_f_g_3_21005->GetXaxis()->SetTitleFont(42);
   Graph_f_g_3_21005->GetYaxis()->SetLabelFont(42);
   Graph_f_g_3_21005->GetYaxis()->SetLabelSize(0.035);
   Graph_f_g_3_21005->GetYaxis()->SetTitleSize(0.035);
   Graph_f_g_3_21005->GetYaxis()->SetTitleOffset(0);
   Graph_f_g_3_21005->GetYaxis()->SetTitleFont(42);
   Graph_f_g_3_21005->GetZaxis()->SetLabelFont(42);
   Graph_f_g_3_21005->GetZaxis()->SetLabelSize(0.035);
   Graph_f_g_3_21005->GetZaxis()->SetTitleSize(0.035);
   Graph_f_g_3_21005->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_f_g_3_21005);
   
   
   TF1 *fun1006 = new TF1("fun","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,367.2098, TF1::EAddToList::kNo);
   fun1006->SetFillColor(19);
   fun1006->SetFillStyle(0);
   fun1006->SetLineColor(4);
   fun1006->SetLineWidth(2);
   fun1006->SetChisquare(5.242303);
   fun1006->SetNDF(10);
   fun1006->GetXaxis()->SetLabelFont(42);
   fun1006->GetXaxis()->SetLabelSize(0.035);
   fun1006->GetXaxis()->SetTitleSize(0.035);
   fun1006->GetXaxis()->SetTitleFont(42);
   fun1006->GetYaxis()->SetLabelFont(42);
   fun1006->GetYaxis()->SetLabelSize(0.035);
   fun1006->GetYaxis()->SetTitleSize(0.035);
   fun1006->GetYaxis()->SetTitleOffset(0);
   fun1006->GetYaxis()->SetTitleFont(42);
   fun1006->SetParameter(0,0.197381);
   fun1006->SetParError(0,0.01237221);
   fun1006->SetParLimits(0,0,0);
   fun1006->SetParameter(1,0.04146112);
   fun1006->SetParError(1,0.1523257);
   fun1006->SetParLimits(1,0,0);
   fun1006->SetParent(gre);
   gre->GetListOfFunctions()->Add(fun1006);
   gre->Draw("p");
   
   TF1 *funmm1007 = new TF1("funmm","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,360, TF1::EAddToList::kDefault);
   funmm1007->SetFillColor(19);
   funmm1007->SetFillStyle(0);
   funmm1007->SetLineColor(4);
   funmm1007->SetLineWidth(2);
   funmm1007->SetLineStyle(2);
   funmm1007->GetXaxis()->SetLabelFont(42);
   funmm1007->GetXaxis()->SetLabelSize(0.035);
   funmm1007->GetXaxis()->SetTitleSize(0.035);
   funmm1007->GetXaxis()->SetTitleFont(42);
   funmm1007->GetYaxis()->SetLabelFont(42);
   funmm1007->GetYaxis()->SetLabelSize(0.035);
   funmm1007->GetYaxis()->SetTitleSize(0.035);
   funmm1007->GetYaxis()->SetTitleOffset(0);
   funmm1007->GetYaxis()->SetTitleFont(42);
   funmm1007->SetParameter(0,0.1850088);
   funmm1007->SetParError(0,0);
   funmm1007->SetParLimits(0,0,0);
   funmm1007->SetParameter(1,-0.1108645);
   funmm1007->SetParError(1,0);
   funmm1007->SetParLimits(1,0,0);
   funmm1007->Draw("same");
   
   TF1 *funpm1008 = new TF1("funpm","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,360, TF1::EAddToList::kDefault);
   funpm1008->SetFillColor(19);
   funpm1008->SetFillStyle(0);
   funpm1008->SetLineColor(4);
   funpm1008->SetLineWidth(2);
   funpm1008->SetLineStyle(2);
   funpm1008->GetXaxis()->SetLabelFont(42);
   funpm1008->GetXaxis()->SetLabelSize(0.035);
   funpm1008->GetXaxis()->SetTitleSize(0.035);
   funpm1008->GetXaxis()->SetTitleFont(42);
   funpm1008->GetYaxis()->SetLabelFont(42);
   funpm1008->GetYaxis()->SetLabelSize(0.035);
   funpm1008->GetYaxis()->SetTitleSize(0.035);
   funpm1008->GetYaxis()->SetTitleOffset(0);
   funpm1008->GetYaxis()->SetTitleFont(42);
   funpm1008->SetParameter(0,0.2097532);
   funpm1008->SetParError(0,0);
   funpm1008->SetParLimits(0,0,0);
   funpm1008->SetParameter(1,-0.1108645);
   funpm1008->SetParError(1,0);
   funpm1008->SetParLimits(1,0,0);
   funpm1008->Draw("same");
   
   TF1 *funmp1009 = new TF1("funmp","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,360, TF1::EAddToList::kDefault);
   funmp1009->SetFillColor(19);
   funmp1009->SetFillStyle(0);
   funmp1009->SetLineColor(4);
   funmp1009->SetLineWidth(2);
   funmp1009->SetLineStyle(2);
   funmp1009->GetXaxis()->SetLabelFont(42);
   funmp1009->GetXaxis()->SetLabelSize(0.035);
   funmp1009->GetXaxis()->SetTitleSize(0.035);
   funmp1009->GetXaxis()->SetTitleFont(42);
   funmp1009->GetYaxis()->SetLabelFont(42);
   funmp1009->GetYaxis()->SetLabelSize(0.035);
   funmp1009->GetYaxis()->SetTitleSize(0.035);
   funmp1009->GetYaxis()->SetTitleOffset(0);
   funmp1009->GetYaxis()->SetTitleFont(42);
   funmp1009->SetParameter(0,0.1850088);
   funmp1009->SetParError(0,0);
   funmp1009->SetParLimits(0,0,0);
   funmp1009->SetParameter(1,0.1937868);
   funmp1009->SetParError(1,0);
   funmp1009->SetParLimits(1,0,0);
   funmp1009->Draw("same");
   
   TF1 *funpp1010 = new TF1("funpp","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,360, TF1::EAddToList::kDefault);
   funpp1010->SetFillColor(19);
   funpp1010->SetFillStyle(0);
   funpp1010->SetLineColor(4);
   funpp1010->SetLineWidth(2);
   funpp1010->SetLineStyle(2);
   funpp1010->GetXaxis()->SetLabelFont(42);
   funpp1010->GetXaxis()->SetLabelSize(0.035);
   funpp1010->GetXaxis()->SetTitleSize(0.035);
   funpp1010->GetXaxis()->SetTitleFont(42);
   funpp1010->GetYaxis()->SetLabelFont(42);
   funpp1010->GetYaxis()->SetLabelSize(0.035);
   funpp1010->GetYaxis()->SetTitleSize(0.035);
   funpp1010->GetYaxis()->SetTitleOffset(0);
   funpp1010->GetYaxis()->SetTitleFont(42);
   funpp1010->SetParameter(0,0.2097532);
   funpp1010->SetParError(0,0);
   funpp1010->SetParLimits(0,0,0);
   funpp1010->SetParameter(1,0.1937868);
   funpp1010->SetParError(1,0);
   funpp1010->SetParLimits(1,0,0);
   funpp1010->Draw("same");
  
   TH1F *funBMK_3_21011 = new TH1F("funBMK_3_2","funBMK_3_21011",100,0,360);
   funBMK_3_21011->SetName("funBMK_3_2");
   funBMK_3_21011->SetTitle("funBMK_3_2");
   funBMK_3_21011->SetBinContent(0,2.228936e-17);
   funBMK_3_21011->SetBinContent(1,0.01142984);
   funBMK_3_21011->SetBinContent(2,0.02282388);
   funBMK_3_21011->SetBinContent(3,0.0341464);
   funBMK_3_21011->SetBinContent(4,0.04536176);
   funBMK_3_21011->SetBinContent(5,0.0564345);
   funBMK_3_21011->SetBinContent(6,0.06732937);
   funBMK_3_21011->SetBinContent(7,0.07801138);
   funBMK_3_21011->SetBinContent(8,0.08844588);
   funBMK_3_21011->SetBinContent(9,0.09859856);
   funBMK_3_21011->SetBinContent(10,0.1084355);
   funBMK_3_21011->SetBinContent(11,0.1179233);
   funBMK_3_21011->SetBinContent(12,0.1270291);
   funBMK_3_21011->SetBinContent(13,0.1357203);
   funBMK_3_21011->SetBinContent(14,0.1439653);
   funBMK_3_21011->SetBinContent(15,0.1517329);
   funBMK_3_21011->SetBinContent(16,0.1589927);
   funBMK_3_21011->SetBinContent(17,0.1656792);
   funBMK_3_21011->SetBinContent(18,0.1702898);
   funBMK_3_21011->SetBinContent(19,0.1741677);
   funBMK_3_21011->SetBinContent(20,0.1773067);
   funBMK_3_21011->SetBinContent(21,0.1797042);
   funBMK_3_21011->SetBinContent(22,0.1813611);
   funBMK_3_21011->SetBinContent(23,0.1822814);
   funBMK_3_21011->SetBinContent(24,0.1824723);
   funBMK_3_21011->SetBinContent(25,0.1819435);
   funBMK_3_21011->SetBinContent(26,0.1807077);
   funBMK_3_21011->SetBinContent(27,0.1787798);
   funBMK_3_21011->SetBinContent(28,0.1761768);
   funBMK_3_21011->SetBinContent(29,0.172918);
   funBMK_3_21011->SetBinContent(30,0.1690243);
   funBMK_3_21011->SetBinContent(31,0.1645182);
   funBMK_3_21011->SetBinContent(32,0.1594239);
   funBMK_3_21011->SetBinContent(33,0.1537666);
   funBMK_3_21011->SetBinContent(34,0.147573);
   funBMK_3_21011->SetBinContent(35,0.1408706);
   funBMK_3_21011->SetBinContent(36,0.1336878);
   funBMK_3_21011->SetBinContent(37,0.1260539);
   funBMK_3_21011->SetBinContent(38,0.1179989);
   funBMK_3_21011->SetBinContent(39,0.1095534);
   funBMK_3_21011->SetBinContent(40,0.1007482);
   funBMK_3_21011->SetBinContent(41,0.09161508);
   funBMK_3_21011->SetBinContent(42,0.08218578);
   funBMK_3_21011->SetBinContent(43,0.07249254);
   funBMK_3_21011->SetBinContent(44,0.06256782);
   funBMK_3_21011->SetBinContent(45,0.05244432);
   funBMK_3_21011->SetBinContent(46,0.04215491);
   funBMK_3_21011->SetBinContent(47,0.03173261);
   funBMK_3_21011->SetBinContent(48,0.02121054);
   funBMK_3_21011->SetBinContent(49,0.01062191);
   funBMK_3_21011->SetBinContent(50,0);
   funBMK_3_21011->SetBinContent(51,-0.01062191);
   funBMK_3_21011->SetBinContent(52,-0.02121054);
   funBMK_3_21011->SetBinContent(53,-0.03173261);
   funBMK_3_21011->SetBinContent(54,-0.04215491);
   funBMK_3_21011->SetBinContent(55,-0.05244432);
   funBMK_3_21011->SetBinContent(56,-0.06256782);
   funBMK_3_21011->SetBinContent(57,-0.07249254);
   funBMK_3_21011->SetBinContent(58,-0.08218578);
   funBMK_3_21011->SetBinContent(59,-0.09161508);
   funBMK_3_21011->SetBinContent(60,-0.1007482);
   funBMK_3_21011->SetBinContent(61,-0.1095534);
   funBMK_3_21011->SetBinContent(62,-0.1179989);
   funBMK_3_21011->SetBinContent(63,-0.1260539);
   funBMK_3_21011->SetBinContent(64,-0.1336878);
   funBMK_3_21011->SetBinContent(65,-0.1408706);
   funBMK_3_21011->SetBinContent(66,-0.147573);
   funBMK_3_21011->SetBinContent(67,-0.1537666);
   funBMK_3_21011->SetBinContent(68,-0.1594239);
   funBMK_3_21011->SetBinContent(69,-0.1645182);
   funBMK_3_21011->SetBinContent(70,-0.1690243);
   funBMK_3_21011->SetBinContent(71,-0.172918);
   funBMK_3_21011->SetBinContent(72,-0.1761768);
   funBMK_3_21011->SetBinContent(73,-0.1787798);
   funBMK_3_21011->SetBinContent(74,-0.1807077);
   funBMK_3_21011->SetBinContent(75,-0.1819435);
   funBMK_3_21011->SetBinContent(76,-0.1824723);
   funBMK_3_21011->SetBinContent(77,-0.1822814);
   funBMK_3_21011->SetBinContent(78,-0.1813611);
   funBMK_3_21011->SetBinContent(79,-0.1797042);
   funBMK_3_21011->SetBinContent(80,-0.1773067);
   funBMK_3_21011->SetBinContent(81,-0.1741677);
   funBMK_3_21011->SetBinContent(82,-0.1702898);
   funBMK_3_21011->SetBinContent(83,-0.1656792);
   funBMK_3_21011->SetBinContent(84,-0.1589927);
   funBMK_3_21011->SetBinContent(85,-0.1517329);
   funBMK_3_21011->SetBinContent(86,-0.1439653);
   funBMK_3_21011->SetBinContent(87,-0.1357203);
   funBMK_3_21011->SetBinContent(88,-0.1270291);
   funBMK_3_21011->SetBinContent(89,-0.1179233);
   funBMK_3_21011->SetBinContent(90,-0.1084355);
   funBMK_3_21011->SetBinContent(91,-0.09859856);
   funBMK_3_21011->SetBinContent(92,-0.08844588);
   funBMK_3_21011->SetBinContent(93,-0.07801138);
   funBMK_3_21011->SetBinContent(94,-0.06732937);
   funBMK_3_21011->SetBinContent(95,-0.0564345);
   funBMK_3_21011->SetBinContent(96,-0.04536176);
   funBMK_3_21011->SetBinContent(97,-0.0341464);
   funBMK_3_21011->SetBinContent(98,-0.02282388);
   funBMK_3_21011->SetBinContent(99,-0.01142984);
   funBMK_3_21011->SetBinContent(100,-2.228936e-17);
   funBMK_3_21011->SetBinContent(101,0);
   funBMK_3_21011->SetBinContent(102,360);
   funBMK_3_21011->SetFillColor(19);
   funBMK_3_21011->SetFillStyle(0);
   funBMK_3_21011->SetLineColor(2);
   funBMK_3_21011->SetLineWidth(2);
   funBMK_3_21011->GetXaxis()->SetLabelFont(42);
   funBMK_3_21011->GetXaxis()->SetLabelSize(0.035);
   funBMK_3_21011->GetXaxis()->SetTitleSize(0.035);
   funBMK_3_21011->GetXaxis()->SetTitleFont(42);
   funBMK_3_21011->GetYaxis()->SetLabelFont(42);
   funBMK_3_21011->GetYaxis()->SetLabelSize(0.035);
   funBMK_3_21011->GetYaxis()->SetTitleSize(0.035);
   funBMK_3_21011->GetYaxis()->SetTitleOffset(0);
   funBMK_3_21011->GetYaxis()->SetTitleFont(42);
   funBMK_3_21011->Draw("Lsame");
   
   Double_t f_g_3_2_fx1007[12] = {
   28.119,
   50.2425,
   73.7667,
   103.4781,
   134.5666,
   164.3223,
   197.3449,
   224.7145,
   257.1572,
   285.6241,
   309.35,
   336.3834};
   Double_t f_g_3_2_fy1007[12] = {
   -0.3824,
   0.162,
   0.1763,
   0.1786,
   0.1254,
   0.0811,
   -0.1069,
   -0.1401,
   -0.2071,
   -0.1902,
   -0.1448,
   -0.6382};
   Double_t f_g_3_2_fex1007[12] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t f_g_3_2_fey1007[12] = {
   0.6852,
   0.0274,
   0.0258,
   0.0317,
   0.0424,
   0.0402,
   0.0427,
   0.0377,
   0.0324,
   0.0246,
   0.0322,
   0.3839};
   gre = new TGraphErrors(12,f_g_3_2_fx1007,f_g_3_2_fy1007,f_g_3_2_fex1007,f_g_3_2_fey1007);
   gre->SetName("f_g_3_2");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);
   gre->SetLineWidth(3);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(2.5);
   
   TH1F *Graph_Graph_f_g_3_210051007 = new TH1F("Graph_Graph_f_g_3_210051007","Graph",100,0,367.2098);
   Graph_Graph_f_g_3_210051007->SetMinimum(-1.20464);
   Graph_Graph_f_g_3_210051007->SetMaximum(0.43984);
   Graph_Graph_f_g_3_210051007->SetDirectory(0);
   Graph_Graph_f_g_3_210051007->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_f_g_3_210051007->SetLineColor(ci);
   Graph_Graph_f_g_3_210051007->GetXaxis()->SetLabelFont(42);
   Graph_Graph_f_g_3_210051007->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_f_g_3_210051007->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_f_g_3_210051007->GetXaxis()->SetTitleFont(42);
   Graph_Graph_f_g_3_210051007->GetYaxis()->SetLabelFont(42);
   Graph_Graph_f_g_3_210051007->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_f_g_3_210051007->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_f_g_3_210051007->GetYaxis()->SetTitleOffset(0);
   Graph_Graph_f_g_3_210051007->GetYaxis()->SetTitleFont(42);
   Graph_Graph_f_g_3_210051007->GetZaxis()->SetLabelFont(42);
   Graph_Graph_f_g_3_210051007->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_f_g_3_210051007->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_f_g_3_210051007->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_f_g_3_210051007);
   
   
   TF1 *fun1008 = new TF1("fun","[0]*TMath::Sin(TMath::DegToRad()*x)/(1+[1]*TMath::Cos(TMath::DegToRad()*x))",0,367.2098, TF1::EAddToList::kNo);
   fun1008->SetFillColor(19);
   fun1008->SetFillStyle(0);
   fun1008->SetLineColor(4);
   fun1008->SetLineWidth(2);
   fun1008->SetChisquare(5.242303);
   fun1008->SetNDF(10);
   fun1008->GetXaxis()->SetLabelFont(42);
   fun1008->GetXaxis()->SetLabelSize(0.035);
   fun1008->GetXaxis()->SetTitleSize(0.035);
   fun1008->GetXaxis()->SetTitleFont(42);
   fun1008->GetYaxis()->SetLabelFont(42);
   fun1008->GetYaxis()->SetLabelSize(0.035);
   fun1008->GetYaxis()->SetTitleSize(0.035);
   fun1008->GetYaxis()->SetTitleOffset(0);
   fun1008->GetYaxis()->SetTitleFont(42);
   fun1008->SetParameter(0,0.197381);
   fun1008->SetParError(0,0.01237221);
   fun1008->SetParLimits(0,0,0);
   fun1008->SetParameter(1,0.04146112);
   fun1008->SetParError(1,0.1523257);
   fun1008->SetParLimits(1,0,0);
   fun1008->SetParent(gre);
   gre->GetListOfFunctions()->Add(fun1008);
   gre->Draw("p");
   CANb_2->Modified();
   CANb->cd();
   
   TPaveText *pt = new TPaveText(0.45,0,0.55,0.08,"br");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetLineWidth(0);
   TText *pt_LaTex = pt->AddText("#phi (#circ)");
   pt->Draw();
   
   pt = new TPaveText(0,0.85,0.049,0.99,"br");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetLineWidth(0);
   pt->SetTextAngle(270);
   pt->Draw();
   
   pt = new TPaveText(0.91,0,0.99,0.08,"br");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetLineWidth(0);
   pt->Draw();
   TLatex *   tex = new TLatex(0.035,0.86,"BSA");
   tex->SetTextSize(0.07);
   tex->SetTextAngle(90);
   tex->SetLineWidth(2);
   tex->Draw();
   CANb->Modified();
   CANb->cd();
   CANb->SetSelected(CANb);
   CANb->SaveAs("figure2.png");
}
