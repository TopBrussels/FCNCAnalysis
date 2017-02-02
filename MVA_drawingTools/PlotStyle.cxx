#include <iostream>

#include "PlotStyle.h"

#include "TROOT.h"

void SetPlotStyle ()
{
  static TStyle* plotStyle = 0;
  if ( plotStyle==0 ) plotStyle = PlotStyle();
  gROOT->SetStyle("PLOT");
  gROOT->ForceStyle();
}

TStyle* PlotStyle()
{
  TStyle *plotStyle = new TStyle("PLOT","Plot style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  plotStyle->SetFrameBorderMode(icol);
  plotStyle->SetFrameFillColor(icol);
  plotStyle->SetCanvasBorderMode(icol);
  plotStyle->SetCanvasColor(icol);
  plotStyle->SetPadBorderMode(icol);
  plotStyle->SetPadColor(icol);
  plotStyle->SetStatColor(icol);
  //plotStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  plotStyle->SetPaperSize(20,26);

  // set margin sizes
  plotStyle->SetPadTopMargin(0.05);
  plotStyle->SetPadRightMargin(0.2);
//  plotStyle->SetPadRightMargin(0.05);
  plotStyle->SetPadBottomMargin(0.16);
  plotStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  plotStyle->SetTitleXOffset(1.4);
  plotStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
//  Double_t tsize=0.05;
  Double_t tsize=0.04;
  plotStyle->SetTextFont(font);

  plotStyle->SetTextSize(tsize);
  plotStyle->SetLabelFont(font,"x");
  plotStyle->SetTitleFont(font,"x");
  plotStyle->SetLabelFont(font,"y");
  plotStyle->SetTitleFont(font,"y");
  plotStyle->SetLabelFont(font,"z");
  plotStyle->SetTitleFont(font,"z");
  
  plotStyle->SetLabelSize(tsize,"x");
  plotStyle->SetTitleSize(tsize,"x");
  plotStyle->SetLabelSize(tsize,"y");
  plotStyle->SetTitleSize(tsize,"y");
  plotStyle->SetLabelSize(tsize,"z");
  plotStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  plotStyle->SetMarkerStyle(20);
  plotStyle->SetMarkerSize(1.2);
  plotStyle->SetHistLineWidth(2.);
  plotStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //plotStyle->SetErrorX(0.001);
  // get rid of error bar caps
  plotStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  plotStyle->SetOptTitle(0);
  //plotStyle->SetOptStat(1111);
  plotStyle->SetOptStat(0);
  //plotStyle->SetOptFit(1111);
  plotStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  plotStyle->SetPadTickX(1);
  plotStyle->SetPadTickY(1);

  plotStyle->SetPalette(1);
   
  return plotStyle;

}

