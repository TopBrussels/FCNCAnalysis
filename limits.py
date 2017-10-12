import numpy as np
from ROOT import TH2F, TCanvas, gStyle, TLatex, TAxis, TLine, TGraphErrors, TGraphAsymmErrors, TLegend, kGreen, kYellow

nchannels = 7
fontsize = 0.04

def limitsZut():

    # put real limits here: lepton+jets, dilepton, combined
                             # 3mu;    1e2mu,        2e1mu,       ,3e,       st ,          TT,    combined
    obs       = np.array( [0.00053, 0.00064, 0.00056, 0.00038, 0.00025, 0.00039, 0.00024] )

    upper2sig = np.array( [0.00074, 0.00074, 0.00084, 0.0013,  0.00049, 0.00054, 0.00033]   )
    upper1sig = np.array( [0.00050, 0.00050, 0.00056, 0.00082, 0.00032, 0.00038, 0.00023] )
    expect    = np.array( [0.00032, 0.00032, 0.00036, 0.00050, 0.00020, 0.00025, 0.00015]                )
    lower1sig = np.array( [0.00021, 0.00021, 0.00023, 0.00031, 0.00013, 0.00017, 0.000097]  )
    lower2sig = np.array( [0.00014, 0.00014, 0.00016, 0.00021, 0.000086, 0.00012, 0.000068]   )

    for ichan in range( nchannels ):
         upper2sig[ichan] = upper2sig[ichan] - expect[ichan]
	 upper1sig[ichan] = upper1sig[ichan] - expect[ichan]
	 lower1sig[ichan] = expect[ichan] - lower1sig[ichan] 
	 lower2sig[ichan] = expect[ichan] - lower2sig[ichan]
	 
    sig_inj   = np.array( [ 0.00027,0.00027,0.00027,0.00027,0.00027,0.00027,0.00027] ) 

    channels = np.array( [19.5,16.5,13.5, 10.5 , 7.5, 4.5, 1.5 ] )
    ey      = np.array( [ 0.7, 0.7, 0.7 , 0.7, 0.7, 0.7, 0.7 ] )
    zero    = np.zeros( nchannels )

    xmin = 0.00005
    xmax = 0.01 

    c,h = draw_canvas_histo( xmin, xmax, "95% CL limit on B(tZu)" )
    c.SetLogx()

    gexpect1sig = TGraphAsymmErrors( nchannels, expect, channels, lower1sig, upper1sig, ey, ey )
    gexpect1sig.SetFillColor( kGreen+2 )
    gexpect1sig.SetLineWidth( 2 )
    gexpect1sig.SetLineStyle( 2 )
    
    gexpect2sig = TGraphAsymmErrors( nchannels, expect, channels, lower2sig, upper2sig, ey, ey )
    gexpect2sig.SetFillColor( kYellow-4 )
    gexpect2sig.SetLineWidth( 2 )
    gexpect2sig.SetLineStyle( 2 )

    gexpect2sig.Draw("2")
    gexpect1sig.Draw("2")

    gobs = TGraphErrors( nchannels, obs, channels, zero, ey )
    gobs.SetMarkerStyle( 21 )
    gobs.SetMarkerSize( 0.5 )
    gobs.SetLineWidth( 2 )
    gobs.Draw("pz")

    gsinj = TGraphErrors( nchannels, sig_inj, channels, zero, ey )
    gsinj.SetLineWidth( 2 )
    gsinj.SetLineColor( 4 )
    gsinj.SetLineStyle( 2 )
    gsinj.Draw("z")

    # dashed line at median expected limits
    l = TLine()
    l.SetLineStyle( 2 )
    l.SetLineWidth( 2 )
    for bin in range( nchannels ):
        l.DrawLine( expect[bin], channels[bin]-ey[bin], expect[bin], channels[bin]+ey[bin] )

    # legend
    x2 = 1-gStyle.GetPadRightMargin()
    y1 = gStyle.GetPadBottomMargin()+gStyle.GetTickLength()+0.02
    y2 = 1-gStyle.GetPadTopMargin()
    #leg = TLegend( x2-0.3, y2-0.21, x2, y2 )
    leg = TLegend( x2-0.34, y1, x2-0.002, y1+0.18 )
    leg.SetFillColor( 10 )
    leg.SetFillStyle(1000)
    leg.AddEntry( gexpect1sig, "Expected #pm1#sigma", "FL" )
    leg.AddEntry( gexpect2sig, "Expected #pm2#sigma", "FL" )
    leg.AddEntry( gsinj,       "exp. JHEP07(2017)003", "L" )
    leg.AddEntry( gobs,        "Observed", "pl" )
    leg.Draw()

    #draw_disclaimer()

    c.RedrawAxis()    
    c.Modified()
    c.Update()
    c.SaveAs( "TOP-17-017_limitsZut.pdf" )


def limitsZct():

    # put real limits here: lepton+jets, dilepton, combined
                             # 3mu;    1e2mu,        2e1mu,       ,3e,       st ,          TT,    combined
    obs       = np.array( [0.0012, 0.00096, 0.00099, 0.00095, 0.0017, 0.00043, 0.00045] )

    upper2sig = np.array( [0.0015, 0.0018, 0.0020, 0.0029,  0.0023, 0.00094, 0.00081]   )
    upper1sig = np.array( [0.0010, 0.0012, 0.0014, 0.0019, 0.0016, 0.00066, 0.00056] )
    expect    = np.array( [0.00070, 0.00079, 0.00089, 0.0012, 0.0010, 0.00044, 0.00037]                )
    lower1sig = np.array( [0.00046, 0.00052, 0.00058, 0.00075, 0.00066, 0.00029, 0.00025]  )
    lower2sig = np.array( [0.00032, 0.00036, 0.00040, 0.00050, 0.00045, 0.00020, 0.00017]   )

    for ichan in range( nchannels ):
         upper2sig[ichan] = upper2sig[ichan] - expect[ichan]
	 upper1sig[ichan] = upper1sig[ichan] - expect[ichan]
	 lower1sig[ichan] = expect[ichan] - lower1sig[ichan] 
	 lower2sig[ichan] = expect[ichan] - lower2sig[ichan]
	 
    sig_inj   = np.array( [ 0.00118,  0.00118, 0.00118,0.00118,0.00118,0.00118,0.00118] ) 

    channels = np.array( [19.5,16.5,13.5, 10.5 , 7.5, 4.5, 1.5 ] )
    ey      = np.array( [ 0.7, 0.7, 0.7 , 0.7, 0.7, 0.7, 0.7 ] )
    zero    = np.zeros( nchannels )

    xmin = 0.00005
    xmax = 0.1 

    c,h = draw_canvas_histo( xmin, xmax, "95% CL limit on B(tZc)" )
    c.SetLogx()

    gexpect1sig = TGraphAsymmErrors( nchannels, expect, channels, lower1sig, upper1sig, ey, ey )
    gexpect1sig.SetFillColor( kGreen+2 )
    gexpect1sig.SetLineWidth( 2 )
    gexpect1sig.SetLineStyle( 2 )
    
    gexpect2sig = TGraphAsymmErrors( nchannels, expect, channels, lower2sig, upper2sig, ey, ey )
    gexpect2sig.SetFillColor( kYellow-4 )
    gexpect2sig.SetLineWidth( 2 )
    gexpect2sig.SetLineStyle( 2 )

    gexpect2sig.Draw("2")
    gexpect1sig.Draw("2")

    gobs = TGraphErrors( nchannels, obs, channels, zero, ey )
    gobs.SetMarkerStyle( 21 )
    gobs.SetMarkerSize( 0.5 )
    gobs.SetLineWidth( 2 )
    gobs.Draw("pz")

    gsinj = TGraphErrors( nchannels, sig_inj, channels, zero, ey )
    gsinj.SetLineWidth( 2 )
    gsinj.SetLineColor( 4 )
    gsinj.SetLineStyle( 2 )
    gsinj.Draw("z")

    # dashed line at median expected limits
    l = TLine()
    l.SetLineStyle( 2 )
    l.SetLineWidth( 2 )
    for bin in range( nchannels ):
        l.DrawLine( expect[bin], channels[bin]-ey[bin], expect[bin], channels[bin]+ey[bin] )

    # legend
    x2 = 1-gStyle.GetPadRightMargin()
    y1 = gStyle.GetPadBottomMargin()+gStyle.GetTickLength()+0.02
    y2 = 1-gStyle.GetPadTopMargin()
    #leg = TLegend( x2-0.3, y2-0.21, x2, y2 )
    leg = TLegend( x2-0.34, y1, x2-0.002, y1+0.18 )
    leg.SetFillColor( 10 )
    leg.SetFillStyle(1000)
    leg.AddEntry( gexpect1sig, "Expected #pm1#sigma", "FL" )
    leg.AddEntry( gexpect2sig, "Expected #pm2#sigma", "FL" )
    leg.AddEntry( gsinj,       "exp. JHEP07(2017)003", "L" )
    leg.AddEntry( gobs,        "Observed", "pl" )
    leg.Draw()

    #draw_disclaimer()

    c.RedrawAxis()    
    c.Modified()
    c.Update()
    c.SaveAs( "TOP-17-017_limitsZct.pdf" )




def draw_canvas_histo( xmin, xmax, title ):
    c = TCanvas( "c", "Canvas", 600, 600 )
    c.Draw()
    
    h = TH2F( "h", "", 10, xmin, xmax, 3*nchannels+2, 0, 3*nchannels+2 )
    h.Draw()
    h.SetStats( 0 )
    h.SetXTitle( title )

    yaxis = h.GetYaxis()
    yaxis.SetLabelSize( 1.5*fontsize )
    yaxis.SetBinLabel( 20, "3#mu" )
    yaxis.SetBinLabel( 17, "1e2#mu" )
    yaxis.SetBinLabel(14, "2e1#mu" )
    yaxis.SetBinLabel(11, "3e" )
    yaxis.SetBinLabel( 8, "STSR" )
    yaxis.SetBinLabel( 5, "TTSR" )
    yaxis.SetBinLabel( 2, "Combined" )

    # separating combined result
    l = TLine()
    l.SetLineStyle( 1 )
    l.DrawLine( xmin, 3, xmax, 3 )

    pub = TLatex();
    pub.SetNDC()
    pub.SetTextFont( 42 )
    pub.SetTextSize( 0.045 )
    pub.SetTextAlign( 11 )
    pub.DrawLatex( gStyle.GetPadLeftMargin()+0.02, 0.84, "#scale[1.2]{#bf{CMS}} #it{Preliminary}" )

    lumi = TLatex();
    lumi.SetNDC()
    lumi.SetTextFont( 42 )
    lumi.SetTextSize( 0.035 )
    lumi.SetTextAlign( 31 )
    lumi.DrawLatex( 1-gStyle.GetPadRightMargin(), 0.92, "36 fb^{-1} (13 TeV)" )


    return c,h

def draw_disclaimer():
    # disclaimer
    t = TLatex();
    t.SetNDC()
    t.SetTextSize( 0.1 )
    t.SetTextAlign( 22 )
    t.SetTextAngle( 45 )
    t.DrawText( 0.5, 0.5, "FAKE VALUES" )
    
def my_style():
    
    gStyle.SetLabelSize( fontsize, "x" );
    gStyle.SetLabelSize( fontsize, "y" );
    gStyle.SetLabelSize( fontsize, "z" );

    gStyle.SetTitleFontSize( 1.5*fontsize );
    gStyle.SetTitleSize( fontsize, "x" );
    gStyle.SetTitleSize( fontsize, "y" );
    gStyle.SetTitleSize( fontsize, "z" );


    gStyle.SetTitleOffset( 1.5, "xy" );
    gStyle.SetTitleFont( 62, "bla" );

    gStyle.SetPadBottomMargin(0.15);
    gStyle.SetPadTopMargin(0.10);
    gStyle.SetPadLeftMargin(0.22);
    gStyle.SetPadRightMargin(0.05);

    gStyle.SetStatX( 0.88 );
    gStyle.SetStatY( 0.87 );
    gStyle.SetNdivisions( 505 );
    gStyle.SetTickLength(0,"Y");
    gStyle.SetEndErrorSize(6)

    gStyle.SetCanvasColor(-1); 
    gStyle.SetPadColor(-1); 
    gStyle.SetFrameFillColor(-1); 
    gStyle.SetTitleFillColor(-1); 
    gStyle.SetFillColor(-1); 
    gStyle.SetFillStyle(4000); 
    gStyle.SetStatStyle(0); 
    gStyle.SetTitleStyle(0); 
    gStyle.SetCanvasBorderSize(0); 
    gStyle.SetFrameBorderSize(0); 
    gStyle.SetLegendBorderSize(0); 
    gStyle.SetStatBorderSize(0); 
    gStyle.SetTitleBorderSize(0); 
    
if __name__ == '__main__':
    my_style()
    limitsZct()
    limitsZut()
    
