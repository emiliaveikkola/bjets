import ROOT as root
from ROOT import TCanvas, TH1D, TFile, TBrowser, TProfile, TColor, TLegend, TColor
from collections import Counter
import numpy as np
import math
#Datafile
f = root.TFile("Muo16_MC.root")
myTree = f.Get("tree")
#Canvas
c1 = TCanvas( 'c1', "", 200, 10, 700, 900 )
c1.SetLogx()
c1.SetWindowSize(1100, 700)
#Jets
rjet = TProfile( 'All b jets', '', 200, 0, 400 )
rjet_bquark = TProfile( 'B fragmentation', '', 200, 0, 400 )
rjet_bhad = TProfile( 'B production fraction', '', 200, 0, 400 )
rjet_bsl = TProfile( 'B to SL BR ', '', 200, 0, 400 )
#Functions for the weights
#B fragmentation
ppar1 = [24.75, -55.70, 47.86, -20.58, 5.085, 0.5978]
p_bquark = np.poly1d(ppar1)
#B to SL BR
b_ppar = [0.00000000000000003964, 0.000000000008662, -0.000000007341,
 0.000001999, -0.0001747 , 0.9964]
b_p = np.poly1d(b_ppar)
bt_ppar = [0.00000000000008169, -0.0000000001078, 0.00000004767,
 -8.380E-06, 0.0005008 , 0.9685]
bt_p = np.poly1d(bt_ppar)
bc_ppar = [0.0000000000002362, -0.0000000003011, 0.0000001280, -0.00002129,
 0.001291 , 1.013]
bc_p = np.poly1d(bc_ppar)
bac_ppar = [0.0000000000002264, -0.0000000003879, 0.0000002135,
 -0.00004485, 0.003478 , 1.340 ]
bac_p = np.poly1d(bc_ppar)
#Loop for the histograms
i=0
for entry in myTree:
    i += 1
    bpt1 = entry.bpt1
    gen_bpt1 = entry.gen_bpt1
    bflav1 = entry.bflav1

    bpt2 = entry.bpt2
    gen_bpt2 = entry.gen_bpt2
    bflav2 = entry.bflav2
    x1 = entry.gen_bXB1
    x2 = entry.gen_bXB2
    fitprob = entry.fitProb
    gen_bFlags1 = entry.gen_bFlags1
    gen_bFlags2 = entry.gen_bFlags2
    gen_bLeadId1 = entry.gen_bLeadId1
    gen_bLeadId2 = entry.gen_bLeadId2
    if (gen_bpt1 == 0) or (gen_bpt2 == 0):
        continue
    #if (i==10000):
        #break
    if fitprob < 0.2:
        continue
    resp1 = bpt1/gen_bpt1
    resp2 = bpt2/gen_bpt2
    had_bpt1 = x1*gen_bpt1
    had_bpt2 = x2*gen_bpt2
    #Functions for B production fraction
    p1_b0 = -0.4147*math.exp(-0.09205*x1)+0.9424
    p2_b0 = -0.4147*math.exp(-0.09205*x2)+0.9424
    p1_b = -0.4066*math.exp(-0.1038*x1)+0.9996
    p2_b = -0.4066*math.exp(-0.1038*x2)+0.9996
    p1_l = 6.943*math.exp(-0.07238*x1)+1.637
    p2_l = 6.943*math.exp(-0.07238*x2)+1.637
    #Without weights
    rjet.Fill(gen_bpt1,resp1,1)
    rjet.Fill(gen_bpt2,resp2,1)
    weight1 = 1
    weight2 = 1
    #B production fraction
    #B0 and B+
    if abs((gen_bLeadId1 == 521)) or abs((gen_bLeadId1 == 511)):
        weight1 = p1_b0
    elif abs((gen_bLeadId2 == 521)) or abs((gen_bLeadId2 == 511)):
        weight2 = p2_b0
    #B0_s and other B mesons
    elif abs((gen_bLeadId1) == 531) or abs((gen_bLeadId1) == 541) or abs((gen_bLeadId1) == 551) or abs((gen_bLeadId1) == 553): weight1 = p1_b

    elif abs((gen_bLeadId2) == 531) or abs((gen_bLeadId2) == 541) or abs((gen_bLeadId2) == 551) or abs((gen_bLeadId2) == 553): weight2 = p2_b
    #Lambda and b baryons
    elif (abs(gen_bLeadId1) == 5122) or (abs(gen_bLeadId1) == 5232) or (abs(gen_bLeadId1) == 5132) or (abs(gen_bLeadId1) == 5332): weight1 = p1_l
    elif (abs(gen_bLeadId2) == 5122) or (abs(gen_bLeadId2) == 5232) or (abs(gen_bLeadId2) == 5132) or (abs(gen_bLeadId2) == 5332):
        weight2 = p2_l
    rjet_bhad.Fill(gen_bpt1,resp1,weight1)
    rjet_bhad.Fill(gen_bpt2,resp2,weight2)
    #B fragmentation
    rjet_bquark.Fill(gen_bpt1,resp1,p_bquark(x1))
    rjet_bquark.Fill(gen_bpt2,resp2,p_bquark(x2))
    flagsIds = [1,9,17,25,33,41,49,57,3,11,19,27,35,43,51,59]
    #Semileptonic
    weight3 = 1
    if gen_bFlags1 in flagsIds:
        if abs((gen_bLeadId1 == 521)):
            weight3 = b_p(had_bpt1)
        elif abs((gen_bLeadId1 == 511)):
            weight3 = bt_p(had_bpt1)
        elif abs((gen_bLeadId1) == 531):
            weight3 = bc_p(had_bpt1)
        elif (abs(gen_bLeadId1) == 5122):
            weight3 = bac_p(had_bpt1)
    weight4 = 1
    if gen_bFlags2 in flagsIds:
        if abs((gen_bLeadId2 == 521)):
            weight4 = b_p(had_bpt2)
        elif abs((gen_bLeadId2 == 511)):
            weight4 = bt_p(had_bpt2)
        elif abs((gen_bLeadId2) == 531):
            weight4 = bc_p(had_bpt2)
        elif (abs(gen_bLeadId2) == 5122):
            weight4 = bac_p(had_bpt2)
    rjet_bsl.Fill(gen_bpt1,resp1,weight3)
    rjet_bsl.Fill(gen_bpt2,resp2,weight4)
#Ratios
h1 = rjet_bquark.ProjectionX()
h2 = rjet.ProjectionX()
h3 = rjet_bhad.ProjectionX()
h4 = rjet_bsl.ProjectionX()
h1.Divide(h2)
h3.Divide(h2)

h4.Divide(h2)
#Markers
h1.SetMarkerStyle(20)
h1.SetMarkerSize(1)
h3.SetMarkerStyle(21)
h3.SetMarkerSize(1)
h4.SetMarkerStyle(22)
h4.SetMarkerSize(1)
#Color
root.gStyle.SetPalette(root.kRainBow)
#Legend
legend = TLegend(0.7,0.75,0.8,0.85)
legend.SetBorderSize(0)
legend.SetTextSize(0.025)
legend.AddEntry(h1,  "B fragmentation", "p")
legend.AddEntry(h3,  "B production fraction", "p")
legend.AddEntry(h4,  "B to SL BR", "p")
#Title and axis setup
setup = TH1D("setup","", rjet.GetXaxis().GetNbins(), 28, 1240)
setup.SetStats(0) #Suppress stat box
setup.SetAxisRange(0.995,1.035,"Y") #Vertical axis limits
setup.SetAxisRange(0,400,"X") #Horizontal axis limits
setup.GetXaxis().SetMoreLogLabels()
setup.GetXaxis().SetNoExponent()
setup.GetXaxis().SetTitleOffset(1.3)
setup.GetYaxis().SetTitle("Response ratio")
setup.GetYaxis().SetTitleSize(0.03)
setup.GetXaxis().SetTitle("p_{T}^{gen} (GeV)")
setup.GetXaxis().SetTitleSize(0.03)
setup.SetTitle(" B jet response")
setup.Draw()
#Plot
h1.Draw("PLC PMC same")
h3.Draw("PLC PMC same")
h4.Draw("PLC PMC same")
legend.Draw("PLC PMC same")
c1.Print("bjetsresponse.pdf")
#Save plot
myfile = TFile( 'bjetsresponse.root', 'RECREATE' )
rjet_bquark.Write()
myfile.Close()
