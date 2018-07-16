import ROOT as r
import re
import numpy as n
import os

Params = ['emax', 'e2max', 'edmax', 'eocore', 'w3st', 'w21st']
Particles = ["pi0","photon"]
############
Colours = {}
Colours["pi0"] = r.kBlue
Colours["photon"] = r.kOrange
############
Var = {}
Var['e2max']="e2max"
Var['edmax']="edmax"
Var['eocore']="eocore"
Var['wxst']="w3st"
Var['emax']="emax"
#Var['w21st']="w21st"

############
Paths = {}
Paths["pi0"] = "/afs/cern.ch/user/r/rastein/public/test/pi0/pi0_20k.root"
Paths["photon"] = "/afs/cern.ch/user/r/rastein/public/test/photon/photon_20k.root"
'''Path1['emax']='/afs/cern.ch/user/r/rastein/public/pionrejection/pion/10GeV/pi0_20k.root'
Path1['e2max']='/afs/cern.ch/user/r/rastein/public/pionrejection/pion/10GeV/pi0_20k.root'
Path1['edmax']='/afs/cern.ch/user/r/rastein/public/pionrejection/pion/10GeV/pi0_20k.root'
Path1['eocore']='/afs/cern.ch/user/r/rastein/public/pionrejection/pion/10GeV/pi0_20k.root'
Path1['wxst']='/afs/cern.ch/user/r/rastein/public/pionrejection/pion/10GeV/pi0_20k.root'
#Path1['w21st']='/afs/cern.ch/user/r/rastein/public/pionrejection/pion/10GeV/pi0_20k.root'
'''
Path2 = '/afs/cern.ch/user/r/rastein/public/pionrejection/photon/10GeV/photon_20k.root'
'''Path2['emax']='/afs/cern.ch/user/r/rastein/public/pionrejection/photon/10GeV/photon_20k.root'
Path2['e2max']='/afs/cern.ch/user/r/rastein/public/pionrejection/photon/10GeV/photon_20k.root'
Path2['edmax']='/afs/cern.ch/user/r/rastein/public/pionrejection/photon/10GeV/photon_20k.root'
Path2['eocore']='/afs/cern.ch/user/r/rastein/public/pionrejection/photon/10GeV/photon_20k.root'
Path2['wxst']='/afs/cern.ch/user/r/rastein/public/pionrejection/photon/10GeV/photon_20k.root'
#Path2['w21st']='/afs/cern.ch/user/r/rastein/public/pionrejection/photon/10GeV/photon_20k.root'
'''

############
Title = {}
Title['emax']='Maximum energy, emax'
Title['e2max']='Second maximum energy, e2max'
Title['edmax']='Energy in cell with minimal deposit between Emax & E2max, edmax'
Title['eocore']='Energy deposited outside the shower core, eocore'
Title['w3st']='Shower width over 3 cells, w3st'
Title['w21st']='Shower width over 21 cells, w21st'



file = r.TFile.Open(Paths["pi0"])
file2 = r.TFile.Open(Paths["photon"])
tree = file.Get("events")
tree2 = file2.Get("events")

canvas = r.TCanvas("Parameters", "Parameters", 900, 900)
canvas.Divide(2,3)
pad1 = r.TPad('pad1','This is pad1',0.02,0.68,0.48,0.98,0)
pad2 = r.TPad('pad2','This is pad2',0.52,0.68,0.98,0.98,0)
pad3 = r.TPad('pad3','This is pad3',0.02,0.36,0.48,0.66,0)
pad4 = r.TPad('pad4','This is pad4',0.52,0.36,0.98,0.66,0)
pad5 = r.TPad('pad5','This is pad5',0.02,0.02,0.48,0.34,0)
pad6 = r.TPad('pad6','This is pad6',0.52,0.02,0.98,0.34,0)
pad1.Draw()
pad2.Draw()
pad3.Draw()
pad4.Draw()
pad5.Draw()
pad6.Draw()

############
col1 = Colours["pi0"]
col2 = Colours["photon"]

for p in Params:
    title = Title[p]
    Entries_pi0 = []
    Entries_photon = []
    particle1 = Particles[0]
    particle2 = Particles[1]
    bins = 48

    norm = int(tree.GetEntries())
    unit = 'Gev'

    if p == "emax":

        pad = pad1
        pad.cd()
        pad.SetGrid()

        for i in tree:
            Entries_pi0.append(i.emax)
        for i in tree2:
            Entries_photon.append(i.emax)

        vals_pi0 = n.array(Entries_pi0)
        vals_photon = n.array(Entries_photon)
        meanval = (n.mean(vals_pi0)+n.mean(vals_photon))/2
        #bins = (len(n.unique(vals_pi0)) + len(n.unique(vals_photon)))/2
        hist_pi01 = r.TH1F(p+particle1,p+particle1,bins,0.0,5.0*meanval)
        hist_photon1 = r.TH1F(p+particle2,p+particle2,bins,0.0,5.0*meanval)
        for i in tree:
            hist_pi01.Fill(i.emax)
        for i in tree2:
            hist_photon1.Fill(i.emax)

        hist_pi01.SetLineWidth(2)
        hist_photon1.SetLineWidth(2)

        hist_pi01.SetLineColor(col1)
        hist_photon1.SetLineColor(col2)

        #hist_pi01.Scale(norm)
        #hist_photon1.Scale(norm)

        hist_pi01.GetYaxis().SetTitleOffset(1.40)
        hist_pi01.GetXaxis().SetTitleOffset(1.40)
        hist_pi01.SetTitle("")
        hist_pi01.GetXaxis().SetTitle(title)
        hist_pi01.GetYaxis().SetTitle("events")

        leg1 = r.TLegend(0.60,0.70,0.86,0.88)
        leg1.SetFillColor(0)
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(10)
        leg1.SetTextSize(0.065)
        leg1.SetTextFont(42)
        leg1.AddEntry(hist_pi01,"pi0")
        leg1.AddEntry(hist_photon1,"photon")


        hist_pi01.Draw()
        hist_pi01.SetDirectory(0)
        hist_photon1.Draw("same")
        hist_photon1.SetDirectory(0)
        leg1.Draw("same")

    elif p == "e2max":

        pad = pad2
        pad.cd()
        pad.SetGrid()

        for i in tree:
            Entries_pi0.append(i.e2max)
        for i in tree2:
            Entries_photon.append(i.e2max)

        vals_pi0 = n.array(Entries_pi0)
        vals_photon = n.array(Entries_photon)
        meanval = (n.mean(vals_pi0)+n.mean(vals_photon))/2
        #bins = (len(n.unique(vals_pi0)) + len(n.unique(vals_photon)))/2
        hist_pi02 = r.TH1F(p+particle1,p+particle2,bins,0.0,5.0*meanval)
        hist_photon2 = r.TH1F(p+particle2,p+particle2,bins,0.0,5.0*meanval)
        for i in tree:
            hist_pi02.Fill(i.e2max)
        for i in tree2:
            hist_photon2.Fill(i.e2max)

        hist_pi02.SetLineWidth(2)
        hist_photon2.SetLineWidth(2)

        hist_pi02.SetLineColor(col1)
        hist_photon2.SetLineColor(col2)

        #hist_pi02.Scale(norm)
        #hist_photon2.Scale(norm)

        hist_pi02.GetYaxis().SetTitleOffset(1.40)
        hist_pi02.GetXaxis().SetTitleOffset(1.40)
        hist_pi02.SetTitle("")
        hist_pi02.GetXaxis().SetTitle(title)
        hist_pi02.GetYaxis().SetTitle("events")

        leg2 = r.TLegend(0.60,0.70,0.86,0.88)
        leg2.SetFillColor(0)
        leg2.SetFillStyle(0)
        leg2.SetLineColor(0)
        leg2.SetShadowColor(10)
        leg2.SetTextSize(0.065)
        leg2.SetTextFont(42)
        leg2.AddEntry(hist_pi02,"pi0")
        leg2.AddEntry(hist_photon2,"photon")

        hist_pi02.Draw()
        hist_pi02.SetDirectory(0)
        hist_photon2.Draw("same")
        hist_photon2.SetDirectory(0)
        leg2.Draw("same")

    elif p == "edmax":

        pad = pad3
        pad.cd()
        pad.SetGrid()

        for i in tree:
            Entries_pi0.append(i.edmax)
        for i in tree2:
            Entries_photon.append(i.edmax)

        vals_pi0 = n.array(Entries_pi0)
        vals_photon = n.array(Entries_photon)
        meanval = (n.mean(vals_pi0)+n.mean(vals_photon))/2
        #bins = (len(n.unique(vals_pi0)) + len(n.unique(vals_photon)))/2
        hist_pi03 = r.TH1F(p+particle1,p+particle1,bins+16,0.0,10.0*meanval)
        hist_photon3 = r.TH1F(p+particle2,p+particle2,bins+16,0.0,10.0*meanval)
        for i in tree:
            hist_pi03.Fill(i.edmax)
        for i in tree2:
            hist_photon3.Fill(i.edmax)

        hist_pi03.SetLineWidth(2)
        hist_photon3.SetLineWidth(2)
        pad.SetLogy()
        hist_pi03.SetLineColor(col1)
        hist_photon3.SetLineColor(col2)

        #hist_pi03.Scale(norm)
        #hist_photon3.Scale(norm)

        hist_photon3.GetYaxis().SetTitleOffset(1.40)
        hist_photon3.GetXaxis().SetTitleOffset(1.40)
        hist_photon3.SetTitle("")
        hist_photon3.GetXaxis().SetTitle(title)
        hist_photon3.GetYaxis().SetTitle("events")

        leg3 = r.TLegend(0.60,0.70,0.86,0.88)
        leg3.SetFillColor(0)
        leg3.SetFillStyle(0)
        leg3.SetLineColor(0)
        leg3.SetShadowColor(10)
        leg3.SetTextSize(0.065)
        leg3.SetTextFont(42)
        leg3.AddEntry(hist_pi03,"pi0")
        leg3.AddEntry(hist_photon3,"photon")

        hist_photon3.Draw()
        hist_photon3.SetDirectory(0)
        hist_pi03.Draw("same")
        hist_pi03.SetDirectory(0)
        leg3.Draw("same")


    elif p == "eocore":

        pad = pad4
        pad.cd()
        pad.SetGrid()

        for i in tree:
            Entries_pi0.append(i.eocore)
        for i in tree2:
            Entries_photon.append(i.eocore)

        vals_pi0 = n.array(Entries_pi0)
        vals_photon = n.array(Entries_photon)
        meanval = (n.mean(vals_pi0)+n.mean(vals_photon))/2
        #bins = (len(n.unique(vals_pi0)) + len(n.unique(vals_photon)))/2
        hist_pi04 = r.TH1F(p+particle1,p+particle1,bins,0.0,2.0*meanval)
        hist_photon4 = r.TH1F(p+particle2,p+particle2,bins,0.0,2.0*meanval)
        for i in tree:
            hist_pi04.Fill(i.eocore)
        for i in tree2:
            hist_photon4.Fill(i.eocore)

        hist_pi04.SetLineWidth(2)
        hist_photon4.SetLineWidth(2)

        #hist_pi04.Scale(norm)
        #hist_photon4.Scale(norm)

        hist_pi04.SetLineColor(col1)
        hist_photon4.SetLineColor(col2)

        hist_photon4.GetYaxis().SetTitleOffset(1.40)
        hist_photon4.GetXaxis().SetTitleOffset(1.40)
        hist_photon4.SetTitle("")
        hist_photon4.GetXaxis().SetTitle(title)
        hist_photon4.GetYaxis().SetTitle("events")

        leg4 = r.TLegend(0.60,0.70,0.86,0.88)
        leg4.SetFillColor(0)
        leg4.SetFillStyle(0)
        leg4.SetLineColor(0)
        leg4.SetShadowColor(10)
        leg4.SetTextSize(0.065)
        leg4.SetTextFont(42)
        leg4.AddEntry(hist_pi04,"pi0")
        leg4.AddEntry(hist_photon4,"photon")

        hist_photon4.Draw()
        hist_photon4.SetDirectory(0)
        hist_pi04.Draw("same")
        hist_pi04.SetDirectory(0)
        leg4.Draw("same")


    elif p == "w3st":

        pad = pad5
        pad.cd()
        pad.SetGrid()

        for i in tree:
            Entries_pi0.append(i.w3st)
        for i in tree2:
            Entries_photon.append(i.w3st)

        vals_pi0 = n.array(Entries_pi0)
        vals_photon = n.array(Entries_photon)
        meanval = (n.mean(vals_pi0)+n.mean(vals_photon))/2
        #bins = (len(n.unique(vals_pi0)) + len(n.unique(vals_photon)))/2
        hist_pi05 = r.TH1F(p+particle1,p+particle1,bins,0.0,2.0*meanval)
        hist_photon5 = r.TH1F(p+particle2,p+particle2,bins,0.0,2.0*meanval)
        for i in tree:
            hist_pi05.Fill(i.w3st)
        for i in tree2:
            hist_photon5.Fill(i.w3st)

        hist_pi05.SetLineWidth(2)
        hist_photon5.SetLineWidth(2)

        #hist_pi05.Scale(norm)
        #hist_photon5.Scale(norm)

        hist_pi05.SetLineColor(col1)
        hist_photon5.SetLineColor(col2)

        hist_photon5.GetYaxis().SetTitleOffset(1.40)
        hist_photon5.GetXaxis().SetTitleOffset(1.40)
        hist_photon5.SetTitle("")
        hist_photon5.GetXaxis().SetTitle(title)
        hist_photon5.GetYaxis().SetTitle("events")

        leg5 = r.TLegend(0.60,0.70,0.86,0.88)
        leg5.SetFillColor(0)
        leg5.SetFillStyle(0)
        leg5.SetLineColor(0)
        leg5.SetShadowColor(10)
        leg5.SetTextSize(0.065)
        leg5.SetTextFont(42)
        leg5.AddEntry(hist_pi05,"pi0")
        leg5.AddEntry(hist_photon5,"photon")



        hist_photon5.Draw()
        hist_photon5.SetDirectory(0)
        hist_pi05.Draw("same")
        hist_pi05.SetDirectory(0)
        leg5.Draw("same")

    elif p == "w21st":

        pad = pad6
        pad.cd()
        pad.SetGrid()

        for i in tree:
            Entries_pi0.append(i.w21st)
        for i in tree2:
            Entries_photon.append(i.w21st)

        vals_pi0 = n.array(Entries_pi0)
        vals_photon = n.array(Entries_photon)
        meanval = (n.mean(vals_pi0)+n.mean(vals_photon))/2
        #bins = (len(n.unique(vals_pi0)) + len(n.unique(vals_photon)))/2
        hist_pi06 = r.TH1F(p+particle1,p+particle1,bins,0.0,2.0*meanval)
        hist_photon6 = r.TH1F(p+particle2,p+particle2,bins,0.0,2.0*meanval)
        for i in tree:
            hist_pi06.Fill(i.w21st)
        for i in tree2:
            hist_photon6.Fill(i.w21st)

        hist_pi06.SetLineWidth(2)
        hist_photon6.SetLineWidth(2)

        hist_pi06.SetLineColor(col1)
        hist_photon6.SetLineColor(col2)

        #hist_pi06.Scale(norm)
        #hist_photon6.Scale(norm)

        hist_photon6.GetYaxis().SetTitleOffset(1.40)
        hist_photon6.GetXaxis().SetTitleOffset(1.40)
        hist_photon6.SetTitle("")
        hist_photon6.GetXaxis().SetTitle(title)
        hist_photon6.GetYaxis().SetTitle("events")

        leg6 = r.TLegend(0.60,0.70,0.86,0.88)
        leg6.SetFillColor(0)
        leg6.SetFillStyle(0)
        leg6.SetLineColor(0)
        leg6.SetShadowColor(10)
        leg6.SetTextSize(0.065)
        leg6.SetTextFont(42)
        leg6.AddEntry(hist_pi06,"pi0")
        leg6.AddEntry(hist_photon6,"photon")

        hist_photon6.Draw()
        hist_photon6.SetDirectory(0)
        hist_pi06.Draw("same")
        hist_pi06.SetDirectory(0)
        leg6.Draw("same")

    else:
        continue

    maxpi0 = n.amax(vals_pi0)
    minpi0 = n.amin(vals_photon)
    leftText = "Unknown"

    r.gStyle.SetStatW(0.30)
    r.gStyle.SetStatH(0.20)
    r.gStyle.SetStatColor(42)



    pad.SetTitle(title)
    pad.SetTicks(1,1)
    #pad.SetLeftMargin(0.14)
    #pad.SetRightMargin(0.08)
    r.gStyle.SetOptStat(0)

    ################



    pad.RedrawAxis()
    pad.GetFrame().SetBorderSize( 12 )
    pad.Modified()


    canvas.Update()
canvas.Print()
'''Text = r.TLatex()

Text.SetNDC()
Text.SetTextAlign(31);
Text.SetTextSize(0.04)

text = '#it{' + leftText +'}'

Text.DrawLatex(0.90, 0.92, text)

rightText = re.split(",", rightText)
text = '#bf{#it{' + rightText[0] +'}}'

Text.SetTextAlign(12);
Text.SetNDC(r.kTRUE)
Text.SetTextSize(0.04)
Text.DrawLatex(0.18, 0.83, text)

text = '#bf{#it{' + rightText[1] +'}}'
Text.SetTextSize(0.065)
Text.DrawLatex(0.18, 0.78, text)'''
