from ROOT import *

gStyle.SetOptStat(0)

gROOT.SetBatch(True)

f = TFile("hitplots.root")

topdir = f.Get('plotter')
evtdir = topdir.Get('evt_2302')

c = TCanvas()
gPad.Print('traces.pdf[')

for wire in evtdir.GetListOfKeys():
    wiredir = evtdir.Get(wire.GetName())

    print wire
    for rangekey in wiredir.GetListOfKeys():
        print '  '+rangekey.GetName()
        rng = wiredir.Get(rangekey.GetName())

        gausdir = rng.Get('gaushit')
        mydir = rng.Get('myhitfinder')

        # Only show cases where the reconstruction is different
        if gausdir and mydir and len(gausdir.GetListOfKeys()) == len(mydir.GetListOfKeys()): continue

        # Skip the very common case where myhitfinder finds a small extra hit
        # in the raised pedestal to the left.
        if gausdir and mydir and len(gausdir.GetListOfKeys()) == 1 and len(mydir.GetListOfKeys()) == 2: continue

        hwire = rng.Get('wire')
        hwire.SetLineColor(kBlack)
        hwire.SetLineWidth(2)
        hwire.SetTitle(wire.GetName())
        hwire.Draw('hist')

        hdig = rng.Get('dig')
        hdig.Scale(.5) # arbitrary, try to make more similar to wires
        hdig.SetLineColor(kGray)
        hdig.SetLineWidth(2)
        hdig.SetTitle(wire.GetName())
        hdig.Draw('hist same')

        hwire.Draw('hist same')


        if gausdir:
            for gkey in gausdir.GetListOfKeys():
                g = gausdir.Get(gkey.GetName())
                g.Draw('c same')
                g.SetLineColor(kRed)

        if mydir:
            for gkey in mydir.GetListOfKeys():
                g = mydir.Get(gkey.GetName())
                g.Draw('c same')
                g.SetLineColor(kBlue)

# Fun in non-batch mode
#        gPad.Update()

        gPad.Print('traces.pdf')

gPad.Print('traces.pdf]')
