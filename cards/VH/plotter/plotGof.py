#!/usr/bin/env python

import ROOT as rt
rt.gSystem.Load('/uscms_data/d3/mwadud/combine_mod/CMSSW_14_1_0_pre4/src/CombineHarvester/CombineTools/test/datacards_VH_paper/feb_19_undecorrelate_jmr_jms/plotter/extra_tools_cc.so')
rt.gROOT.SetBatch()
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)
rt.gStyle.SetLineScalePS(1)
rt.gStyle.SetCanvasPreferGL(1)


## https://www.physics.ucla.edu/~cousins/stats/cousins_saturated.pdf

print("Starting GoF plotter @ " + rt.getCurrentTime() + '\n') 

options 		= rt.parseOptions('plotGoFOptions.txt', '===', 0)

toysGoFtree		= rt.openTree(options.get('toysGoFfile'), 'limit')
toysHist 		= rt.TH1D('toysHist', ';' + options.get("xtitle") + ';Number of toys', options.getIntList('toyHistBinning')[0], options.getDoubleList('toyHistBinning')[1], options.getDoubleList('toyHistBinning')[2])
drawString  	= 'min( max(' + rt.to_string_with_precision(toysHist.GetBinCenter(1),6) + ', limit), ' + rt.to_string_with_precision(toysHist.GetBinCenter(toysHist.GetNbinsX()),6) + ')'
nToys 			= toysGoFtree.Draw(drawString + ' >> toysHist', options.get('cuts'), 'goff')
print('N toys = %d' % nToys)
toysHist.SetLineColor(options.getTColFromHex('toysLineColor'))
toysHist.SetLineWidth(options.getInt('toysLineWidth'))
toysHist.GetXaxis().CenterTitle()
toysHist.GetYaxis().CenterTitle()
toysHist.GetXaxis().SetTitleSize(options.getDouble('canvasAxisTitleSize'))
toysHist.GetYaxis().SetTitleSize(options.getDouble('canvasAxisTitleSize'))
toysHist.GetXaxis().SetLabelSize(options.getDouble('canvasAxisLabelSize'))
toysHist.GetYaxis().SetLabelSize(options.getDouble('canvasAxisLabelSize'))
toysHist.GetYaxis().SetLabelSize(options.getDouble('canvasAxisLabelSize'))
toysHist.GetXaxis().SetTitleOffset(options.getDouble('canvasXtitleOffset'))
toysHist.GetYaxis().SetTitleOffset(options.getDouble('canvasYtitleOffset'))

obsGoFtree 		= rt.openTree(options.get('obsGoFfile'), 'limit')
nObs 			=  obsGoFtree.Draw('limit', options.get('cuts'), 'goff')
print('N observed = %d' % nObs)
assert(nObs == 1)
obsChi2 		= obsGoFtree.GetV1()[0]
print('Observed chi2 = %f' % obsChi2)
obsGr  			= rt.TArrow(obsChi2, 0, obsChi2, options.getFloat('obsMaxHratio')*toysHist.GetMaximum(), options.getFloat('obsArrowSize'), '<|')
obsGr.SetFillColor(options.getTColFromHex('obsLineColor'))
obsGr.SetLineColor(options.getTColFromHex('obsLineColor'))
obsGr.SetLineWidth(options.getInt('obsLineWidth'))

pValHist 		= rt.copyHistSubRange(toysHist, obsChi2, options.getDoubleList('toyHistBinning')[2], 'pValHist')
pValHist.SetFillColor(options.getTColFromHex('obsLineColor'))
pValHist.SetFillStyle(options.getInt('pValHistFillStyle'))
pValHist.SetLineWidth(0)

canvas 			= rt.TCanvas('canvas', '', options.getInt('canvasX'), options.getInt('canvasY'))
canvas.SetMargin(options.getDouble('canvasMarginL'), options.getFloat('canvasMarginR'), options.getFloat('canvasMarginB'), options.getFloat('canvasMarginT'))
canvas.SetFillStyle(4000)
canvas.SetFillColor(0)
canvas.SetFrameFillStyle(4000)
canvas.SetTicks(1,1)
canvas.Draw()

toysHist.Draw('HIST')
pValHist.Draw('HIST SAME')
toysHist.Draw('HIST SAME')
obsGr.Draw()

toysHist.GetXaxis().SetRangeUser(toysHist.GetXaxis().GetBinLowEdge(1), toysHist.GetXaxis().GetBinLowEdge(toysHist.GetNbinsX()-2))

rt.cmsLegend(canvas, options.getList('cmsLegend', ';')[0], options.getList('cmsLegend', ';')[1], options.getFloat('cmsTitleScale'));

legend 			= rt.TLegend(options.getDouble('legx1'), options.getDouble('legy1'), options.getDouble('legx2'), options.getDouble('legy2'))
legend.SetHeader(options.get("legendHeader"), 'C')
legHeader 		= legend.GetListOfPrimitives().First();
legHeader.SetTextSize(options.getDouble('legHeaderTextSize'))
legend.SetTextSize(options.getDouble('legTextSize'))
legend.SetNColumns(1)
legend.SetFillStyle(4000)
legend.SetLineWidth(0)


toysLeg 		= str(nToys) + ' toys : #mu=' + rt.to_string_with_precision(toysHist.GetMean(), 1) + ', #sigma=' + rt.to_string_with_precision(toysHist.GetStdDev(), 1)
legend.AddEntry(toysHist, toysLeg, 'L')

pVal 			= toysHist.Integral(toysHist.FindBin(obsChi2), options.getIntList('toyHistBinning')[0])/toysHist.Integral()
obsLeg 			= 'Obs. = ' + rt.to_string_with_precision(obsChi2, 2)
legend.AddEntry(obsGr, obsLeg, 'L')

pValHistLeg  	= 'p-value = ' + rt.to_string_with_precision(pVal, 2)
legend.AddEntry(pValHist, pValHistLeg, 'F')

canvas.RedrawAxis()
canvas.Update()
canvas.Modified()

legend.Draw()

saveName 		= options.get('writeDir') + '/' + options.get('tag')
canvas.SaveAs(saveName + '.png')
canvas.SaveAs(saveName + '.pdf')

rt.clearHeap()