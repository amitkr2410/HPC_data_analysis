import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt

# define format for the plots
import matplotlib as mpl

mpl.rc('font', family='serif')
mpl.rcParams['text.usetex'] = True
params = {'text.latex.preamble' : r'\usepackage{amsmath} \usepackage{lmodern} \boldmath'}
mpl.rcParams.update(params)

mpl.rcParams['figure.figsize'] = [6., 4.5]
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 7
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.labelsize'] = 17
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.labelsize'] = 17
mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rcParams['ytick.minor.visible'] = False
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams['legend.fontsize'] = 17
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['font.size'] = 17
# mpl.rcParams['font.weight'] = 'black'
# mpl.rcParams['axes.labelweight'] = 'black'
mpl.rcParams['savefig.format'] = "pdf"

########################################################################################################################################################
# Set Color, Line Style, and Markers
##--

color_theory = ['red','blue','green']
line_theory = ['solid','dashed','dotted']

color_data = ['black']#['indianred','cadetblue','darkolivegreen','peru']
marker_data = ['o','s','^','X'] ## o: circle, s: square, ^:triangle, X:cross
########################################################################################################################################################
########################################################################################################################################################
# Set File Path and labels
##--
file_denom = '../PartonicJet/ColorlessJetCrossSection_R_0p4_AntikT_ppPbPb_Jetscape_STAR_200GeV_EtaJetMin_0p0_EtaJetMax_1p0_FS_allHadron_Cent_0-10_Matter_pp_BGS0_pTTrackMin_0p0GeV.txt'

files_theory = ['../PartonicJet/ColorlessJetCrossSection_R_0p4_AntikT_ppPbPb_Jetscape_STAR_200GeV_EtaJetMin_0p0_EtaJetMax_1p0_FS_allParton_Cent_0-10_Matter_pp_BGS0_pTTrackMin_0p0GeV.txt']
                

label_theory = [r'$\mathrm{JETSCAPE~3.0}$']


files_data = ['']
#ATLAS PLB 790, 108 (2019)
label_data = ['$\mathrm{ATLAS~[PLB~790,~108~(2019)]}$']
OutputFilename='../Plots/RatioPlot_DJCS_PartonicJet_200GeV_STAR_R0p4_eta1p0.pdf'
########################################################################################################################################################          

########################################################################################################################################################
# Set Multiplication Factors
##--
factor_theory = [1,1,1]
factor_data = factor_theory
########################################################################################################################################################

########################################################################################################################################################
# Fucntions to load ascii format files
##--

def GetTheory(filename):
    data = np.loadtxt(filename, comments='#')
    x = data[:,0]
    y = data[:,1]
    xerr = data[:,2]    
    yerr = data[:,3]
    xstep = np.append(x-xerr,x[-1]+xerr[-1])
    ystep = np.append(y,y[-1])    

    return x, y, xerr, yerr, xstep, ystep


def GetExp(filename):
    data = np.loadtxt(filename, comments='#')
    x = data[:,0]
    y = data[:,1]
    xerrl = data[:,2]    
    xerrh = data[:,3]        
    yerr = data[:,4]
    ysysl = np.append(data[:,5],data[-1,5])
    ysysh = np.append(data[:,6],data[-1,6])
    xstep = np.append(x-xerrl,x[-1]+xerrh[-1])
    ystep = np.append(y,y[-1])    

    return x, y, 0.5*(xerrl+xerrh), yerr, xstep, ystep, ysysl, ysysh


##---
def GetRatioTheoryToTheory(num_file, den_file):
    x_num, y_num, xerr_num, yerr_num, xstep_num, ystep_num = GetTheory(num_file)
    x_den, y_den, xerr_den, yerr_den, xstep_den, ystep_den = GetTheory(den_file)  
  
    return x_num, y_num/y_den, xerr_num, yerr_num/y_den, xstep_num, ystep_num/ystep_den


##---
def GetRatioExpToTheory(num_file, den_file):
    x_num, y_num, xerr_num, yerr_num, xstep_num, ystep_num, ysysl_num, ysysh_num = GetExp(num_file)
    x_den, y_den, xerr_den, yerr_den, xstep_den, ystep_den = GetTheory(den_file)  
  
    return x_num, y_num/y_den, xerr_num, yerr_num/y_den, xstep_num, ystep_num/ystep_den, ysysl_num/ystep_den, ysysh_num/ystep_den

########################################################################################################################################################

########################################################################################################################################################
# Fucntions to add JETSCAPE Logo
##--
from reportlab.pdfgen import canvas
from PyPDF2 import PdfFileWriter, PdfFileReader
import sys

def add_jetscape_logo(filename, x0, y0, dx, tag=''):

    input_filename = filename+'.pdf'

    dy = dx*(67.0/100.0)
    c = canvas.Canvas('temp.pdf')
    c.drawImage('JetscapeLogo.jpg', x0,y0,dx,dy)
    c.save()

    output = PdfFileWriter()
    input1 = PdfFileReader(open(input_filename, "rb"))
    watermark = PdfFileReader(open("temp.pdf", "rb"))

    input_page = input1.getPage(0)
    input_page.mergePage(watermark.getPage(0))
    output.addPage(input_page)

    output_filename = filename+tag+'.pdf'
    # finally, write "output" to document-output.pdf
    outputStream = open(output_filename, "wb")
    output.write(outputStream)
    outputStream.close()

########################################################################################################################################################
########################################################################################################################################################
# Main Plot Code 1. Multiple Raw Graph
##--
import os
fig, axs = plt.subplots(2, 1, figsize=(6, 6.5), sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0,'height_ratios': (5,2.5)})


## Log Scale for xy axes
axs[0].set_yscale('log')

## Domain and Range
x_min = 50
x_max = 300

y0_min = 1.0e-13
y0_max = 5.0e-3

y1_min = 0.5
y1_max = 1.25
axs[0].set_xlim(x_min,x_max)
axs[0].set_ylim(y0_min,y0_max)
axs[1].set_ylim(y1_min,y1_max)

axs[1].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(25))
axs[1].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
axs[1].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))

## Label for xy axes
axs[1].set_xlabel(r'$p_{\mathrm{T}}^{\mathrm{jet}}~\mathrm{(GeV)}$', fontsize=18)
axs[0].set_ylabel(r'$\frac{d^2 \sigma}{dy dp_{\mathrm{T}}^{\mathrm{jet}}}~\mathrm{(mb/GeV)}$', fontsize=18)
axs[1].set_ylabel(r'$\mathrm{Ratio~to}$'+'\n'+'$\mathrm{PYTHIA~8.230}$', fontsize=12)


############################################################################################################################################
## Plot Theory Lines
for i, file in enumerate(files_theory):
  x, y, xerr, yerr, xstep, ystep = GetTheory(file)
  # Steps
  axs[0].step(xstep, factor_theory[i]*ystep, where='post', color=color_theory[i],linestyle=line_theory[i])
  # Error bars for stat errors  
  axs[0].errorbar(x, factor_theory[i]*y, factor_theory[i]*yerr, marker="", linestyle="none", color=color_theory[i])

## Plot Exp Data Points
#for i, file in enumerate(files_data):
#  x, y, xerr, yerr, xstep, ystep, ysysl, ysysh = GetExp(file)
  # Markers with error bars for stat errors
  #axs[0].errorbar(x, factor_theory[i]*y, yerr=factor_theory[i]*yerr, linestyle="none", color = color_data[i], marker=marker_data[i])
  # Shades for sys errors
  #axs[0].fill_between(xstep, factor_theory[i]*(ystep-ysysl), factor_theory[i]*(ystep+ysysh), step='post', alpha=0.2, color = color_data[i])
############################################################################################################################################
## Plot Theory Lines
for i, file in enumerate(files_theory):
  x, y, xerr, yerr, xstep, ystep = GetRatioTheoryToTheory(file,file_denom)
  # Steps
  axs[1].step(xstep, factor_theory[i]*ystep, where='post', color=color_theory[i],linestyle=line_theory[i])
  # Error bars for stat errors  
  axs[1].errorbar(x, factor_theory[i]*y, factor_theory[i]*yerr, marker="", linestyle="none", color=color_theory[i])

# Plot Exp Data Points
#for i, file in enumerate(files_data):
#  x, y, xerr, yerr, xstep, ystep, ysysl, ysysh = GetRatioExpToTheory(file,file_denom)
  # Markers with error bars for stat errors
 # axs[1].errorbar(x, factor_theory[i]*y, yerr=factor_theory[i]*yerr, linestyle="none", color = color_data[i], marker=marker_data[i])
  # Shades for sys errors
#  axs[1].fill_between(xstep, factor_theory[i]*(ystep-ysysl), factor_theory[i]*(ystep+ysysh), step='post', alpha=0.2, color = color_data[i])
############################################################################################################################################
# Legends
handles = []
labels = []
#for i, file in enumerate(files_data):
#  dp = axs[0].errorbar(0, 0, 0, linestyle="none", color = color_data[i], marker=marker_data[i]) 
#  dsys = axs[0].fill(np.NaN, np.NaN, alpha=0.2, color = color_data[i])
#  handles.append((dp[0],dsys[0]))
#  labels.append(label_data[i])

for i, file in enumerate(files_theory):
  tl = axs[0].errorbar(0, 0, color=color_theory[i],linestyle=line_theory[i])
  handles.append(tl[0]) 
  labels.append(label_theory[i])   
axs[0].legend(handles,labels,ncol=1,loc='lower left',edgecolor='none', frameon=True, facecolor='none', handletextpad=0.4, handleheight=1.8, labelspacing=0.05, bbox_to_anchor=(0, 0), borderaxespad=0.5, handlelength=1.6, fontsize=11)

## Text
axs[0].text(1.2*x_min,0.3*y0_max, '$\mathrm{pp~}\sqrt{s}=5.02\mathrm{~TeV}$\n'+r'$\mathrm{anti\text{-}}k_\mathrm{T}$, $R\!=\!0.4$'+'\n'+'$|y_{jet}|<0.3$', horizontalalignment='left', verticalalignment='top')

## Generate PDF File
plt.tight_layout()
plt.savefig('plot_template_ratio2')
## Add Logos
add_jetscape_logo('plot_template_ratio2',308,385,80,'_logo')
os.remove('plot_template_ratio2.pdf')
########################################################################################################################################################

########################################################################################################################################################
########################################################################################################################################################
# Main Plot Code 2. Single Raw Graph
##--
fig, axs = plt.subplots(1, 1, figsize=(6,3.5), sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})

## Log Scale for xy axes
# axs.set_xscale('log')
# axs.set_yscale('log')

## Domain and Range
x_min = 5
x_max = 60
y_min = 0.7
y_max = 3.5
axs.set_xlim(x_min,x_max)
axs.set_ylim(y_min,y_max)

axs.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
axs.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
axs.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.4))
axs.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.2))


## Label for xy axes
axs.set_xlabel(r'$p_{\mathrm{T}}^{\mathrm{jet}}~\mathrm{(GeV)}$', fontsize=19)
axs.set_ylabel(r'$\frac{\mathrm{Partonic}}{\mathrm{Hadronic}}~\left[\frac{d^2 \sigma}{d\eta dp_{\mathrm{T}}^{\mathrm{jet}}}\right]$', fontsize=19)


## Plot Theory Lines
for i, file in enumerate(files_theory):
  x, y, xerr, yerr, xstep, ystep = GetRatioTheoryToTheory(file,file_denom)
  # Steps
  axs.step(xstep, factor_theory[i]*ystep, where='post', color=color_theory[i],linestyle=line_theory[i])
  # Error bars for stat errors  
  axs.errorbar(x, factor_theory[i]*y, factor_theory[i]*yerr, marker="", linestyle="none", color=color_theory[i])

# Plot Exp Data Points
#for i, file in enumerate(files_data):
#  x, y, xerr, yerr, xstep, ystep, ysysl, ysysh = GetRatioExpToTheory(file,file_denom)
#  # Markers with error bars for stat errors
#  axs.errorbar(x, factor_theory[i]*y, yerr=factor_theory[i]*yerr, linestyle="none", color = color_data[i], marker=marker_data[i])
  # Shades for sys errors
#  axs.fill_between(xstep, factor_theory[i]*(ystep-ysysl), factor_theory[i]*(ystep+ysysh), step='post', alpha=0.2, color = color_data[i])


## Legentds
handles = []
labels = []
#for i, file in enumerate(files_data):
#  dp = axs.errorbar(0, 0, 0, linestyle="none", color = color_data[i], marker=marker_data[i]) 
#  dsys = axs.fill(np.NaN, np.NaN, alpha=0.2, color = color_data[i])
#  handles.append((dp[0],dsys[0]))
#  labels.append(label_data[i])

for i, file in enumerate(files_theory):
  tl = axs.errorbar(0, 0, color=color_theory[i],linestyle=line_theory[i])
  handles.append(tl[0]) 
  labels.append(label_theory[i])   
axs.legend(handles,labels,ncol=1,loc='lower left',edgecolor='none', frameon=True, facecolor='none', handletextpad=0.4, handleheight=1.8, labelspacing=0.05, bbox_to_anchor=(0.05, 0), borderaxespad=0.5, handlelength=1.6, fontsize=14)

## Text
axs.text(1.4*x_min+20,0.95*y_max, '$\mathrm{pp},~\sqrt{s}=200\mathrm{~GeV}$\n'+r'$\mathrm{anti\text{-}}k_\mathrm{T}$, $R\!=\!0.4$'+'\n'+'$|\eta_\mathrm{jet}|<1.0$', horizontalalignment='left', verticalalignment='top')

## Generate PDF File
plt.tight_layout()
plt.savefig('temp2')
## Add Logos
#add_jetscape_logo('temp2',290,160,80,'_logo')
#os.remove('temp.pdf')
#RenameFileName='mv temp2_logo.pdf' + '  ' + OutputFilename
RenameFileName='mv temp2.pdf' + '  ' + OutputFilename
os.system(RenameFileName)
#os.remove('temp2.pdf')
CommandOpenFile='open ' + OutputFilename
os.system(CommandOpenFile)
########################################################################################################################################################
