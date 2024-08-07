import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt

# define format for the plots
import matplotlib as mpl
#serif
mpl.rc('font', family='times', weight='bold')
mpl.rcParams['text.usetex'] = True
params = {'text.latex.preamble' : r'\usepackage{amsmath} \usepackage{lmodern} \boldmath'}
mpl.rcParams.update(params)

mpl.rcParams["font.weight"]="bold"
mpl.rcParams["axes.labelweight"]="bold"
mpl.rcParams['figure.figsize'] = [6., 4.5]
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 7
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rcParams['ytick.minor.visible'] = False
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams['legend.fontsize'] = 15
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['font.size'] = 15
# mpl.rcParams['font.weight'] = 'black'
# mpl.rcParams['axes.labelweight'] = 'black'
mpl.rcParams['savefig.format'] = "pdf"


# In[3]:


########################################################################################################################################################
# Set Color, Line Style, and Markers
##--

color_theory = ['red','blue','green']
line_theory = ['solid','dashed','dotted']

color_data = ['black','magenta']
marker_data = ['o','s','^','X'] ## o: circle, s: square, ^:triangle, X:cross
########################################################################################################################################################


# In[4]:


########################################################################################################################################################
# Set File Path and labels
##--
file_pp_theory = '../ColorlessSingleHadronYield_pp_Jetscape_CMS_5020GeV_EtaMax_1p0_FS_ChargedHadron_BGS0.txt'

files_aa_theory = ['../Tc150MeV_170MeV_tau0/ColorlessSingleHadronYield_PbPb_Jetscape_CMS_5020GeV_EtaMax_1p0_FS_ChargedHadron_Cent_0-10_MatterLBT_withRecoil_Type5_Tc150MeV_BGS1.txt',     '../ColorlessSingleHadronYield_PbPb_Jetscape_CMS_5020GeV_EtaMax_1p0_FS_ChargedHadron_Cent_0-10_MatterLBT_withRecoil_Type5_alphas_0p3_Q0_2GeV_BGS1.txt', '../Tc150MeV_170MeV_tau0/ColorlessSingleHadronYield_PbPb_Jetscape_CMS_5020GeV_EtaMax_1p0_FS_ChargedHadron_Cent_0-10_MatterLBT_withRecoil_Type5_Tc170MeV_BGS1.txt']
                   
#$\hat{q}=\hat{q}^{fix}_{HTL}f(Q^{2}), \mathrm{JS(MATTER+LBT) }$
label_theory = ['$T_{c}=150~\mathrm{MeV}$', '$T_{c}=160~\mathrm{MeV}$', '$T_{c}=170~\mathrm{MeV}$']

files_data = ['../../ExperimentData_All/ExpData_CMS_5020_ChargedHadron_Raa_Eta1p0_Cent0-10.dat']

label_data = ['$\mathrm{CMS~[JHEP~1704,~039~(2017)]}$']
OutputFilename='../Plots/ChargedParticle_RAA_Type5_different_Tc_5020GeV_CMS_0-10_eta1p0.pdf'
x_min =8.0
x_max =450
y_min = 0
y_max = 1.7
########################################################################################################################################################          


# In[5]:


########################################################################################################################################################
# Set Multiplication Factors
##--
factor_theory = [1,1,1]
factor_data = factor_theory
########################################################################################################################################################


# In[6]:


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
def GetRatioTheoryToTheory(num_file, den_file, error=False ):
    x_num, y_num, xerr_num, yerr_num, xstep_num, ystep_num = GetTheory(num_file)
    x_den, y_den, xerr_den, yerr_den, xstep_den, ystep_den = GetTheory(den_file)  
    yerr = yerr_num/y_den
    if error:
        yerr = RatioError(y_num,yerr_num,y_den,yerr_den)
    return x_num, y_num/y_den, xerr_num, yerr, xstep_num, ystep_num/ystep_den


##---
def GetRatioExpToTheory(num_file, den_file):
    x_num, y_num, xerr_num, yerr_num, xstep_num, ystep_num, ysysl_num, ysysh_num = GetExp(num_file)
    x_den, y_den, xerr_den, yerr_den, xstep_den, ystep_den = GetTheory(den_file)  
  
    return x_num, y_num/y_den, xerr_num, yerr_num/y_den, xstep_num, ystep_num/ystep_den, ysysl_num/ystep_den, ysysh_num/ystep_den


##---
def RatioError(v1,e1,v2,e2):
  #v1, e1: numerator value and error
  #v2, e2: denominator value and error  
  error1 = e1/v2
  error2 = (e2/v2)*(v1/v2)
  error = np.sqrt(error1*error1+error2*error2)
  return error
########################################################################################################################################################


# In[7]:


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


# In[34]:


########################################################################################################################################################
# Main Plot Code RAA
##--
import os

fig, axs = plt.subplots(1, 1, figsize=(6,4.5), sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})

## Log Scale for xy axes
axs.set_xscale('log')
# axs.set_yscale('log')

axs.set_xlim(x_min,x_max)
axs.set_ylim(y_min,y_max)

## horizontal line at y=1
#axs.axhline(1, color = "black", linewidth=0.2, alpha=0.5)

axs.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
axs.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))

## Label for xy axes
axs.set_xlabel(r'$p_{\mathrm{T}}~\mathrm{(GeV)}$', fontsize=18)
axs.set_ylabel(r'$R^{\mathrm{ch}}_{\mathrm{AA}}$', fontsize=16)



## Plot Theory Lines
for i, file in enumerate(files_aa_theory):
  x, y, xerr, yerr, xstep, ystep = GetRatioTheoryToTheory(file,file_pp_theory,error=True)
  # Steps
  axs.step(xstep, factor_theory[i]*ystep, where='post', color=color_theory[i],linestyle=line_theory[i])
  # Error bars for stat errors  
  axs.errorbar(x, factor_theory[i]*y, factor_theory[i]*yerr, marker="", linestyle="none", color=color_theory[i])


# Plot Exp Data Points
for i, file in enumerate(files_data):
  x, y, xerr, yerr, xstep, ystep, ysysl, ysysh = GetExp(file)
  # Markers with error bars for stat errors
  axs.errorbar(x, factor_theory[i]*y, yerr=factor_theory[i]*yerr, linestyle="none", color = color_data[i], marker=marker_data[i])
  # Shades for sys errors
  axs.fill_between(xstep, factor_theory[i]*(ystep-ysysl), factor_theory[i]*(ystep+ysysh), step='post', alpha=0.2, color = color_data[i])

## Legentds
handles = []
labels = []
for i, file in enumerate(files_data):
  dp = axs.errorbar(0, 0, 0, linestyle="none", color = color_data[i], marker=marker_data[i]) 
  dsys = axs.fill(np.NaN, np.NaN, alpha=0.2, color = color_data[i])
  handles.append((dp[0],dsys[0]))
  labels.append(label_data[i])

for i, file in enumerate(files_aa_theory):
  tl = axs.errorbar(0, 0, color=color_theory[i],linestyle=line_theory[i])
  handles.append(tl[0]) 
  labels.append(label_theory[i])   
axs.legend(handles,labels,ncol=1,loc='center left',edgecolor='none', frameon=True, facecolor='none', handletextpad=0.4, handleheight=1.8, labelspacing=0.2, bbox_to_anchor=(0, 0.62), borderaxespad=0.8, handlelength=1.6, fontsize=14)


## Text
axs.text(x_min+1, y_max-0.3, '$\mathrm{PbPb~(0\\text{-}10\%),~\sqrt{\it{s}_{NN}}=5.02~TeV}$\n'+ '$|\eta|<1.0$' , horizontalalignment='left', verticalalignment='bottom')
axs.text(x_max-80, y_min+0.08,'$\\hat{q}=\\hat{q}^{run}_{HTL}f(Q^{2})$\n $\\alpha^{fix}_{s}\!\!=\!\!0.3,Q_{sw}\!\!=\!\!\mathrm{2~GeV}$ \n $~\mathrm{JS(MATTER+LBT)}$',horizontalalignment='right', verticalalignment='bottom')
## Generate PDF File
plt.tight_layout()
plt.savefig('temp2')
## Add Logos
add_jetscape_logo('temp2',320,240,80,'_logo')
os.remove('temp.pdf')
RenameFileName='mv temp2_logo.pdf' + '  ' + OutputFilename
os.system(RenameFileName)
os.remove('temp2.pdf')
CommandOpenFile='open ' + OutputFilename
os.system(CommandOpenFile)
########################################################################################################################################################


# In[ ]:



