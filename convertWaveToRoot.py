import readTrc
import matplotlib.pyplot as plt

import ROOT

from ROOT import TFile, TTree
from array import array

import sys
import os.path


input_folder  = "./data/";
#output_folder = "./root_files/NewAmplifier_MPTR/S12572-015C/";
output_folder = "./root_files/";
  
  
  
file_header = "_S12572_3015C_"
#NPULSES = 1;

      

#if (sys.argc>1):
input_folder = sys.argv[1]
file_header = sys.argv[2]
caption = sys.argv[3]


print "input folder: " + sys.argv[1]
print "file header: " + sys.argv[2]
print "caption: " + sys.argv[3]

#if (sys.argc>2) NPULSES = sys.argv[2]

  
#std::cout << "file_header: " << file_header << std::endl;
#std::cout << "NPULSES to be processed: " << NPULSES << std::endl;
output_file_name = output_folder + caption + ".root"
output_file = TFile(output_file_name, "RECREATE")



n = array('i', [ 0 ] )
waveTree = TTree("waveTree", "waveTree");

maxN = 80000
index = array('i', [ 0 ])

#t_time_l = array('f', maxN*[ 0. ] )
#t_amp_l  = array('f', maxN*[ 0. ] )
#t_time = array('f', maxN*[ 0. ] )
#t_amp  = array('f', maxN*[ 0. ] )

t_time_l = ROOT.std.vector(float)()
t_amp_l  = ROOT.std.vector(float)()
t_time   = ROOT.std.vector(float)()
t_amp    = ROOT.std.vector(float)()

waveTree._t_time_l = t_time_l
waveTree._t_amp_l = t_amp_l
waveTree._t_time = t_time
waveTree._t_amp = t_amp

waveTree.Branch("index", index, 'index/I' )
#waveTree.Branch("t_time_l", t_time_l, 't_time_l[80000]/F' )
#waveTree.Branch("t_amp_l", t_amp_l, 't_amp_l[80000]/F' )
#waveTree.Branch("t_time", t_time, 't_time[80000]/F' )
#waveTree.Branch("t_amp", t_amp, 't_amp[80000]/F' )

waveTree.Branch("t_time_l", t_time_l )
waveTree.Branch("t_amp_l", t_amp_l )
waveTree.Branch("t_time", t_time)
waveTree.Branch("t_amp", t_amp)


NWAVEFORMS = 20000

for iWave in range (0, NWAVEFORMS):
    
    #fName = "C3-cross-talk-00000.trc"
    fName_l = "%s%s%s%05d.trc" % (input_folder, "C2", file_header, iWave)
    fName   = "%s%s%s%05d.trc" % (input_folder, "C3", file_header, iWave)
    #fName   = "%s%s%s%05d.trc" % (input_folder, "C2", file_header, iWave)
    index[0] = iWave
    
    #print fName;
    #reading SiPM waveforms
    if os.path.isfile(fName_l) and os.path.isfile(fName) :
    
        print "reading SiPM and LASER waveforms (" + str(index) + "): " + fName + " :: " + fName_l
        datX_l, datY_l, m = readTrc.readTrc( fName_l )
        datX, datY, m = readTrc.readTrc( fName )
    
        for i in range(len(datX)):
            
            t_time_l.push_back(datX_l[i])
            t_amp_l.push_back(datY_l[i])
            t_time.push_back(datX[i])
            t_amp.push_back(datY[i])
    
    elif os.path.isfile(fName) :
    
        print "reading only SiPM waveform (" + str(index) + "): " + fName
        datX, datY, m = readTrc.readTrc( fName )
    
        for i in range(len(datX)):
            
            t_time_l.push_back(0.)
            t_amp_l.push_back(0.)
            t_time.push_back(datX[i])
            t_amp.push_back(datY[i])
            
    
    #plt.plot(t_time, t_amp)
    #plt.show()
    
    
    
    
    waveTree.Fill()
    
    t_time_l.clear()
    t_amp_l.clear()
    t_time.clear()
    t_amp.clear()
 
    
    
output_file.Write()
output_file.Close()




