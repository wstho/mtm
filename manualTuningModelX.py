
'''
############################################################################################################################################################
updated Aug. 2020: stimulation to match Atherton et al. 2005 protocol (figure 2A). WST

############################################################################################################################################################
     
'''
from Neuron_model_extended import NeuronModel
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from modelSetUp import packageModel
from neuron import h
import bluepyopt.ephys as ephys


modelFile = "modelParameter_update_28_11.txt"
morphologyFile = "morphology/Mailly2003F4Cclean.swc"
parameterFile = "parameters.json"
mechanismsFile = "mechanisms.json"

#read trace from paper
comp = pd.read_csv("/Users/wsthompson/Desktop/Karolinska/SNr/Traces/Atherton2005_2B1.csv", header = None)
comp.rename(columns={0: "x", 1: "y"}, inplace = True, errors = 'raise')
comp.sort_values('x', axis = 0, inplace = True)

#align x-axis
x0 = comp.loc[0,'x']
comp.x -= x0
x_sm = comp.x
y_sm = comp.y


packageModel(modelFile) # Packaging the manual tuning file modelParameter.txt into the two files for NeuroModel: parameters.json and mechanisms.json

infoCell = dict()

modelCell=NeuronModel(param_file=parameterFile,morph_file=morphologyFile,mech_file=mechanismsFile,cell_name="TEST")

simulator=ephys.simulators.NrnSimulator(cvode_active=False) # In here is neuron.h., simulator.neuron.h = neuron.h. in common NEURON+python syntax

modelCell.instantiate(sim=simulator)  

infoCell.update({"Cell": modelCell})

infoCell.update({"simulator": simulator})

infoCell.update({"synapses":list()})
    
vSave = simulator.neuron.h.Vector()

for isec, sec in enumerate(modelCell.icell.soma):

	 for seg in sec:
                  vSave.record(getattr(seg,'_ref_v'))
                  infoCell.update({"soma_access":seg})
               
                  spikeTime = simulator.neuron.h.Vector()
                  recordingSpikingActivity = simulator.neuron.h.NetCon(getattr(seg,'_ref_v'),None, sec = sec)
                  recordingSpikingActivity.threshold = 0
                  recordingSpikingActivity.record(spikeTime)
                                         
                  infoCell.update({"spike_con": recordingSpikingActivity})
                  infoCell.update({"spike_train": spikeTime})

                  infoCell.update({"soma_voltage":vSave})
 


for sec in simulator.neuron.h.allsec():
    
    for seg in sec:
             '''
             Here we are iterating through the whole morphology, which is divided into compartments (sec) and compartments are made up of segments (seg).
             We can place current clamp, voltage clamp and receptors on the segments to stimulate the neuron

             For optimisation of the neuron we use current clamps in the soma and stimulate with protocols consisting of step currents __|¨¨|__ 
             with a holding current in pA and a stimulation current in pA. E.g. -300 pA during whole protocol (3000 ms)  and 100pA (1000 ms)

             
             '''
           
            
             if "soma" in sec.name():
                      '''
                      If you want to place point processes on the soma. You can place it here

                      
                      '''

             elif "dend" in sec.name():
	             
                      '''
                      If you want to place point processes on the dendrite. You can place it here.

                      By using print(dir(seg)) you can see what is availble to do on the current seg
                      
                      For example this is the output I got:

                      'area', 'ca_ch', 'ca_ion', 'cm', 'diam', 'k_ion', 'kir2_ch', 'na2_ch', 'na_ion', 'node_index', 'pas', 'point_processes', 'ri', 'sec', 'v', 'volume', 'x'
                      If you for example add print(seg.point_processes()) to the code, you can see what point processes you have placed here.

                      You can try to add synapses only on distal dendrites, you can get information on the x,y,z. Soma is not at origo, need to correct for this!
                      sec.x3d(i),
                      sec.y3d(i),
                      sec.z3d(i),


                      '''
                      '''
                      expSynapse=simulator.neuron.h.Exp2Syn(seg)
                      expSynapse.tau1=2 # You can change the parameter of the synapse - we have models for this for AMPA and NMDA receptor and GABA_A receptors, but this is a simpler model
                      expSynapse.tau2=10 # you can also change syn.e, which is the reversal time
                      
                      ExtrinsicSpikeTimes = [500,600,700,800]  # The stimulation times only have to be part of a list (in ms)
                      VecStim = simulator.neuron.h.VecStim()
                      convertedVector = simulator.neuron.h.Vector(ExtrinsicSpikeTimes)
                      VecStim.play(convertedVector)
                      
                      ncToSynapse = simulator.neuron.h.NetCon(VecStim,expSynapse)
                      ncToSynapse.delay=1
                      ncToSynapse.threshold=0
                      ncToSynapse.weight[0]=0.1 # You can determine the weight
		      '''
                      '''
                      dendriteResponse = simulator.neuron.h.Vector()
                      dendriteResponse.record(getattr(seg,'_ref_v'))
                      infoCell["synapses"].append(dendriteResponse)
                      '''
                      
tSave = simulator.neuron.h.Vector()
tSave.record(simulator.neuron.h._ref_t)

simulator.neuron.h.tstop=1000
simulator.neuron.h.run()


cellVoltage=infoCell["soma_voltage"]  
arr = np.array(cellVoltage)
minV = np.amin(arr)
  
plt.figure(0)
plt.plot(tSave,cellVoltage,label="Simulated",c='black')
plt.xlabel("Time (ms)")
plt.ylabel("Membrane potential (mV)")
plt.title("SNr GABAergic Neuron (Spontaneous)")
plt.plot(x_sm, y_sm, c ='red', alpha = 0.6, label = 'Atherton 2005')
plt.axhline(y=minV, c = 'black', alpha = 0.5, dashes=[5, 2])
plt.legend(loc = 1)
plt.savefig('Test.png')
plt.show()

for dendriteLocation in infoCell["synapses"]:
         plt.figure(1)
         plt.plot(tSave,dendriteLocation)

plt.show()
                      
                      
                      
