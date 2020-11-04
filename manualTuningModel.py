
'''
############################################################################################################################################################

The model is defined by its morphology file and model file (which contains the ion channel models, their placements and conductances.)

Installing Python+NEURON - https://neuron.yale.edu/neuron/getstd

Check out the python file,Neuron_model_extended.py, in the folder to see how the files are loaded and 
https://github.com/BlueBrain/BluePyOpt/tree/master/bluepyopt/ephys (if you are interested in how it is loaded into the simulator, more specifically the mechanisms.py, morphology.py and parameters.py)

Neuron Tutorial
https://neuron.yale.edu/neuron/static/docs/neuronpython/index.html

                            
############################################################################################################################################################
     
'''
from Neuron_model_extended import NeuronModel
import numpy as np
import json
import matplotlib.pyplot as plt
from modelSetUp import packageModel
from neuron import h
import bluepyopt.ephys as ephys


modelFile = "modelParameter.txt"
morphologyFile = "morphology/MTC180800A-IDB-cor-rep.swc"
parameterFile = "parameters.json"
mechanismsFile = "mechanisms.json"

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
                      IClampHolding = simulator.neuron.h.IClamp(seg)
                      IClampHolding.delay = 0
                      IClampHolding.dur = 2000
                      IClampHolding.amp = 0.2 # nA, but experiments usually use pA

                      IClampStimulation = simulator.neuron.h.IClamp(seg)
                      IClampStimulation.delay = 200
                      IClampStimulation.dur = 1000
                      IClampStimulation.amp = 0.05 # nA, but experiments usually use pA
                  
                      
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
                      dendriteResponse = simulator.neuron.h.Vector()
                      dendriteResponse.record(getattr(seg,'_ref_v'))
                      infoCell["synapses"].append(dendriteResponse)
                      '''
                      
tSave = simulator.neuron.h.Vector()
tSave.record(simulator.neuron.h._ref_t)

simulator.neuron.h.tstop=2000
simulator.neuron.h.run()


cellVoltage=infoCell["soma_voltage"]    
plt.figure(0)
plt.plot(tSave,cellVoltage,label="membrane potential",c='black')
plt.xlabel("Time (ms)")
plt.ylabel("Membrane potential (mV)")
plt.title("TestCell Simulation")

for dendriteLocation in infoCell["synapses"]:
         plt.figure(1)
         plt.plot(tSave,dendriteLocation)

plt.show()
                      
                      
                      
