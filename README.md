# LArSoftDumper
LArSoft analyzer to dump data to train neural networks into a ROOT tree

The analyzer currently dumps:

- Images from 2 different sources:
 
  - TPC (electrons, charged particles) Wire planes, in the side laterals. (using hits extracted from raw waveforms)
  
  - Optical system (photons), form the optical readouts, behind the wire planes(using ophits extracted from raw waveforms).

- Interaction identification, in current studies, we consider(values->interaction dictionary in interaction.png):
 
  - Signal: Charge Current Quasi Elastic - CCQE (W+- boson interactions), labeled by a value of 1001
 
  - Background: everything else (for example, Neutral-CQS)

- Simulated photons arriving at Optical channels as a vector v(number of channels) with the total number of photons in the event for each channel.

We only consider photons from the TPB coated PMTs, sensible to VUV(direct) and Visible(reflected by cathode in the middle of the volume) photons.
