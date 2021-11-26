# LArSoftDumper
LArSoft analyzer to dump data to train neural networks into a ROOT tree

For each event, the analyzer currently dumps:

- Images from 2 different sources:
 
  - 1 from TPC (electrons, charged particles) Wire planes, in the side laterals. (using hits extracted from raw waveforms)
  
  - 1 from Optical system (photons), form the optical readouts, behind the wire planes(using ophits extracted from raw waveforms).

- A number representing the neutrino interaction, in current studies, we consider(values->interaction dictionary in interaction.png):
 
  - Signal: Charge Current Quasi Elastic - CCQE (W+- boson interactions), labeled by a value of 1001
 
  - Background: everything else (for example, Neutral-CQS)

- A vector v(number of channels) with the total number of simulated photons in the event arriving at each Optical channel.

We only consider photons detected by the TPB coated PMTs, sensible to VUV(direct) and Visible(reflected by cathode in the middle of the volume) photons.
For the Optical System Image, as SBND has 2 different volumes/optical&charge planes, in the current approach, we add opposing optical channels.

![Alt text](images/sbnd_tpc.png?raw=true "Title")