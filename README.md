# Simulation-of-HL-ATLASexperiment
Being ATLAS one of the biggest collaborations and experiments when it come to particle physics research, the following upgrades being undertaken at CERN to reach High Luminosity conditions in the following years will generate a higher amount of data and not all of it can be stored. 


This increase in luminosity will produce 200 proton-proton collisions every 25ns, which translates to an outstanding 25 billion collisions a year, producing 5 times more Highs Bosons than nowadays, over 15 million a year after the upgrade is finished in 2025.  
 
This increase in the number of collisions will lead to a significant increase in the amount of data produced which leads to a data storage problem. Unfortunately, the experiments, such as ATLAS, can only record a 1000 of these collisions every second. 

In our duty as scientists, me and my team at university, we have developed a software simulating the conditions of the HL-LHC at ATLAS and implemented a “track trigger” algorithm that differentiates the common collisions (background) from the relevant ones (signal) that could lead to discover unknown predicted particles. 


- This repository contains the JAVA software developed to simulate the events from particles at High Luminosity and improved data-recording capabilities in real-time, improving detection and storage efficiency.


In order to improve the efficiency of the HLLHC, it was key to start by coding the scenario of the detectors and the generation of muons from the proton-proton collision before developing the track trigger algorithm. 

# Each part of the code focuses on the following described aspects: 
 
 ## • Particle.java 
 
   Initialising the properties of the particles to detect, muons in our simulation, the Particle.java class stores the type, position, momentum and direction of the particle. At the same time it creates a space in the code to store the trajectory, later on initialised in RunSimulation.java and creates space to record the reduction of energy from LossEnergy.java but without altering the direction. 
 
 
## •  EnergyLoss.java

  Implements a modified Bethe-Bloch equation to correct the energy as particles pass through matter. It is then stored in Particle.java, 
 
## •  MCS.java 

  Simulates high-energy collisions of composite protons through a Monte Carlo Simulation algorithm. The data for the calculations in these simulation steps is obtained from the properties of the experiment expressed in Geometry.java. 
 
## •  Geometry.java

  Generates a multiple scatter of particles on nuclei. This amount of generated scattering is then used to generate a gaussian-like distribution.  
 There is also a translation into code of the characteristics of the cylindrical detector which is expressed in polar coordinates. 
 
## •  ParticleTracker.java 
 
 This section of the software structure also implements the EnergyLoss.java class as well as the multiple scattering initialised on Geometry.java. 
 
 Following to that, there is an RK4 integrator to record the propagation of the particle generated, expressed in cartesian. 
 
 Finally, this ParticleTracker.java class finds the E (Electric) and B (magnetic) fields at coordinates x, y, z and relativistic speed (t, x, y, z). 
 
## •  RunSimulation.java 

  Presenting the main part of the software developed, it is in charge of running the simulation produced for 1000 events, generating particles in all random directions at 2500 MeV. 
 Records the trajectory of those interesting, high perpendicular momentum, particles between the tracking layers 9A & 9B of the detector set. 
 Accordingly to the results obtained at this stage, we implement an algorithm stating an smearing correction of 0.9 to bring the data closer to a real-life (non-ideal) event, and consequently finds the perpendicular momentum of the particles. 
 
## •  Histogram.java 

  This final class used in our software structure records all the results produced from the RunSimulation.java class and prints them in .csv files ready for analysis. 
 
 
