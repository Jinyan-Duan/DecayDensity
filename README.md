
     =========================================================
     Geant4 - an Object-Oriented Toolkit for Simulation in HEP
     =========================================================

                            Hadr07
                            ------

   Survey energy deposition and particle's flux from an hadronic cascade.
   Use PhysicsConstructor objects rather than predefined G4 PhysicsLists.
   Show how to plot a depth dose profile in a rectangular box.    

	
 1- MATERIALS AND GEOMETRY DEFINITION

  The geometry consists of a stack of one or several blocks of homogenous
  material, called absorbers.

  A minimum of 4 parameters define the geometry :
     - the number of absorbers (NbOfAbsor)	
     - the material of each absorber,
     - the thickness of each absorber,
     - the tranverse dimension of the stack (sizeYZ)

  In addition a transverse uniform magnetic field can be applied.
      eg: /globalField/setValue 0 0 5 tesla

  The absorber is surrounded by a World volume (vacuum)

  A function, and its associated UI command, allows to build a material
  directly from a single isotope.

  The default geometry is built in DetectorConstruction, but the above parameters 
  can be changed interactively via commands defined in DetectorMessenger.

  To be identified by the ThermalScattering module, the elements composing a
  material must have a specific name (see G4ParticleHPThermalScatteringNames.cc)
  Examples of such materials are build in Hadr06/src/DetectorConstruction.
 	
 2- PHYSICS LIST
   
  "Full" set of physics processes are registered, but via PhysicsConstructor
  objects rather than complete pre-defined G4 physics lists. This alternative 
  way gives more freedom to register physics.
  
  Physics constructors are either constructors provided in Geant4 (with G4 prefix)
  or 'local'. They include : HadronElastic, HadronInelastic, IonsInelastic, GammaNuclear,
  RadioactiveDecay and Electomagnetic.
  (see geant4/source/physics_lists/constructors)

  HadronElasticPhysicsHP include a model for thermalized neutrons, under the control of a command
  defined in NeutronHPMesseger.

  GammmaNuclearPhysics is a subset of G4BertiniElectroNuclearBuilder.

  ElectromagneticPhysics is a simplified version of G4EmStandardPhysics.

  Several hadronic physics options are controlled by environment variables.
  To select them, see Hadr07.cc

 3- AN EVENT : THE PRIMARY GENERATOR
 
  The primary kinematic consists of a single particle starting at the
  left face of the box. The type of the particle and its energy are set 
  in the PrimaryGeneratorAction class, and can be changed via the G4 
  build-in commands of G4ParticleGun class (see the macros provided with 
  this example).

  In addition one can choose randomly the impact point of the incident
  particle. The corresponding interactive command is built in
  PrimaryGeneratorMessenger class.

 A RUN is a set of events.
 
 4- PHYSICS

  The program computes the energy deposited in each absorber,
  and the flux of particles emerging in the world.
  Processes invoked and particles generated are listed.

 5- HISTOGRAMS
         
  The test has several built-in 1D histograms, which are managed by
  G4AnalysisManager and its Messenger. The histos can be individually 
  activated with the command :
  /analysis/h1/set id nbBins  valMin valMax unit 
  where unit is the desired unit for the histo (MeV or keV, etc..)
  (see the macros xxxx.mac).

            1     "total energy deposited in absorber 1
            2     "total energy deposited in absorber 2
            ...........................................
            9     "total energy deposited in absorber 9
            10    "Edep (MeV/mm) profile along absorbers"
   
  One can control the name of the histograms file with the command:
  /analysis/setFileName  name  (default Hadr07)

  It is possible to choose the format of the histogram file : root (default),
  xml, csv, by using namespace in HistoManager.hh

  It is also possible to print selected histograms on an ascii file:
  /analysis/h1/setAscii id
  All selected histos will be written on a file name.ascii (default Hadr07) 
  
 6- TRACKING and STEP MAX
 
  Hadr07 computes the distribution of energy deposited along the trajectory of 
  the incident particle : the so-called longitudinal energy profile,
  or depth dose distribution (histogram 10).
  The energy deposited (edep) is randomly distribued along the step (see
  SteppingAction).
  
  In order to control the accuracy of the deposition, the maximum  step size 
  of charged particles is computed automatically from the binning of 
  histogram 10.
     
  As an example, this limitation is implemented as a 'full' process :
  see StepMax class and its Messenger. The 'StepMax process' is registered
  in the Physics List, via a physicsConstructor object (a builder).
     
  StepMax is evaluated in the StepMax process. 
  A boolean UI command allows to deactivate this mechanism.
  Another UI command allows to define directly a stepMax value.
 

 7- VISUALIZATION
 
   The Visualization Manager is set in the main().
   The initialisation of the drawing is done via the commands
   /vis/... in the macro vis.mac. To get visualisation:
   > /control/execute vis.mac
 	
   The tracks are drawn at the end of event, and erased at the end of run.   
   gamma green   
   neutron yellow
   negative particles (e-, ...) red
   positive particles (e+, ions, ...) blue
 	
 8- HOW TO START ?
 
   Execute Hadr07 in 'batch' mode from macro files :
 	% Hadr07   run1.mac
 		
   Execute Hadr07 in 'interactive mode' with visualization :
 	% Hadr07
	Idle> control/execute vis.mac
 	....
 	Idle> type your commands
 	....
 	Idle> exit
	
 Macros provided in this example:
  - Na22.mac: multilayers. Radioactive source
  - alpha.mac: alpha (400 MeV). Limit the step size from histo 10
  - ionC12.mac: C12 (2.4 GeV). Limit the step size from histo 10
  - water.mac: e- (4 MeV) in Water
    
 Macros to be run interactively:
  - proton.mac: proton (1 GeV). Multilayers
  - vis.mac: To activate visualization
