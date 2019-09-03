import java.io.*;
import java.util.Scanner;
import java.util.Random;
import java.util.Arrays;
import java.util.Collections; 

class RunSimulation
{
    static Random randGen = new Random();
    static Scanner keyboard = new Scanner(System.in);

    // maximum allowed number of steps before simulation is aborted
    static final int numSteps = 1000;
    
    static int verbose = 0;
    
    public static Particle[] GetParticles(double[] pos0, double startMomentum)
    {
    // example to simulate just one muon starting at p0
    // with a total momentum startMomentum and theta=startAngle
    // we follow the particle physics "convention"
    // to have the z-axis in the (approximate) direction of the beam
    
       Particle [] Particles_gen = new Particle[1];
    
       Random randGen = new Random();
       
      double a2 = randGen.nextDouble();
      double Phi=a2*(2*Math.PI);
      
      
       double px=2500*Math.cos(Phi); // change total momentum
       double py=2500*Math.sin(Phi);
       double pz=0;
       
       //double px = 2500;
       //double py=2500;
       //double pz=0; 
    
    
       double [] p0 = {px,py,pz};
    
       Particles_gen[0] = new Particle("muon", pos0, p0, numSteps);

       return Particles_gen;
    }
    
    private static double ndetector(double sigma)
    {
        Random randGen = new Random();
        double detector = sigma*randGen.nextGaussian();
        
        return detector;
    }
    
    
     
    
    private static void Verbose (String x){
    if (verbose == 1){
        System.out.println(x);
    }
    
    }
    
    public static void main (String [] args ) throws IOException
     {
      // set the seed for the random number generator so it always
      // produces the same sequence of random numbers
      randGen.setSeed(7894694);
      
       //double px = 2500;
      //double py=2500;
       //double pz=0; 
    
      
      double startMomentum = 2500;
      Random randGen = new Random();
     
      /*System.out.println("Type in total (maximum) simulation time s [e.g. 1E-9]"); //change this 
      double simTime = keyboard.nextDouble(); */
      
      double simTime = 1e-8;
      
      double detectorLength = 6.2; // in m
 
      int numberOfEvents = 1000;

      // Define the genotrical properties of the experiment
      Geometry Experiment = SetupExperiment(detectorLength);

      // initial time and x,y,z position (m) of particles
      double[] pos0 = {0., 0., 0., 0.};

      // setup histograms for analysis
      
      Histogram phi9a = new Histogram(50, -0.1, 0.5, "phi9a");//change this 
      Histogram phi9b = new Histogram(50, -0.1, 0.5, "phi9b");//change this
      
     
      Histogram PP = new Histogram (100, startMomentum*0.9,startMomentum*1, "perpendicular momentum");
      Histogram PP12 = new Histogram (100, startMomentum*0.9,startMomentum*1, "perpendicular momentum at 1 and 2");
      
      //smearing code
      
      Histogram PP2 = new Histogram (100, startMomentum*0.9869,startMomentum*0.9878, "perpendicular momentum with smearing");
      
      double array[];
      array = new double [numberOfEvents];
      double array2[];
      array2 = new double [numberOfEvents];
      
     
      // start of main loop: run the simulation numberOfEvents times
      for (int nev = 0; nev < numberOfEvents; nev++) {

        if (nev % 1000 == 0) {
          //System.out.println("Simulating event " + nev);
        }

        // get the particles of the event to simulate
        Particle [] Particles_gen = GetParticles(pos0, startMomentum);
        
        // simulate propagation of each generated particle,
        // store output in Particles_sim
        Particle [] Particles_sim = new Particle[Particles_gen.length];
        
        for (int ip = 0; ip < Particles_gen.length; ip++) {
          // create an instance of the particle tracker
          // usage: particleTracker(particle, time to track (s), number of time steps)
          particleTracker tracker = new particleTracker(Particles_gen[ip], simTime, numSteps);

          // some output (need to disable before running large numbers of events!)
          //System.out.println("Simulating particle " + ip + " of event " + nev);
          //Particles_gen[ip].print();

          Particles_sim[ip] = tracker.track(Experiment);
        }
        // end of simulated particle propagation

        
        // write scatter plot for event 0, particle to disk into file "sample.csv"
        /*if (nev == 0) {
          Particles_sim[0].DumpTXYZPxPyPz("sample.csv");
        }*/
        
        // simulate detection of each particle in each element
        Particle [] Particles_det = Experiment.detectParticles(Particles_sim);

        // this is just for dumping the simulation to the screen
        for (int ip = 0; ip < Particles_sim.length; ip++) {
            for (int idet = 1; idet < Experiment.getNshapes(); idet++) {
              double [] txyz = Particles_det[ip].getPosition(idet);
                //System.out.println("Particle " + ip + " detection in volume " + idet);
                //System.out.println("[t,x,y,z] = [" + txyz[0] + ", " + txyz[1] + ", " +
                          // txyz[2] + ", " + txyz[3] + "]");
            }
        }
        // end of simulated particle detection


        
        // analyse event (= 1 particle) and fill histograms
        // this is to show, how to access the variables
        // and do some calculations with the results of the simulation

        // retrieve initial generated particle momentum
        double [] gen_mom = Particles_gen[0].getLastMomentum();
        double gen_mommag = Math.sqrt(gen_mom[0]*gen_mom[0]+gen_mom[1]*gen_mom[1]
                      +gen_mom[2]*gen_mom[2]);

        // retrieve simulated particle position and momentum at the end of the simulation
        double [] sim_mom = Particles_sim[0].getLastMomentum();
        double sim_mommag = Math.sqrt(sim_mom[0]*sim_mom[0]+sim_mom[1]*sim_mom[1]
                      +sim_mom[2]*sim_mom[2]);
        double [] sim_pos = Particles_sim[0].getLastPosition();
        
        // fill histograms
        //gen_startMom.fill(gen_mommag);
        //sim_endMom.fill(sim_mommag);

        // calculate theta angles, defined here as angle w.r.t. the z-plane
        // theta ~ acos(Z/sqrt(X^2+Y^2+Z^2))

        // generated - this should be equal to startAngle given at the start
        double gen_Theta = Math.acos(gen_mom[2]/gen_mommag);
        //gen_startTheta.fill(gen_Theta);

        // after simulation - muon will have scattered around a bit
        double sim_Theta = Math.acos(sim_mom[2]/sim_mommag);
        //sim_endTheta.fill(sim_Theta);

        // after detection: uses initial position (assumed to be known with pos0) and detected positions
        // the first detector 9A has volume number 10
        double [] det_pos = Particles_det[0].getPosition(10);
        double det1Theta = Math.acos((det_pos[3]-pos0[3])/
                     Math.sqrt(Math.pow(det_pos[1]-pos0[1], 2)
                           +Math.pow(det_pos[2]-pos0[2], 2)
                           +Math.pow(det_pos[3]-pos0[3], 2)));
                        
        double [] det_mom = Particles_det[0].getMomentum(10); //px, py, pz
        double phi9A = Math.atan2(det_pos[2],(det_pos[1]));
        
        
        
        phi9a.fill(phi9A);

        // the second detector 9B has volume number 11
        double [] det_pos2 = Particles_det[0].getPosition(11);
        double det2Theta = Math.acos((det_pos[3]-pos0[3])/
                     Math.sqrt(Math.pow(det_pos[1]-pos0[1], 2)
                           +Math.pow(det_pos[2]-pos0[2], 2)
                           +Math.pow(det_pos[3]-pos0[3], 2)));
                           
        double [] det_mom2 = Particles_det[0].getMomentum(11); //px, py, pz
        double phi9B = Math.atan2(det_pos2[2],(det_pos2[1]));
        
        double R9A = 0.9;
        double R9B = 0.91;
        
        double delta = Math.atan(((phi9B-phi9A)*R9B)/(R9B-R9A));
        
        double PerpP = (0.3*4*0.905/(2*delta) )*1000;
        
        //smeared
        
        double phi9A2 = Math.atan2(det_pos[2]+ ndetector(0.0000001),(det_pos[1]+ ndetector(0.0000001)));
        
        double phi9B2 = Math.atan2(det_pos2[2]+ ndetector(0.0000001),(det_pos2[1]+ ndetector(0.0000001)));
        
        double delta2 = Math.atan(((phi9B2-phi9A2)*R9B)/(R9B-R9A));
        
        double PerpP2 = (0.3*4*0.905/(2*delta2) )*1000;
        
        array[nev] = PerpP;
        
        
        //track trigger
        if (PerpP >= 2468.2) { 
        //Particles_sim[0].DumpTXYZPxPyPz("sample_" +nev+ ".csv");
        
        double arrayx[];
        arrayx = new double [8];
        double arrayy[];
        arrayy = new double [8];
        
        double [] det_posUno = Particles_det[0].getPosition(2);
        double [] det_posDos = Particles_det[0].getPosition(3);
        double [] det_posTres = Particles_det[0].getPosition(4);
        double [] det_posCuatro= Particles_det[0].getPosition(5);
        double [] det_posCinco = Particles_det[0].getPosition(6);
        double [] det_posSeis = Particles_det[0].getPosition(7);
        double [] det_posSiete = Particles_det[0].getPosition(8);
        double [] det_posOcho = Particles_det[0].getPosition(9);
        
        
        arrayx[0]=det_posUno[1];
        arrayy[0]=det_posUno[2];
        
        arrayx[1]=det_posDos[1];
        arrayy[1]=det_posDos[2];
        
        arrayx[2]=det_posTres[1];
        arrayy[2]=det_posTres[2];
        
        arrayx[3]=det_posCuatro[1];
        arrayy[3]=det_posCuatro[2];
        
        arrayx[4]=det_posCinco[1];
        arrayy[4]=det_posCinco[2];
        
        arrayx[5]=det_posSeis[1];
        arrayy[5]=det_posSeis[2];
        
        arrayx[6]=det_posSiete[1];
        arrayy[6]=det_posSiete[2];
        
        arrayx[7]=det_posOcho[1];
        arrayy[7]=det_posOcho[2];
        
        FileWriter file = new FileWriter("tracked particle"+nev+".csv");     // this creates the file with the given name
        PrintWriter outputFile = new PrintWriter(file);
        
        // now make a loop to write the contents of each bin to disk, one number at a time
        for (int n = 0; n < 8; n++) {
            
            outputFile.println(arrayx[n] + "," + arrayy[n]);
        }
        outputFile.close();
        
        
        }
        
        
        
        //Verbose(PerpP);
        
        
        
        array2[nev] = Math.pow(PerpP-startMomentum,2);
        
        PP.fill(PerpP);
        
        PP2.fill(PerpP2);
        
        
        phi9b.fill(phi9B);
        
        // end of analysis
        
        
        
            

      }
      // end of main event loop
      //processs to calculate the standard deviation
     double sum = 0; // initialize sum 
         
      // Iterate through all elements and add them to sum 
      for (int nev = 0; nev < numberOfEvents; nev++) {
            sum +=  array2[nev]; 
      } 
       
      double sgm = Math.sqrt(sum/numberOfEvents);
      System.out.println("standard deviation is " + sgm); 
      
      
      
     
     
      PP.writeToDisk("perpendicularmomentum.csv");
      
      PP2.writeToDisk("perpendicularmomentum2.csv");
      phi9a.writeToDisk("phi9a.csv");
      phi9b.writeToDisk("phi9b.csv");
      System.out.println("Execution Complete");
     }

    
    public static Geometry SetupExperiment (double detectorLength)
    {
        /*
         example setup the experiment:
         As stated below:
         Assume that on their way particles encounter 
         –the wall of a Berilium beam pipe (r = 35mm, thickness = 3 mm)
         –Silicon detector layer 1 (r = 45 mm,  thickness = 0.5 mm)
         –Silicon detector layer 2 (r = 80 mm ,  thickness = 0.5 mm)
         –Silicon detector layer 3 (r = 120 mm ,  thickness = 0.5 mm)
         –Silicon detector layer 4 (r = 180 mm ,  thickness = 0.5 mm)
         –Silicon detector layer 5 (r = 300 mm ,  thickness = 0.5 mm)
         –Silicon detector layer 6 (r = 400 mm ,  thickness = 0.5 mm)
         –Silicon detector layer 7 (r = 500 mm ,  thickness = 0.5 mm)
         –Silicon detector layer 8 (r = 700 mm ,  thickness = 0.5 mm)
         –Silicon trigger layer 9A (r = 900 mm ,  thickness = 0.5 mm)
         –Silicon trigger layer 9B (r = 910 mm ,  thickness = 0.5 mm)
         first object defines the world 
         then adds the beryllium beam pipe 
         then all 10 detectors 
         These measurements have been converted to m
        */
        
    Geometry Experiment = new Geometry(randGen, 0.0005);
    
    Experiment.AddCylinder(0, -detectorLength/2,  //world
                 0.9105, detectorLength/2,
                 0., 0., 0.);
                 
    //may need a vacuum inside the beam pipe to be defined as our first medium
                 
    Experiment.AddCylinder(0.035, -detectorLength/2,  // beryllium beam pipe
                 0.038, detectorLength/2,
                1.85, 4, 9.012); // density in g/cm3, 4 is Be atomic number
    
    Experiment.AddCylinder(0.045, -detectorLength/2, //DETECTOR 1
                 0.0455,  detectorLength/2,
                 2.33, 14, 28.085); // density in g/cm3, 14 is Si atomic number
    
    Experiment.AddCylinder( 0.080, -detectorLength/2, //DETECTOR 2
                 0.0805, detectorLength/2,
                 2.33, 14, 28.085);
                 
    Experiment.AddCylinder(0.120, -detectorLength/2,  //DETECTOR 3
                 0.1205, detectorLength/2,
                 2.33, 14, 28.085); 
                 
    Experiment.AddCylinder(0.180, -detectorLength/2,  //DETECTOR 4
                 0.1805, detectorLength/2,
                 2.33, 14, 28.085);
    
    Experiment.AddCylinder(0.300, -detectorLength/2,  //DETECTOR 5
                 0.3005, detectorLength/2,
                 2.33, 14, 28.085); 
                 
        Experiment.AddCylinder(0.400, -detectorLength/2,  //DETECTOR 6
                 0.4005, detectorLength/2,
                 2.33, 14, 28.085);
                 
    Experiment.AddCylinder(0.500, -detectorLength/2,  //DETECTOR 7
                 0.5005, detectorLength/2,
                 2.33, 14, 28.085);
                 
    Experiment.AddCylinder(0.700, -detectorLength/2,  //DETECTOR 8
                 0.7005, detectorLength/2,
                 2.33, 14, 28.085);
                 
    Experiment.AddCylinder(0.900, -detectorLength/2,  //DETECTOR 9A
                 0.9005, detectorLength/2,
                 2.33, 14, 28.085);
                 
    Experiment.AddCylinder(0.910, -detectorLength/2,  //DETECTOR 9B
                 0.9105, detectorLength/2,
                 2.33, 14, 28.085);
                 
    //Experiment.Print();

    return Experiment;
    }

}