import java.io.*;
import java.util.Random;

class Geometry
{
    // Simple class to define the Geometry of an experiment
    // the "type" is a number defining the basic geometric object
    // the "rho_Z_A" array stores material information: 
    //    density, atomic number (nucleus charge), mass number (protons+neutrons)
    // the "shape" array stores the geometrical information
    // * cuboid: type=1,
    //   shape = 6 values representing two corners of the internal diagonal (x0,y0,z0), (x1,y1,z1)
    // every basic object can be identified by a number
    // the first object is taken as to be the "world",
    //    i.e. we can use it to abort the simulation
    // maxShapes is the maximum allowed number of basic objects, extend if needed

    // the class also impelments thee additional helper functions:
    // * doEloss -- calculate EnergyLoss via associated class and apply to Particle
    // * doMultScatter -- calculate MultipleScattering via associated class and apply to Particle
    // * detectParticle -- simulate a detected particle position with resolution smearing applied
    
    private Random randGen;
    private static final int maxShapes = 100;
    private int nshapes;
    private int [] type;
    private double [][] rho_Z_A;
    private double [][] shapes;
    
    private EnergyLoss [] Eloss;
    private MCS [] MultScatter;

    private double minfeaturesize;
    
    public Geometry(Random rg, double featuresize)
    {
        randGen = rg;
        nshapes = 0;
        type = new int[maxShapes];
        rho_Z_A = new double[maxShapes][3];
        shapes = new double[maxShapes][];
    
        Eloss = new EnergyLoss[maxShapes];
        MultScatter = new MCS[maxShapes];

        minfeaturesize = featuresize;
    }

    public int getNshapes() { return nshapes; }
     

    public int AddCylinder(double r0,double z0, 
                           double r1, double z1,
                           double rho,double Z, double A)
             {
                 if (nshapes >= maxShapes) {
                     return -1;
        }        
        // in cylindrical: x axis is rcos(theta), y axis is rsin(theta) z axis is normal
        type[nshapes] = 1;
        shapes[nshapes] = new double[4];
        shapes[nshapes][0] = r0;
        shapes[nshapes][1] = z0;
        shapes[nshapes][2] = r1;
        shapes[nshapes][3] = z1;    
        
        rho_Z_A[nshapes][0] = rho;
        rho_Z_A[nshapes][1] = Z;
        rho_Z_A[nshapes][2] = A;
                  
        Eloss[nshapes] = new EnergyLoss();
        MultScatter[nshapes] = new MCS();

        nshapes++;
        return (nshapes-1);
     }
    
    public void Print() // change this
    {
        System.out.println("stored " + getNshapes() + " objects.");
      for (int i = 0; i < nshapes; i++) {
        if (i == 0) {
            System.out.println("Maximum size of experiment given by object 0:");
        }
        if (type[i] == 1) {
            System.out.println("Geometry object #" + i + " = cylinder.");
            System.out.printf("   corners (%f, %f) - (%f, %f)%n",
                  shapes[i][0], shapes[i][1], shapes[i][2],
                  shapes[i][3]);
        }
        System.out.printf("   material rho = %.3f g/cm^3, Z = %f, A = %f%n",
        rho_Z_A[i][0], rho_Z_A[i][1], rho_Z_A[i][2]);
      }
      System.out.println("When scanning for volume transitions, the smallest feature size discovered will be " + minfeaturesize + "m.");
    }

    public boolean isInVolume(double x, double y, double z, int id)
    {
      // test if point (x,y,z) is in volume with identifier id
      // abort if being asked for an undefined volume
      
      double r=Math.sqrt(x*x+y*y);
      
      //need to change for a cylinder
      if (id >= getNshapes()) {
        return false;
      }
    
      if (type[id] == 1) {
        // cylinder
        return ( shapes[id][0] <= r
             && shapes[id][1] <= z
             && r <= shapes[id][2] 
             && z <= shapes[id][3]
             );
      }
      return false;
    }
    
   
    public int getVolume(double x, double y, double z)
    {
      // cycle through volumes in opposite order
      for (int i = getNshapes()-1; i >= 0; i--) {
        if (isInVolume(x, y, z, i)) {
        return i;
        }
      }
      // if we arrived here, we are outside everything
      return -1;
     } 

    public boolean isInVolume(Particle p, int id)
    {
      // test if particle p is currently in volume with identifier id
      double [] txyz = p.getLastPosition(); // read the current position: t,x,y,z
      return isInVolume(txyz[1], txyz[2], txyz[3], id);
    }
    
    public int getVolume(Particle p)
    {
      // get the highest volume number the particle is in
      double [] txyz = p.getLastPosition(); // read the current position: t,x,y,z
      return getVolume(txyz[1], txyz[2], txyz[3]);
    }
    
    public void doEloss(Particle p)
    {
      // CHECK AND COMPLETE HERE

      // dist gives the length of the last step
      double dist = p.getLastDistance();
      // determine which volume we are in
      int volume = getVolume(p);
      
      
      
      if (volume >= 0) {
    
      double lostE = Eloss[volume].getEnergyLoss(p, rho_Z_A[volume][0],rho_Z_A[volume][1],rho_Z_A[volume][2])*dist;
      if (lostE > p.getLastE() - p.getMass()) {
          lostE = p.getLastE() - p.getMass();
      }
      p.reduceEnergy(lostE);
      }
    
      /*
      if (volume >= 0) {
        // use Eloss[volume] to get dE/dx and dist to calculate dE
        double lostE =dist*Eloss[volume].getEnergyLoss(p, rho_Z_A[volume][0],rho_Z_A[volume][1],rho_Z_A[volume][2]);
        p.reduceEnergy(lostE);
      } */
    }
    
    public void doMultScatter(Particle p)
    {
      
      // use MultScatter[volume] to determine theta0
      // draw two Gaussian Random numbers
      int volume = getVolume(p);
       if (volume >= 0) {
        double theta0 = MultScatter[volume].getTheta0(p,rho_Z_A[volume][0],rho_Z_A[volume][1],rho_Z_A[volume][2]);
        
        double gausx = theta0*randGen.nextGaussian();
        double gausy = theta0*randGen.nextGaussian();

        p.applySmallRotation(gausx, gausy);
      }
     }

    public Particle[] detectParticles(Particle [] psim)
    {
      Particle [] pdet = new Particle[psim.length];
      for (int ip = 0; ip < psim.length; ip++) {
        pdet[ip] = new Particle(psim[ip].getType(), psim[ip].getPosition(0),
                    psim[ip].getMomentum(0), getNshapes());
        for (int idet = 1; idet < getNshapes(); idet++) {
          
          // average over all simulated particle positions in each volume
          int ncross = 0;
          double [] txyz_pxpypz = new double[7];
    
          for (int npoint = 0; npoint <= psim[ip].getSteps(); npoint++) {
            double [] pos = psim[ip].getPosition(npoint);
            if (getVolume(pos[1], pos[2], pos[3]) == idet) {
              double [] mom = psim[ip].getMomentum(npoint);
              for (int i = 0; i < 4; i++) {
                txyz_pxpypz[i] += pos[i];
              }
              for (int i = 0; i < 3; i++) {
                txyz_pxpypz[i+4] += mom[i];
              }
              ncross++;
            }
          }
          if (ncross > 0) {
            for (int i = 0; i < 7; i++) {
                txyz_pxpypz[i] /= ncross;
            }
          }
          pdet[ip].updateParticle(txyz_pxpypz);
        }
      }
      return pdet;
    }

    public double scanVolumeChange(Particle p, int lastVolume)
    {
        // scan the straight line between last and new position, if we missed a
         // feature of the experiment
         double dist = p.getLastDistance();
         int nsteps = (int) (dist/minfeaturesize) + 1;

         if (nsteps <= 1) {
        // step was small enough, continue
        return 1.;
       }
    
       double [] pos = p.getPosition(p.getSteps()-1);
       double [] end = p.getLastPosition();
       double [] delta = new double[4];
       for (int i = 0; i < 4; i++) {
        delta[i] = (end[i]-pos[i])/nsteps;
       }
       for (int n = 1; n <= nsteps; n++) {
        if (getVolume(pos[1]+n*delta[1], pos[2]+n*delta[2], pos[3]+n*delta[3]) != lastVolume) {
            if (n == 1) {
            return (n-0.5) / nsteps;
            } else {
            return (n-1.0) / nsteps;
           }
         }
       }
       return 1.;
    }
}