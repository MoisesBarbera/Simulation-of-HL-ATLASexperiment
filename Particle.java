import java.io.*;

class Particle
{
    // a simple class to store properties of a particle as well as its history
    
    //Particle specific constants
    private String ptype; // type of particle: currently "electron" or "proton" or "muon"
    private double m0;    // rest mass (MeV)
    private int charge;   // particle charge (in units of elementary charge)
    
    // Fourposition and fourmomenta: from start to end
    private double [][] position;
    private double [][] momentum;
    private int laststep; // how many steps the particle has made

    // constructor to initialise the particle
    // and make space to store the trajectory
    public Particle(String particleType, double [] position0, double [] momentum0, int maxSteps)
    {
        position = new double[maxSteps+1][4];
        momentum = new double[maxSteps+1][3];
        laststep = 0;
        position[0] = position0;
        momentum[0] = momentum0;
        ptype = particleType;
	if (ptype == "electron") {
	    m0 = 0.51099895;
	    charge = -1;
	} else if (ptype == "muon") {
	    m0 = 105.658;
	    charge = -1;
	} else if (ptype == "proton") {
	    m0 = 938.27231;
	    charge = 1;
	}
    }

    // return properties of particle
    public int getCharge() {return charge;}
    public double getMass() {return m0;}
    public int getSteps() {return laststep;}
    public String getType() {return ptype;}

    // Lorentz calculations for beta and gamma
    public double getBeta() { return getLastMomentumMag()/getLastE(); }
    public double getGamma() { return getLastE()/getMass(); }

    // calculate distance in space of last step
    public double getLastDistance()
    {
	double [] last = getLastPosition();
	double [] previous = getPosition(laststep-1);
	double dist = 0.;
	for (int i = 1; i < 4; i++) {
	    double delta = last[i]-previous[i];
	    dist += delta*delta;
	}
	return Math.sqrt(dist);
    }


    // set the next position+momentum
    public void updateParticle(double[] state)
    {
        laststep = laststep+1;
        for(int i=0; i<4; i++){
            position[laststep][i] = state[i];
        }
        for(int i=0; i<3; i++){
            momentum[laststep][i] = state[i+4];
        }
    }

    public void undoLastStep()
    {
        laststep = laststep-1;
    }

    // return the last position and momentum stored
    public double[] getPosition(int step) {return position[step];}
    public double[] getLastPosition() {return position[laststep];}

    public double[] getMomentum(int step) {return momentum[step];}
    public double[] getLastMomentum() {return momentum[laststep];}

    // return components of last position and momentum stored
    public double getLastT() {return position[laststep][0];}
    public double getLastX() {return position[laststep][1];}
    public double getLastY() {return position[laststep][2];}
    public double getLastZ() {return position[laststep][3];}
    public double getLastPx() {return momentum[laststep][0];}
    public double getLastPy() {return momentum[laststep][1];}
    public double getLastPz() {return momentum[laststep][2];}

    // return position and momentum of step i
    public double[] getState(int i)
    {
	double[] txyz_pxpypz = new double[7];
	for(int j=0; j<4; j++){
            txyz_pxpypz[j] = position[i][j];
	}
	for(int j=0; j<3; j++){
            txyz_pxpypz[j+4] = momentum[i][j];
	}
	return txyz_pxpypz;
    }
    public double[] getLastState()
    {
	return getState(laststep);
    } 
    
    public double getMomentumMag(int step)
    {
	double[] txyz_pxpypz = getState(step);
	double p = 0.;
	for (int j = 0; j < 3; j++) {
	    p += txyz_pxpypz[j+4]*txyz_pxpypz[j+4];
	}
	p = Math.sqrt(p);
	return p;
    }    
    public double getLastMomentumMag()
    {
	return getMomentumMag(laststep);
    }

    public double getE(int step)
    {
	double p = getMomentumMag(step);
	return Math.sqrt(m0*m0+p*p);
    }    

    public double getLastE()
    {
	return getE(laststep);
    }


    // the following methods return one big array
    // with all positions/momenta
    public double[] getAllPosMom(int coo) {
	double [] all = new double[laststep];
	for (int istep = 0; istep < laststep; istep++) {
	    if (coo < 4)
		all[istep] = position[istep][coo];
	    else if (coo < 7)
		all[istep] = momentum[istep][coo-4];
	    else if (coo == 7)
		all[istep] = getMomentumMag(istep);
	    else if (coo == 8)
		all[istep] = getE(istep);
	}
	return all;
    }
    public double[] getAllT() {return getAllPosMom(0);}
    public double[] getAllX() {return getAllPosMom(1);}
    public double[] getAllY() {return getAllPosMom(2);}
    public double[] getAllZ() {return getAllPosMom(3);}
    public double[] getAllPx() {return getAllPosMom(4);}
    public double[] getAllPy() {return getAllPosMom(5);}
    public double[] getAllPz() {return getAllPosMom(6);}
    public double[] getAllP() {return getAllPosMom(7);}
    public double[] getAllE() {return getAllPosMom(8);}

    public void print()
    {
        System.out.println("E = " + getLastE()
			   + " MeV, px = " + getLastPx()
			   + " MeV, py = " + getLastPy()
			   + " MeV, pz = " + getLastPz()
			   + " MeV, mass = " + getMass() + " MeV"
			   );
        System.out.println("t = " + getLastT()
			   + " s, x = " + getLastX()
			   + " m, y = " + getLastY()
			   + " m, z = " + getLastZ()
			   + " m, charge = " + getCharge()
			   );
    }
    
    public void applySmallRotation(double dtheta_xz, double dtheta_yz)
    {
        // calculate new momentum direction
        // approximates this into two sequential x-z and y-z changes
        double p, theta;
	
	double[] mom = getLastMomentum();
	
        p = Math.sqrt(mom[0]*mom[0] + mom[2]*mom[2]);
        theta = Math.atan2(mom[0], mom[2]) + dtheta_xz;
        momentum[laststep][0] = p*Math.sin(theta);
        momentum[laststep][2] = p*Math.cos(theta);

        p = Math.sqrt(mom[1]*mom[1] + mom[2]*mom[2]);
        theta = Math.atan2(mom[1], mom[2]) + dtheta_yz;
        momentum[laststep][1] = p*Math.sin(theta);
        momentum[laststep][2] = p*Math.cos(theta);
    }
    
    public void reduceEnergy(double Eloss)
    {
        // reduce energy of particle while keeping direction the same
	double E = getLastE();
        if (Eloss > E) {
	    momentum[laststep][0] = 0.;
	    momentum[laststep][1] = 0.;
	    momentum[laststep][2] = 0.;
        } else {
	    double pnew = Math.sqrt(Math.pow(E-Eloss, 2) - Math.pow(getMass(),2));
	    double factor = pnew/getLastMomentumMag();
	    momentum[laststep][0] *= factor;
	    momentum[laststep][1] *= factor;
	    momentum[laststep][2] *= factor;
        }
    }

    public void DumpTXYZPxPyPz(String filename) throws IOException
    {
        // write time and x,y,z coordinates as well as px,py,pz to CSV file
        FileWriter file = new FileWriter(filename);     // this creates the file with the given name
        PrintWriter outputFile = new PrintWriter(file); // this sends the output to file1

        // now make a loop to write the contents of each bin to disk, one number at a time
        for (int n = 0; n < getSteps()+1; n++) {
            double[] Y = getState(n);
            // comma separated values
            outputFile.println(Y[0] + "," + Y[1] + "," + Y[2] + "," + Y[3] + ","
			       + Y[4] + "," + Y[5] + "," + Y[6]);
        }
        outputFile.close(); // close the output file
        return;
    }

    
}