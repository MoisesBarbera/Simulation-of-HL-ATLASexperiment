import java.io.*;
import java.util.Scanner;

public class EnergyLoss
{

   public double K = 0.307075;

   // iron details

   //public double rho = 7.87;
   //public double Z=26;
   //public double A=55.845;
   public double em = 0.511;




   public double getEnergyLoss (Particle particle, double rho, double Z, double A)
   {
       if (rho==0 || Z==0 || A==0){
           double bb=0;
           return bb;
       }
       else  {
             int charge =particle.getCharge();
             double Mass =particle.getMass();


             double Beta = particle.getBeta();
             double gamma = particle.getGamma();
             double I = 0.0000135*Z;
             double w=((2*em*Beta*Beta*gamma*gamma)/(1+((2*gamma*em)/(Mass+((em/Mass)*(em/Mass))))));

             double bb=(100)*( K*charge*charge*rho*(Z/A)*(1/(Beta*Beta))*((0.5*Math.log((2*em*Beta*Beta*gamma*gamma*w)/(I*I)))-(Beta*Beta)));
             return bb;
       }
   }

}