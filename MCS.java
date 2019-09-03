import java.io.*;
import java.util.Scanner;

public class MCS
{

   public double K = 0.307075;

   // iron details

   //public double rho = 7.87;
   //public double Z=26;
   //public double A=55.845;
   public double em = 0.511;




   public double getX0 (double rho, double Z, double A)
   {
     if (rho==0 || Z==0 || A==0){
           double bb=0;
           return bb;
       }
       else  {double x0 = ((716.4*A)/((rho*Z*(Z+1)*(Math.log((287/(Math.sqrt(Z))))))))/100;
               return x0;
     }
  }

   public double getTheta0 (Particle particle, double rho, double Z, double A)
   {
       if (rho==0 || Z==0 || A==0){
           double bb=0;
           return bb;
       }
       else  {
           double x = 0.01;
           int charge =particle.getCharge();
           double Beta = particle.getBeta();
           double p = particle.getLastMomentumMag();
           double I = 0.0000135*Z;
           double x0 = ((716.4*A)/((rho*Z*(Z+1)*(Math.log((287/(Math.sqrt(Z))))))))/100;
           double thasd = -(13.6*charge*(Math.sqrt(x/x0)*(1+(0.038*Math.log(x/x0)))))/(Beta*p);
           return thasd;
   }
  }

}