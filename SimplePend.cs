//============================================================================
// SimplePend.cs Defines a class for simulation simple pendulum
//============================================================================

using System;

namespace Sim
{
    public class SimplePend
    {
        private double len = 1.1; // pendulum length
        private double g = 9.81; // gravitational field strength
        int n = 2; // number of states
        private double[] x; // array of states
        private double[] f; // right side of equation evaluated

        //--------------------------------------------------------------------
        // constructor
        //--------------------------------------------------------------------
        public SimplePend()
        {
            x = new double[n];
            f = new double[n];
            
            x[0] = 1.0;
            x[1] = 0.0;
        }
        
        //------------------------------------------------------------------------
        // step: perform one integration step via Euler's method
        //------------------------------------------------------------------------
        public void step(double dt)
        {
            rhsFunc(x,f);

            int i;
            for (i=0;i<n;++i)
            {
                x[i] = x[i] + f[i]*dt;
            }
        }

        //------------------------------------------------------------------------
        // Fourth-Order Runge-Kutta Method
        //------------------------------------------------------------------------
        public void runge_kutta(double dt)
        {
            rhsFunc(x,f);

            double[,] k = new double[2,4]; // k values for theta and omega
            
            k[0,0] = x[1];
            k[0,1] = x[1]+0.5*k[0,0]*dt;
            k[0,2] = x[1]+0.5*k[0,1]*dt;
            k[0,3] = x[1]+k[0,2]*dt;
            k[1,0] = x[1];
            k[1,1] = x[1]+0.5*k[0,0]*dt;
            k[1,2] = x[1]+0.5*k[0,1]*dt;
            k[1,3] = x[1]+k[0,2]*dt;

            int i;
            for (i=0;i<n;++i)
            {
                x[i] = x[i]+(1.0/6.0)*(k[i,0]+2.0*k[i,1]+2.0*k[i,2]+k[i,3])*dt;
            }
        }

        //------------------------------------------------------------------------
        // rhsFunc: function to calculate rhs of pendulum ODEs
        //------------------------------------------------------------------------
        public void rhsFunc(double[] st, double[] ff)
        {
            ff[0] = st[1];
            ff[1] = -(g/len)*Math.Sin(st[0]);
        }

        //------------------------------------------------------------------------
        // Getters and Setters
        //------------------------------------------------------------------------
        public double L
        {
            get {return(len);}

            set
            {
                if(value > 0.0)
                    len = value;
            }
        }

        public double G
        {
            get {return(g);}
            
            set
            {
                if(value >= 0.0)
                    g = value;
            }
        }

        public double theta
        {
            get {return x[0];}

            set {x[0] = value;}
        }

        public double thetaDot
        {
            get {return x[1];}

            set {x[1] = value;}
        }
    }
}