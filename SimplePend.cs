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
            double[] k1 = new double[n]; // k1 values for theta and omega
            double[] k2 = new double[n]; // k2 values for theta and omega
            double[] k3 = new double[n]; // k3 values for theta and omega
            double[] k4 = new double[n]; // k4 values for theta and omega
            double[] mid1 = new double[n];
            double[] mid2 = new double[n];
            double[] end = new double[n];
            
            rhsFunc(x,f);
            k1 = f;

            int i;
            for (i=0;i<n;++i)
            {
                mid1[i] = x[i]+0.5*k1[i]*dt;
            }
            rhsFunc(mid1,k2);
            for (i=0;i<n;++i)
            {
                mid2[i] = x[i]+0.5*k2[i]*dt;
            }
            rhsFunc(mid2,k3);
            for (i=0;i<n;++i)
            {
                end[i] = x[i]+k3[i]*dt;
            }
            rhsFunc(end,k4);

            for (i=0;i<n;++i)
            {
                x[i] = x[i]+(1.0/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*dt;
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