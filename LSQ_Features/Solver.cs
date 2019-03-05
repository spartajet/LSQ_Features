using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace LSQ_Features
{
    abstract class Solver
    {
        public static double IterationEndValue = 0.00001;
        public static string Output_Format = "F6";
        public double x, y, z, i, j, k;
        public static Solver FromType(SolverType solverType)
        {
            switch (solverType)
            {
                case SolverType.Line:
                    return (new LineSolver());
                case SolverType.Plane:
                    return (new PlaneSolver());
                case SolverType.Cylinder:
                    return (new CylinderSolver());
            }

            throw (new ApplicationException(String.Format("Unknown solver type: '{0}'.", solverType)));
        }
        public abstract void Estimate(List<Point3D> datas);

        public abstract override string ToString();
      
    }
}
