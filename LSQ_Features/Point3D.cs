using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;

namespace LSQ_Features
{
    class Point3D
    {
        public double x
        {
            get
            {
                return pos[0];
            }
        }
        public double y
        {
            get
            {
                return pos[1];
            }
        }
        public double z
        {
            get
            {
                return pos[2];
            }
        }

        internal Vector<double> _pos;
        public Vector<double> pos
        {
            get
            {
                return _pos.SubVector(0,3);
            }
        }
      
        public Point3D(double[] input)
        {
            _pos = new DenseVector(
                new double[] { input[0], input[1], input[2] ,1}
                );
        }
        public Point3D(Vector<double> v)
            :this(v.ToArray())
        {
        }
        public Point3D(string str)
            : this(str
                .Split(' ')
                .Select(n => Convert.ToDouble(n))
                .ToArray())
        {
        }
        public static Point3D operator -(Point3D a,Point3D b)
        {
            return new Point3D(a.pos.Subtract(b.pos));
        }

        public static Point3D operator *(Point3D a, DenseMatrix b)
        {
            return new Point3D(a.pos * b);
        }
        public void Transform(DenseMatrix trans)
        {
            this._pos = (this.pos * trans);
        }
        public static Point3D operator *(DenseMatrix a, Point3D b)
        {
            return new Point3D(a * b.pos);
        }
    }
}
