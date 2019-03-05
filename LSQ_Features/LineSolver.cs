using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.LinearAlgebra.Double.Factorization;

namespace LSQ_Features
{
    class LineSolver :Solver
    {
        public override void Estimate(List<Point3D> datas)
        {
            double sum_x = 0;
            double sum_y = 0;
            double sum_z = 0;
            foreach (Point3D temp in datas)
            {
                sum_x += temp.x;
                sum_y += temp.y;
                sum_z += temp.z;
            }
            sum_x /= datas.Count;
            sum_y /= datas.Count;
            sum_z /= datas.Count;

            DenseMatrix jacobian = new DenseMatrix(datas.Count, 3);
            foreach (Point3D temp in datas)
            {
                Vector<double> gradient = new DenseVector(3);
                gradient[0] = temp.x - sum_x;
                gradient[1] = temp.y - sum_y;
                gradient[2] = temp.z - sum_z;
                jacobian.SetRow(datas.IndexOf(temp), gradient);
            }
            Svd svd = jacobian.Svd(true);
            // get matrix of left singular vectors with first n columns of U
            Matrix<double> U1 = svd.U().SubMatrix(0, datas.Count, 0, 3);
            // get matrix of singular values
            Matrix<double> S = new DiagonalMatrix(3, 3, svd.S().ToArray());
            // get matrix of right singular vectors
            Matrix<double> V = svd.VT().Transpose();

            Vector<double> parameters = new DenseVector(3);
            parameters = V.Column(0);
            x = sum_x;
            y = sum_y;
            z = sum_z;
            i = parameters[0];
            j = parameters[1];
            k = parameters[2];
        }
        public override string ToString()
        {
            string temp = "x:  " + x.ToString(Solver.Output_Format) + "\r\n" +
                "y:  " + y.ToString(Solver.Output_Format) + "\r\n" +
                "z:  " + z.ToString(Solver.Output_Format) + "\r\n" +
                "a:  " + i.ToString(Solver.Output_Format) + "\r\n" +
                "b:  " + j.ToString(Solver.Output_Format) + "\r\n" +
                "c:  " + k.ToString(Solver.Output_Format);
            return temp;
        }
    }
}
