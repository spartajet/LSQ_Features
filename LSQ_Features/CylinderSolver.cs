using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics;

namespace LSQ_Features
{
    class CylinderSolver : Solver
    {
        double radius_result;

        List<double> init_paras(List<Point3D> datas)
        {
            List<double> paras = new List<double>();
            //the paras should contain x0,y0,z0,a,b,c & radius;
            double x0, y0, z0, a, b, c, radius;

            double k, p_A, p_B, p_C, p_D, p_E, p_F, p_G, p_H, p_I, p_J;
            //init the vector
            #region init vector
            DenseMatrix jacobian = new DenseMatrix(datas.Count, 9);
            Vector B = new DenseVector(datas.Count);
            for (int i = 0; i < datas.Count; ++i)
            {
                double x, y, z;
                Point3D temp = datas[i];
                x = temp.x;
                y = temp.y;
                z = temp.z;
                Vector<double> line = new DenseVector(9);
                line[0] = Math.Pow(y, 2);
                line[1] = Math.Pow(z, 2);
                line[2] = x * y;
                line[3] = x * z;
                line[4] = y * z;
                line[5] = x;
                line[6] = y;
                line[7] = z;
                line[8] = 1;
                jacobian.SetRow(i, line);

                B[i] = -Math.Pow(x, 2);
            }
            Vector<double> P = jacobian.QR(MathNet.Numerics.LinearAlgebra.Generic.Factorization.QRMethod.Thin).Solve(B);

            k = 2.0d / (1.0d + P[0] + P[1]);
            p_A = k;
            p_B = k * P[0];
            p_C = k * P[1];
            p_D = k * P[2];
            p_E = k * P[3];
            p_F = k * P[4];
            p_G = k * P[5];
            p_H = k * P[6];
            p_I = k * P[7];
            p_J = k * P[8];

            if (
               P[2].AlmostEqualInDecimalPlaces(0, 3) &&
               P[3].AlmostEqualInDecimalPlaces(0, 3) &&
               P[4].AlmostEqualInDecimalPlaces(0, 3)
               )
            {
                if (
                P[0].AlmostEqualInDecimalPlaces(1, 3))
                {
                    a = 0;
                    b = 0;
                    c = 1;
                }
                else if (P[1].AlmostEqualInDecimalPlaces(1, 3))
                {
                    a = 0;
                    b = 1;
                    c = 0;
                }
                else
                {
                    a = 1;
                    b = 0;
                    c = 0;
                }
            }
            else
            {

                double v_a, v_b, v_c;

                if (p_C < p_B && p_C < p_A)
                {
                    v_c = Math.Pow(1 - p_C, 0.5);
                    v_a = p_E / (-2 * v_c);
                    v_b = p_F / (-2 * v_c);
                }
                else
                {
                    if (p_A < p_B)
                    {
                        v_a = Math.Pow(1 - p_A, 0.5);
                        v_b = k * P[2] / (-2 * v_a);
                        v_c = k * P[3] / (-2 * v_a);
                    }
                    else
                    {
                        v_b = Math.Pow(1 - p_B, 0.5);
                        v_a = p_D / (-2 * v_b);
                        v_c = p_F / (-2 * v_b);
                    }
                }


                Vector<double> vec = (new DenseVector(new double[] { v_a, v_b, v_c })).Normalize(2.0);
                a = vec[0];
                b = vec[1];
                c = vec[2];
            }
            #endregion
            //init the zero position
            #region init the point
            DenseMatrix point_A = new DenseMatrix(4, 3, new double[]
            {
                -2*(b*b+c*c),   2*a*b,    2*a*c,
                2*a*b,    -2*(a*a+c*c),   2*b*c,
                2*a*c,    2*b*c,    -2*(a*a+b*b),
                a,         b,                 c
            }
            );
            Vector<double> point_B = new DenseVector(new double[]{
                    p_G,p_H,p_I,1*(a*datas[0].x+b*datas[0].y+c*datas[0].z)});

            DenseMatrix test_A = new DenseMatrix(3, 3, new double[]
            {
                -2*(b*b+c*c),   2*a*b,    2*a*c,
                2*a*b,    -2*(a*a+c*c),   2*b*c,
                //2*a*c,    2*b*c,    -2*(a*a+b*b),
                a,         b,                 c
            }
            );
            Vector<double> test_B = new DenseVector(new double[]{
                    p_G,p_H,1*(a*datas[0].x+b*datas[0].y+c*datas[0].z)});
            Vector<double> point_P;
            if (!c.AlmostEqual(1, 3))
            {
                point_P =
                 point_A.Transpose().Multiply(point_A).Cholesky().Solve(
                 point_A.Transpose().Multiply(point_B)
                 );
            }
            else
            {
                point_P = test_A.Transpose().Multiply(test_A).Cholesky().Solve(test_A.Transpose()).Multiply(test_B);
            }

            x0 = point_P[0];
            y0 = point_P[1];
            z0 = point_P[2];
            #endregion

            //init the radius
            #region init the radius
            radius = 0;
            foreach (var temp in datas)
            {
                double u = c * (temp.y - y0) - b * (temp.z - z0);
                double v = c * (temp.y - y0) - b * (temp.z - z0);
                double w = c * (temp.y - y0) - b * (temp.z - z0);
                radius += Math.Sqrt(u * u + v * v + w * w);
            }
            radius /= datas.Count;
            #endregion
            paras.Add(x0);
            paras.Add(y0);
            paras.Add(z0);
            paras.Add(a);
            paras.Add(b);
            paras.Add(c);
            paras.Add(radius);
            return paras;
        }

        public override void Estimate(List<Point3D> datas)
        {

            double c2, s2, c1, s1;
            double x0, y0, z0, a, b, c, radius;
            double deltavalue_old;
            double deltavalue_new;



            #region try to fit by gauss-newton method


            var current_para = init_paras(datas);
            deltavalue_old = getDeltaValue(current_para, datas);
            while (true)
            {
                x0 = current_para[0];
                y0 = current_para[1];
                z0 = current_para[2];
                a = current_para[3];
                b = current_para[4];
                c = current_para[5];
                radius = current_para[6];

                Vector<double> pos_old = new DenseVector(new double[]{
                        x0,y0,z0
                        });
                Vector<double> vec_old = new DenseVector(new double[]{
                        a,b,c
                        });
                //step I
                #region step 1

                Point3D zero_new = new Point3D(pos_old);
                List<Point3D> transfer_data = new List<Point3D>();
                #endregion
                //step II
                #region step 2
                if (a.AlmostEqual(1))
                {
                    s1 = 0;
                    c1 = 1;
                    s2 = -1;
                    c2 = 0;
                }
                else
                {
                    c1 = c / Math.Sqrt(Math.Pow(b, 2) + Math.Pow(c, 2));
                    s1 = -b / Math.Sqrt(Math.Pow(b, 2) + Math.Pow(c, 2));
                    c2 = (c * c1 - b * s1) / Math.Sqrt(
                        Math.Pow(a, 2) + Math.Pow(c * c1 - b * s1, 2)
                        );
                    s2 = -a / Math.Sqrt(
                        Math.Pow(a, 2) + Math.Pow(c * c1 - b * s1, 2)
                        );
                }
                DenseMatrix U1 = new DenseMatrix(3, 3, new double[]
                       {
                           c2,0,s2,
                           0,1,0,
                           -s2,0,c2,
                        }
                          );
                DenseMatrix U2 = new DenseMatrix(3, 3, new double[]
                       {
                           1,0,0,
                           0,c1,s1,
                           0,-s1,c1,
                       }
                );

                DenseMatrix U = U2 * U1;

                //move all the points to the new positions.
                foreach (Point3D temp in datas)
                {
                    transfer_data.Add((temp - zero_new) * U);
                }
                x0 = 0;
                y0 = 0;
                z0 = 0;
                a = 0;
                b = 0;
                c = 1;

                #endregion

                //step III
                //calculate the jacobian matrix & deviations
                #region step 3
                DenseMatrix J = new DenseMatrix(transfer_data.Count, 5);
                Vector<double> dev = new DenseVector(transfer_data.Count);

                foreach (Point3D temp in transfer_data)
                {
                    double u, v, w;
                    double x, y, z;
                    double r;
                    x = temp.x;
                    y = temp.y;
                    z = temp.z;

                    //in theory, but x0 = y0 = z0 =0
                    //          and  a = b = 0 , c=1
                    //u = c * (y - y0) - b * (z - z0);
                    //v = a * (z - z0) - c * (x - x0);
                    //w = b * (x - x0) - a * (y - y0);
                    u = y;
                    v = -x;
                    w = 0;
                    //double rr = Math.Sqrt(Math.Pow(u, 2) + Math.Pow(v, 2) + Math.Pow(w, 2))
                    r = Math.Sqrt(u * u + v * v);

                    Vector<double> single = new DenseVector(
                        new double[]{
                            -x/r, -y/r,-x*z/r,-y*z/r,-1
                        });
                    J.SetRow(transfer_data.IndexOf(temp), single);
                    dev[transfer_data.IndexOf(temp)] = -(r - radius);
                }
                #endregion

                //step IV
                #region step 4
                Vector<double> P =
               J.Transpose().Multiply(J).Cholesky().Solve(
               J.Transpose().Multiply(dev)
               );
                #endregion
                //step V
                #region step 5
                Vector<double> pos_off =
                    new DenseVector(
                        new double[]{
                            P[0],P[1],-P[0]*P[2]-P[1]*P[3]
                        }
                        );

                Vector<double> pos_new = pos_old + pos_off * U.Transpose();

                Vector<double> vec_off = new DenseVector(
                        new double[]{
                            P[2],P[3],1
                        }
                        );

                Vector<double> vec_new = vec_off * U.Transpose();

                double radius_new = radius + P[4];
                #endregion
                List<double> paras_new = new List<double>();
                paras_new.AddRange(pos_new);
                paras_new.AddRange(vec_new);
                paras_new.Add(radius_new);
                deltavalue_new = getDeltaValue(paras_new, datas);

                double value = Judge_out(paras_new, current_para);
                double value2 = Math.Abs(deltavalue_new - deltavalue_old);
                if (value < Solver.IterationEndValue || value2 < Solver.IterationEndValue)
                {
                    break;
                }
                else
                {
                    current_para = paras_new;
                    deltavalue_old = deltavalue_new;
                }
            }
            this.x = current_para[0];
            this.y = current_para[1];
            this.z = current_para[2];
            this.i = current_para[3];
            this.j = current_para[4];
            this.k = current_para[5];
            this.radius_result = current_para[6];

            #endregion
        }

        private double getDeltaValue(List<double> para, List<Point3D> datas)
        {
            double x0 = para[0];
            double y0 = para[1];
            double z0 = para[2];
            double a = para[3];
            double b = para[4];
            double c = para[5];
            double radius = para[6];
            double result = 0;
            foreach (Point3D temp in datas)
            {
                double u, v, w;
                double x, y, z;
                double r;
                x = temp.x;
                y = temp.y;
                z = temp.z;

                u = c * (y - y0) - b * (z - z0);
                v = a * (z - z0) - c * (x - x0);
                w = b * (x - x0) - a * (y - y0);
                //u = c * y - b * z;
                //v = a * z - c * x;
                //w = b * x - a * y;
                // r = Math.Sqrt(y * y + x * x);
                r = Math.Sqrt(Math.Pow(u, 2) + Math.Pow(v, 2) + Math.Pow(w, 2))
                    /
                    Math.Sqrt(Math.Pow(a, 2) + Math.Pow(b, 2) + Math.Pow(c, 2));
                result += r - radius;
            }
            return result / datas.Count;
        }

        private double Judge_out(List<double> paras_new, List<double> current_para)
        {
            DenseVector a = DenseVector.OfEnumerable(paras_new);
            DenseVector b = DenseVector.OfEnumerable(current_para);
            return (a - b).Norm(2);
        }
        public override string ToString()
        {
            string temp = "x:  " + x.ToString(Solver.Output_Format) + "\r\n" +
                           "y:  " + y.ToString(Solver.Output_Format) + "\r\n" +
                           "z:  " + z.ToString(Solver.Output_Format) + "\r\n" +
                           "a:  " + i.ToString(Solver.Output_Format) + "\r\n" +
                           "b:  " + j.ToString(Solver.Output_Format) + "\r\n" +
                           "c:  " + k.ToString(Solver.Output_Format) + "\r\n" +
                           "diameter:  " + (radius_result * 2).ToString(Solver.Output_Format);
            return temp;
        }
    }
}
