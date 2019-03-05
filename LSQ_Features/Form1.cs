using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using System.Diagnostics;

namespace LSQ_Features
{
    public partial class Form1 : Form
    {
        SolverType type;
        Solver solver;
        Stopwatch timer;
        bool status_points
        {
            get;
            set;

        }
        List<Point3D> points;
        public Form1()
        {
            InitializeComponent();
            status_points = false;
            timer = new Stopwatch();
            points = new List<Point3D>();
            comboBox1.SelectedIndex = 2;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            status_points = false;
            points.Clear();
            var abc = new OpenFileDialog();
            abc.InitialDirectory = System.AppDomain.CurrentDomain.BaseDirectory;
            abc.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
            abc.FilterIndex = 2;
            abc.RestoreDirectory = true;
            if (abc.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    if (File.Exists(abc.FileName))
                    {
                        StreamReader myStream = new StreamReader(abc.FileName);
                        string buf;
                        while ((buf = myStream.ReadLine()) != null)
                        {
                            points.Add(new Point3D(buf));
                        }
                        if (points.Count > 0)
                        {
                            status_points = true;
                        }
                    }
                    else
                    {
                        throw (new Exception(" file didn't exist"));
                    }

                }
                catch (Exception ex)
                {
                    MessageBox.Show("Error: Could not read file from disk. Original error: " + ex.Message);
                }
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            if (!status_points)
            {
                MessageBox.Show("Please read a Point files first");
                return;
            }
            timer.Restart();
            solver.Estimate(points);
            timer.Stop();
            label1.Text = solver.ToString() +"\r\ntime cost:   " 
                + timer.ElapsedMilliseconds.ToString();
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            switch (comboBox1.SelectedIndex)
            {
                case 0:
                    type = SolverType.Line;
                    break;
                case 1:
                    type = SolverType.Plane;
                    break;
                case 2:
                    type = SolverType.Cylinder;
                    break;
                default:
                    break;
            }
            solver = Solver.FromType(type);
        }
    }
}
