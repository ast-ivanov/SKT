using System;
using System.Collections.Generic;
using Mke;

namespace SKT
{
    struct Coordinate
    {
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
    }
    class Program
    {
        const int Kx = 10, Kz = 10, Ky = 1;
        const int K = Kx * Ky * Kz;
        const int n = 1000;
        const double x_left = 0, x_right = 100;
        const double z_bottom = 0, z_top = 100;
        const double hx = (x_right - x_left) / Kx, hy = 1, hz = (z_top - z_bottom) / Kz;
        const double mes = hx * hy * hz;
        const double alpha = 1e-2;
        static double[] p_pract = new double[K];
        static Coordinate[] Sources;
        static void Main(string[] args)
        {
#region Filling
            P_Filling();
            Sources_Filling();
            calc_g_pract();
            calcA();
            calcB();
#endregion

            double[,] I = new double[K, K];
            for (int i = 0; i < K; i++)
            {
                I[i, i] = alpha; 
            }
            
            
            SLAE.SumTwoMatrixes(I, A);
            double[] p = SLAE.Gauss(A, b);


            for (int i = 0; i < Kz; i++)
            {
                for (int j = 0; j < Kx; j++)
                {
                    Console.Write("{0:0.##}\t",p[i * Kz + j]);
                }
                Console.WriteLine();
            }
        }

        static void P_Filling()
        {
            // p_pract[(Kx / 2 - 1) + (Kz / 2 - 1) * Kx] = 10;
            // p_pract[(Kx / 2 - 1) + (Kz / 2) * Kx] = 10;
            // p_pract[(Kx / 2) + (Kz / 2 - 1) * Kx] = 10;
            // p_pract[(Kx / 2) + (Kz / 2) * Kx] = 10;
            p_pract[(Kx - 2) + 0 * Kx] = 10;
            p_pract[(Kx - 2) + 1 * Kx] = 10;
            p_pract[(Kx - 1) + 0 * Kx] = 10;
            p_pract[(Kx - 1) + 1 * Kx] = 10;
            for (int i = 0; i < Kz; i++)
            {
                for (int j = 0; j < Kx; j++)
                {
                    Console.Write("{0:0.##}\t",p_pract[i * Kz + j]);
                }
                Console.WriteLine();
            }
        }

        static double[] dg_pract = new double[n];
        static void calc_g_pract()
        {
            for (int i = 0; i < n; i++)
            {
                for (int k = 0; k < K; k++)
                {
                    dg_pract[i] += p_pract[k] * dg(k, i);
                }
            }
        }

        static void Sources_Filling()
        {
            double step = (x_right - x_left) / (n - 1);
            Sources = new Coordinate[n];
            for (int i = 0; i < n; i++)
            {
                Sources[i] = new Coordinate { X = i * step, Y = 0, Z = 0 };                
            }
        }

        static double r(int cell_number, int source_number)
        {
            double x = hx * (cell_number % Kx + 0.5);
            double y = hy * (cell_number % Ky + 0.5);
            double z = hz * (cell_number % Kz + 0.5);
            return Math.Sqrt(Math.Pow(x - Sources[source_number].X, 2) + Math.Pow(y - Sources[source_number].Y, 2) + Math.Pow(z - Sources[source_number].Z, 2));
        }
        static double dz(int cell_number, int source_number)
        {
            double dz = hz * (cell_number % Kz + 0.5);
            return dz - Sources[source_number].Z;
        }
        static double dg(int cell_number, int source_number)
        {
            return mes / (4 * Math.PI * Math.Pow(r(cell_number, source_number), 3)) * dz(cell_number, source_number);
        }

        static double[,] A = new double[K, K];
        static void calcA()
        {
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        A[i, j] += dg(i, k) * dg(j, k);
                    }
                }
            }
        }

        static double[] b = new double[K];
        static void calcB()
        {
            for (int k = 0; k < K; k++)
            {
                for (int i = 0; i < n; i++)
                {
                    b[k] += dg(k, i) * dg_pract[i];
                }
                
            }
        }

        static List<int>[] neighbourIndexes = new List<int>[K];//хранить индексы соседних ячеек (массив списков)
        static void fillNeighbourIndexes()
        {
            for(int j = 0; j < K; j++)
                neighbourIndexes[j] = new List<int>();
            int i = 0;
            for (int z = 0; z < Kz; z++)
            {
                for(int y = 0; y < Ky; y++)
                {
                    for (int x = 0; x < Kx; x++)
                    {
                        if (x != 0) neighbourIndexes[i].Add(i - 1);
                        if (x != Kx - 1) neighbourIndexes[i].Add(i + 1);
                        if (y != 0) neighbourIndexes[i].Add(i - Kx);
                        if (y != Ky - 1) neighbourIndexes[i].Add(i + Kx);
                        if (z != 0) neighbourIndexes[i].Add(i - Kx * Ky);
                        if (z != Kz - 1) neighbourIndexes[i].Add(i + Kx * Ky);
                        i++;
                    }
                }
            }
        }
    }
}
