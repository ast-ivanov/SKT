using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Mke
{
    delegate ref double Matrix(int i, int j);
    class SLAE
    {
        public static double[] Gauss(double[,] matrix, double[] b)
        {
            int i, j, k;
            double sum;
            int n = matrix.GetLength(0);
            double[] RightPart = new double[n];
            for (int ii = 0; ii < n; ii++)
                RightPart[ii] = b[ii];

            double t;
            double[] x = new double[n];

            double max, d;
            double str;
            int m, imax = 0;
            for (k = 0; k < n - 1; k++)
            {
                max = 0;
                for (m = k; m < n; m++)
                {
                    if (Math.Abs(matrix[m, k]) > Math.Abs(max)) { max = matrix[m, k]; imax = m; }
                }
                if (max == 0) Console.WriteLine("деление на 0");
                else
                {
                    for (i = 0; i < n; i++)
                    {
                        str = matrix[k, i]; matrix[k, i] = matrix[imax, i]; matrix[imax, i] = str;
                    }
                    d = RightPart[k]; RightPart[k] = RightPart[imax]; RightPart[imax] = d;
                }

                for (i = k + 1; i < n; i++)
                {
                    t = matrix[i, k] / matrix[k, k];
                    RightPart[i] = RightPart[i] - t * RightPart[k];
                    for (j = k + 1; j < n; j++)
                    {
                        matrix[i, j] = matrix[i, j] - t * matrix[k, j];
                    }
                }

            }

            x[n - 1] = RightPart[n - 1] / matrix[n - 1, n - 1];
            for (k = n - 2; k >= 0; k--)
            {
                sum = 0;
                for (j = k + 1; j < n; j++) sum += matrix[k, j] * x[j];
                x[k] = (RightPart[k] - sum) / matrix[k, k];
            }

            return x;
        }

        public static void SumTwoMatrixes(double[,] A1, double[,] A2)
        {
            var n = A1.GetLength(0);
            for(int i=0; i<n; i++)
            {
                for(int j = 0; j<n; j++)
                {
                    A2[i, j] += A1[i, j];
                }
            }
        }
    }
}