using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.Numerics;

namespace Numeric
{
	/// <summary>
	/// ルンゲクッタ法
	/// </summary>
	class RungeKutta
	{
		/// <summary>
		/// 時間ステップ幅とステップ数を元にRungeKuttaクラスを初期化する
		/// </summary>
		/// <param name="timestep">ステップ幅</param>
		/// <param name="element">ステップ数</param>
		public RungeKutta(double timestep, int element)
		{
			Element = element;
			TimeStep = timestep;
			Eqnum = 0;
			time = new double[Element];
			for (int i = 0; i < Element; i++)
			{
				time[i] = i * TimeStep;
			}
			Eq = new List<SystemEq>();
		}

		public int Element { get; private set; }	// 時間節点の数

		public delegate double SystemEq(double t, double[] x);

		/// <summary>
		/// 微分方程式を登録する
		/// </summary>
		/// <param name="CompEq">微分方程式</param>
		public void RegistrateMethod(SystemEq eq)
		{
			Eq.Add(eq);
			Eqnum++;
		}

		/// <summary>
		/// 微分方程式を消去する
		/// </summary>
		public void ClearRegistratedMethod()
		{
			Eq.Clear();
			Eqnum = 0;
		}

		/// <summary>
		/// セットされた連立常微分方程式を初期条件を元に解く
		/// </summary>
		/// <param name="initial_conditions">初期条件ベクトル</param>
		public void ODE4(params double[] initial_conditions)
		{
			x = new double[Eqnum, Element];
			try
			{
				setInitialCondition(initial_conditions);
			}
			catch (ArgumentException) { throw; }

			unsafe
			{
				double* k1 = stackalloc double[Eqnum];
				double* k2 = stackalloc double[Eqnum];
				double* k3 = stackalloc double[Eqnum];
				double* k4 = stackalloc double[Eqnum];
				double[] temp_x = new double[Eqnum];

				int loop = Element - 1;
				int eqnum = Eqnum;
				for (int i = 0; i < loop; ++i)
				{
					for (int j = 0; j < eqnum; ++j) { temp_x[j] = x[j, i]; }
					for (int j = 0; j < eqnum; ++j) { k1[j] = TimeStep * Eq[j](time[i], temp_x); }

					for (int j = 0; j < eqnum; ++j) { temp_x[j] = x[j, i] + k1[j] * 0.5; }
					for (int j = 0; j < eqnum; ++j) { k2[j] = TimeStep * Eq[j](time[i] + TimeStep * 0.5, temp_x); }

					for (int j = 0; j < eqnum; ++j) { temp_x[j] = x[j, i] + k2[j] * 0.5; }
					for (int j = 0; j < eqnum; ++j) { k3[j] = TimeStep * Eq[j](time[i] + TimeStep * 0.5, temp_x); }

					for (int j = 0; j < eqnum; ++j) { temp_x[j] = x[j, i] + k3[j]; }
					for (int j = 0; j < eqnum; ++j) { k4[j] = TimeStep * Eq[j](time[i] + TimeStep, temp_x); }

					for (int j = 0; j < eqnum; ++j)
					{
						x[j, i + 1] = x[j, i] + (k1[j] + 2 * (k2[j] + k3[j]) + k4[j]) * 0.1666666666666667;
					}
				}
			}
		}

		/// <summary>
		/// 解を取得する
		/// </summary>
		/// <returns>解</returns>
		public double[,] GetSolution()
		{
			return x;
		}

		/// <summary>
		/// 解に対応した時刻列を取得する
		/// </summary>
		/// <returns>時刻列</returns>
		public double[] GetTime()
		{
			return time;
		}

		/// <summary>
		/// 指定された独立変数の解の最大値を取得する
		/// </summary>
		/// <param name="line"></param>
		/// <returns></returns>
		public double GetMaxValue(int line)
		{
			if ((line < 0) || (Eqnum <= line) || (x.GetLength(1) < 2))
			{
				throw new ArgumentException();
			}
			double max = x[line, 0];
			for (int i = 1; i < x.GetLength(1); i++)
			{
				if (x[line, i] > max)
				{
					max = x[line, i];
				}
			}
			return max;
		}

		/// <summary>
		/// 指定された独立変数の解の最小値を取得する
		/// </summary>
		/// <param name="line"></param>
		/// <returns></returns>
		public double GetMinValue(int line)
		{
			if ((line < 0) || (Eqnum <= line) || (x.GetLength(1) < 2))
			{
				throw new ArgumentException();
			}
			double min = x[line, 0];
			for (int i = 1; i < x.GetLength(1); i++)
			{
				if (x[line, i] < min)
				{
					min = x[line, i];
				}
			}
			return min;
		}

		/// <summary>
		/// 初期条件のセット
		/// </summary>
		/// <param name="ic">初期条件</param>
		private void setInitialCondition(double[] ic)
		{
			if (ic.Length != Eqnum)
			{
				throw new ArgumentException();
			}
			for (int i = 0; i < Eqnum; i++)
			{
				x[i, 0] = ic[i];
			}
		}

		private List<SystemEq> Eq;
		private double TimeStep { get; set; }		// 時間刻み[s]
		private int Eqnum { get; set; }				// 独立変数の数
		private double[] time;						// 時刻列
		private double[,] x;						// 解
	}


	namespace LinearAlgebra
	{
		public static class Vector
		{
			public static void Init(ref double[] vector)
			{

				for (int i = 0; i < vector.Length; i++)
				{
					vector[i] = 0;
				}
			}

			public static void Init(ref double[] vector, double c)
			{
				int i = 0;

				for (i = 0; i < vector.Length; i++)
				{
					vector[i] = c;
				}
			}


			public static int Max(ref int[] vector)
			{
				int i = 0;
				int temp_max = vector[0];

				for (i = 1; i < vector.GetLength(0); i++)
				{
					if (vector[i] > temp_max)
					{
						temp_max = vector[i];
					}
				}
				return temp_max;
			}

			public static int AbsMax(ref int[] vector)
			{
				int i = 0;
				int temp_absmax = vector[0];

				for (i = 1; i < vector.GetLength(0); i++)
				{
					if (Math.Abs(vector[i]) > temp_absmax)
					{
						temp_absmax = Math.Abs(vector[i]);
					}
				}
				return temp_absmax;
			}

			public static double AbsMax(ref double[] vector)
			{
				int i = 0;
				double temp_absmax = vector[0];

				for (i = 1; i < vector.GetLength(0); i++)
				{
					if (Math.Abs(vector[i]) > temp_absmax)
					{
						temp_absmax = Math.Abs(vector[i]);
					}
				}
				return temp_absmax;
			}

			public static int Min(ref int[] vector)
			{
				int i = 0;
				int temp_min = vector[0];

				for (i = 1; i < vector.GetLength(0); i++)
				{
					if (vector[i] < temp_min)
					{
						temp_min = vector[i];
					}
				}
				return temp_min;
			}
		}

		public static class Matrix
		{
			public static Boolean InitUnit(ref double[,] matrix)
			{
				int i = 0, j = 0;
				int row = matrix.GetLength(0);
				int col = matrix.GetLength(1);

				if (row == col)
				{
					for (i = 0; i < row; i++)
					{
						for (j = 0; j < col; j++)
						{
							if (i == j)
							{
								matrix[i, j] = 1;
							}
							else
							{
								matrix[i, j] = 0;
							}
						}
					}
					return true;
				}
				else
				{
					return false;
				}
			}

			public static void InitZero(ref double[,] matrix)
			{
				int i = 0, j = 0;
				int row = matrix.GetLength(0);
				int col = matrix.GetLength(1);

				for (i = 0; i < row; i++)
				{
					for (j = 0; j < col; j++)
					{
						matrix[i, j] = 0;
					}
				}
			}

			public static double[,] Sum(ref double[,] matrix1, ref double[,] matrix2)
			{
				int i = 0, j = 0;
				int row = matrix1.GetLength(0);
				int col = matrix1.GetLength(1);
				double[,] sum_rt = new double[row, col];

				for (i = 0; i < row; i++)
				{
					for (j = 0; j < col; j++)
					{
						sum_rt[i, j] = matrix1[i, j] + matrix2[i, j];
					}
				}

				return (double[,])sum_rt.Clone();
			}

			public static double[,] Mult(ref double[,] m1, ref double[,] m2)
			{
				int i = 0, j = 0, k = 0;
				int m1_l = 0, m1_r = 0, m2_l = 0, m2_r = 0;
				double[,] dummy_rt = new double[,] { { 0, 0 }, { 0, 0 } };

				m1_l = m1.GetLength(0);
				m1_r = m1.GetLength(1);
				m2_l = m2.GetLength(0);
				m2_r = m2.GetLength(1);

				if (m1_r != m2_l)
				{
					return (double[,])dummy_rt.Clone();	//とりあえず零を返す
				}
				else
				{
					double[,] comp = new double[m1_l, m2_r];

					for (i = 0; i < m1_l; i++)
					{
						for (j = 0; j < m2_r; j++)
						{
							comp[i, j] = 0;		//初期化しつつ・・・
							for (k = 0; k < m1_r; k++)
							{
								comp[i, j] += m1[i, k] * m2[k, j];
							}
						}
					}
					return (double[,])comp.Clone();
				}
			}

			public static double[,] ConstMult(ref double[,] matrix, double c)
			{
				int i = 0, j = 0;
				int row = matrix.GetLength(0);
				int col = matrix.GetLength(1);
				double[,] matrix_rt = new double[row, col];

				for (i = 0; i < row; i++)
				{
					for (j = 0; j < col; j++)
					{
						matrix_rt[i, j] = matrix[i, j] * c;
					}
				}
				return (double[,])matrix_rt.Clone();
			}

			public static double[,] LinearCombination(double c1, ref double[,] matrix1, double c2, ref double[,] matrix2)
			{
				int i = 0, j = 0;
				int row1 = matrix1.GetLength(0);
				int col1 = matrix1.GetLength(1);
				int row2 = matrix2.GetLength(0);
				int col2 = matrix2.GetLength(1);
				double[,] dummy_rt = new double[,] { { 0, 0 }, { 0, 0 } };

				if ((row1 == row2) && (col1 == col2))
				{
					double[,] matrix_rt = new double[row1, col1];

					for (i = 0; i < row1; i++)
					{
						for (j = 0; j < col1; j++)
						{
							matrix_rt[i, j] = c1 * matrix1[i, j] + c2 * matrix2[i, j];
						}
					}

					return (double[,])matrix_rt.Clone();
				}
				else
				{
					return (double[,])dummy_rt.Clone();
				}
			}

			public static double[,] Trans(ref double[,] m)
			{
				int i = 0, j = 0;
				int row = m.GetLength(0);
				int col = m.GetLength(1);
				double[,] trans_rt = new double[col, row];

				for (i = 0; i < col; i++)
				{
					for (j = 0; j < row; j++)
					{
						trans_rt[i, j] = m[j, i];
					}
				}
				return (double[,])trans_rt.Clone();
			}

			public static double Det2(ref double[,] mat)
			{
				return mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0];
			}

			public static double[,] Inv2(ref double[,] mat)
			{
				double det = 0;
				double[,] temp = new double[2, 2];

				det = Matrix.Det2(ref mat);
				temp[0, 0] = mat[1, 1] / det;
				temp[0, 1] = -mat[0, 1] / det;
				temp[1, 0] = -mat[1, 0] / det;
				temp[1, 1] = mat[0, 0] / det;

				return (double[,])temp.Clone();
			}

			public static double Det3(ref double[,] mat)
			{
				return (mat[0, 0] * (mat[1, 1] * mat[2, 2] - mat[1, 2] * mat[2, 1])
							+ mat[1, 0] * (mat[2, 1] * mat[0, 2] - mat[2, 2] * mat[0, 1])
							+ mat[2, 0] * (mat[0, 1] * mat[1, 2] - mat[0, 2] * mat[1, 1])
						);
			}

			public static double Determinant(ref double[,] mat)		//再帰で定義通り計算（超絶遅い）
			{
				int i = 0;
				double[,] submat;
				double cofactor = 0;
				double det = 0;

				for (i = 0; i < mat.GetLength(0); i++)
				{
					if (mat.GetLength(0) > 2)
					{
						if (mat[i, 0] == 0)
						{
							det += 0;
						}
						else
						{
							submat = new double[mat.GetLength(0) - 1, mat.GetLength(0) - 1];
							MinerMat(ref mat, ref submat, i, 0);
							cofactor = Determinant(ref submat);
							if (i % 2 == 0)
							{
								det += mat[i, 0] * cofactor;
							}
							else
							{
								det += mat[i, 0] * cofactor * (-1);
							}
						}
					}
					else
					{
						det = Det2(ref mat);
					}
				}
				return det;
			}

			private static void MinerMat(ref double[,] mat, ref double[,] submat, int ri, int rj)
			{
				int i = 0, j = 0;

				for (i = 0; i < mat.GetLength(0); i++)
				{
					for (j = 0; j < mat.GetLength(1); j++)
					{
						if ((i < ri) && (j < rj))
						{
							submat[i, j] = mat[i, j];
						}
						else if ((i == ri) || (j == rj))
						{
							//何もしない
						}
						else if ((i < ri) && (j > rj))
						{
							submat[i, j - 1] = mat[i, j];
						}
						else if ((i > ri) && (j < rj))
						{
							submat[i - 1, j] = mat[i, j];
						}
						else
						{
							submat[i - 1, j - 1] = mat[i, j];
						}
					}
				}
			}

			public static double Det(ref double[,] mat)		//http://thira.plavox.info/blog/2008/06/_c.html　からコピペ
			{
				int i = 0, j = 0, k = 0;
				double det = 1.0;
				double buf = 0;
				int n = mat.GetLength(0);  //配列の次数

				//三角行列を作成
				for (i = 0; i < n; i++)
				{
					for (j = 0; j < n; j++)
					{
						if (i < j)
						{
							buf = mat[j, i] / mat[i, i];
							for (k = 0; k < n; k++)
							{
								mat[j, k] -= mat[i, k] * buf;
							}
						}
					}
				}
				//対角部分の積
				for (i = 0; i < n; i++)
				{
					det *= mat[i, i];
				}

				return det;
			}

			private static Boolean CheckSymmetric(ref double[,] mat, Boolean disp)
			{
				int i = 0, j = 0;
				Boolean temp = true;

				if (disp == true)
					Console.Write("Symmetric matrix? ...... ");

				for (i = 0; i < mat.GetLength(0); i++)
				{
					for (j = i + 1; j < mat.GetLength(1); j++)
					{
						if (Math.Abs(mat[i, j] - mat[j, i]) > 1e-8)
						{
							temp = false;
						}
					}
				}

				if (temp == true)
				{
					if (disp == true)
						Console.WriteLine("OK!!");
					return true;
				}
				else
				{
					if (disp == true)
						Console.WriteLine("⊂ミ⊃＾ω＾ ）⊃　ｱｳｱｳ!!");
					return false;
				}
			}

			private static Boolean CheckSymmetric(ref SymmetricBandMatrix mat, Boolean disp)
			{
				int i = 0, j = 0;
				Boolean temp = true;

				if (disp == true)
					Console.Write("Symmetric matrix? ...... ");

				for (i = 0; i < mat.GetLength(0); i++)
				{
					for (j = i + 1; j < mat.GetLength(1); j++)
					{
						if (Math.Abs(mat[i, j] - mat[j, i]) > 1e-5)
						{
							temp = false;
						}
					}
				}

				if (temp == true)
				{
					if (disp == true)
						Console.WriteLine("OK!!");
					return true;
				}
				else
				{
					if (disp == true)
						Console.WriteLine("⊂ミ⊃＾ω＾ ）⊃　ｱｳｱｳ!!");
					return false;
				}
			}

			private static Boolean CheckSingular(ref double[,] mat, Boolean disp)
			{
				//	int i = 0, j = 0;
				Boolean temp = true;

				return temp;
			}

			private static void Swap<Type>(ref Type lhs, ref Type rhs)
			{
				Type temp;
				temp = lhs;
				lhs = rhs;
				rhs = temp;
			}

			public static Boolean SwapRows(ref double[,] matrix, int r1, int r2)
			{
				int i = 0;
				int row = matrix.GetLength(0);
				int col = matrix.GetLength(1);

				if ((r1 <= row) && (r2 <= row))
				{
					for (i = 0; i < col; i++)
					{
						Swap<double>(ref matrix[r1, i], ref matrix[r2, i]);
					}
					return true;
				}
				else
				{
					return false;
				}

			}

			public static Boolean SwapRows<Type>(ref Type[,] matrix, int r1, int r2)
			{
				int i = 0;
				int row = matrix.GetLength(0);
				int col = matrix.GetLength(1);

				if ((r1 <= row) && (r2 <= row))
				{
					for (i = 0; i < col; i++)
					{
						Swap<Type>(ref matrix[r1, i], ref matrix[r2, i]);
					}
					return true;
				}
				else
				{
					return false;
				}

			}

			public static Boolean SwapCols<Type>(ref Type[,] matrix, int c1, int c2)
			{
				int i = 0;
				int row = matrix.GetLength(0);
				int col = matrix.GetLength(1);

				if ((c1 <= col) && (c2 <= col))
				{
					for (i = 0; i < row; i++)
					{
						Swap<Type>(ref matrix[i, c1], ref matrix[i, c2]);
					}
					return true;
				}
				else
				{
					return false;
				}

			}

			public static double LUdecomp(double[,] m, double[,] mlu, double[] pivot)	//行列 m を LU 分解し、結果を mlu と pivot に入れ、また m の行列式を返す
			{
				int i = 0, j = 0, k = 0, ipiv = 0;
				double det = 0, max = 0, tmp = 0, weight = 0;

				mlu = (double[,])m.Clone();	/* m を mlu にコピー */

				for (i = 0; i < pivot.Length; i++)	/* pivot の初期化 */
				{
					pivot[i] = i;
				}

				det = 1; /* 行列式初期値 */
				for (i = 0; i < pivot.Length; i++)	/* pivot 選択 */
				{
					max = 0;
					ipiv = (int)pivot[i];
					for (j = i; j < pivot.Length; j++)	/* 行の重み weight を調べる */
					{
						weight = 0;
						for (k = j; k < pivot.Length; k++)	/* j行の絶対値最大の要素を探す */
						{
							tmp = Math.Abs(m[j, k]);
							if (tmp > weight)
							{
								weight = tmp;
							}
						}
						if (weight == 0)	/* 行列 m は正則でない */
						{
							Console.WriteLine("This matrix is singular![in LUdecomp()]");
							det = 0;
							return det;
						}
						tmp = Math.Abs(mlu[j, i]) / weight;	/* 重みづけした絶対値を max と比較 */
						if (tmp > max)
						{
							max = tmp;
							ipiv = j;
						}
					}
					if (ipiv != i)	/* 選択した pivot 行と入れ換え、それを記憶する */
					{
						SwapRows(ref mlu, i, ipiv);
						Swap(ref pivot[i], ref pivot[ipiv]);
						det = -det;    /* 行を交換すれば行列式の符号が変わる */
					}
					det *= mlu[i, i];	/* 行列式は U の対角成分の積 */
					if (det == 0)	/* 行列 m は正則でない */
					{
						Console.WriteLine("This matrix is singular![in LUdecomp()]");
						det = 0;
						return det;
					}
					for (j = i + 1; j < pivot.Length; j++)	/* LU分解 */
					{
						mlu[j, i] /= mlu[i, i];
						for (k = i + 1; k < pivot.Length; k++)
						{
							mlu[j, k] -= mlu[j, i] * mlu[i, k];
						}
					}
				}
				return det;
			}

			public static void InvMatLU(double[,] mlu, double[,] minv, double[] pivot)	//LU分解済みの行列mlu[][]の逆行列minv[]を求める
			{
				int i = 0, j = 0, k = 0;
				int kpiv = 0;

				for (i = 0; i < pivot.Length; i++)	//minv を pivot に従って単位行列の列を入れ換えたものにしておく
				{
					for (j = 0; j < pivot.Length; j++)
					{
						if ((int)pivot[i] == j)
						{
							minv[i, j] = 1;
						}
						else
						{
							minv[i, j] = 0;
						}
					}
				}

				for (k = 0; k < pivot.Length; k++)	//逆行列を求める
				{
					kpiv = (int)pivot[k];
					for (i = 0; i < pivot.Length; i++)
					{
						for (j = i + 1; j < pivot.Length; j++)
						{
							minv[j, kpiv] -= mlu[j, i] * minv[i, kpiv];
						}
					}
					for (i = pivot.Length - 1; i >= 0; i--)
					{
						if (mlu[i, i] == 0)	//行列 m は正則でない
						{
							Console.WriteLine("This matrix is singular![in InvMatLU()]");
							String[] args = Environment.GetCommandLineArgs();	// args[0] is the program name and, args[1] is the first argument.Test for a command-line argument.
							try
							{
								int exitCode = int.Parse(args[1]);
								Console.WriteLine("it exits with code: 0x{0:X8}.", exitCode);
								Environment.Exit(exitCode);
							}
							catch { }
							finally
							{
								Environment.Exit(0);
							}
						}

						for (j = i + 1; j < pivot.Length; j++)
						{
							minv[i, kpiv] -= mlu[i, j] * minv[j, kpiv];
						}
						minv[i, kpiv] /= mlu[i, i];
					}
				}
			}

			public static void LDLTCholeskyDecomp(ref double[,] mat, ref UnitLowerTriMatrix l, ref double[] d)
			{
				int i = 0, j = 0, k = 0;
				double temp = 0;

				for (j = 0; j < mat.GetLength(1); j++)		//i:行，j:列
				{
					for (i = j; i < mat.GetLength(0); i++)
					{
						if (i == j)		//対角項
						{
							d[j] = mat[j, j];
							for (k = 0; k < j; k++)
							{
								d[j] -= d[k] * l[j, k] * l[j, k];
							}
						}
						else		//非対角項
						{
							temp = mat[i, j];
							for (k = 0; k < j; k++)
							{
								temp -= d[k] * l[i, k] * l[j, k];
							}
							l[i, j] = temp / d[j];
						}
					}
				}
			}

			public static void LDLTCholeskyDecomp(ref SymmetricBandMatrix mat, ref UnitLowerTriMatrix l, ref double[] d)
			{
				int i = 0, j = 0, k = 0;
				double temp = 0;

				for (j = 0; j < mat.GetLength(1); j++)		//i:行，j:列
				{
					for (i = j; i < mat.GetLength(0); i++)
					{
						if (i == j)		//対角項
						{
							d[j] = mat[j, j];
							for (k = 0; k < j; k++)
							{
								d[j] -= d[k] * l[j, k] * l[j, k];
							}
						}
						else		//非対角項
						{
							temp = mat[i, j];
							for (k = 0; k < j; k++)
							{
								temp -= d[k] * l[i, k] * l[j, k];
							}
							l[i, j] = temp / d[j];
						}
					}
				}
			}

			public static void LDLTCholeskyDecomp(ref SymmetricMatrix mat, out UnitLowerTriMatrix l, out double[] d)
			{
				int i = 0, j = 0, k = 0;
				double temp = 0;
				l = new UnitLowerTriMatrix(mat.GetLength(0));
				d = new double[mat.GetLength(0)];

				for (j = 0; j < mat.GetLength(1); j++)		//i:行，j:列
				{
					for (i = j; i < mat.GetLength(0); i++)
					{
						if (i == j)		//対角項
						{
							d[j] = mat[j, j];
							for (k = 0; k < j; k++)
							{
								d[j] -= d[k] * l[j, k] * l[j, k];
							}
						}
						else		//非対角項
						{
							temp = mat[i, j];
							for (k = 0; k < j; k++)
							{
								temp -= d[k] * l[i, k] * l[j, k];
							}
							l[i, j] = temp / d[j];
						}
					}
				}
			}

			public static void LDLTSolve(ref SymmetricMatrix a, ref double[] y, out double[] x)
			{
				int i = 0, j = 0;
				int n = y.Length;
				double temp = 0;

				x = new double[n];
				UnitLowerTriMatrix l = new UnitLowerTriMatrix(n);	//LDLT分解のL行列（下三角項のみ）
				double[] d = new double[n];							//LDLT分解のD行列（対角項のみ）

				LDLTCholeskyDecomp(ref a, out l, out d);

				for (i = 0; i < n; i++)
				{
					temp = 0;
					for (j = 0; j < i; j++)
					{
						temp += l[i, j] * x[j];
					}
					x[i] = (y[i] - temp) / l[i, i];
				}
				for (i = 0; i < n; i++)
				{
					x[i] /= d[i];
				}

				for (i = n - 1; i >= 0; i--)
				{
					temp = 0;
					for (j = n - 1; j > i; j--)
					{
						temp += l[j, i] * x[j];
					}
					x[i] = (y[i] - temp) / l[i, i];
				}
			}
		}

		public class BandMatrix		//帯行列を通常の行列のインデックスを用いて取り扱うクラス
		{
			public BandMatrix(int dim, int bandwidth)
			{
				Bandwidth = bandwidth;
				mat = new double[dim, Bandwidth];
				size.Add(dim);
				size.Add(Bandwidth);
			}

			private double[,] mat;
			private ArrayList size = new ArrayList();
			private int bandwidth;
			public int Bandwidth { get { return bandwidth; } private set { this.bandwidth = value; } }

			public double this[int i, int j]
			{
				set
				{
					if (Math.Abs(i - j) >= (Bandwidth + 1) / 2)
					{
						//何もしない
					}
					else if (i < (Bandwidth + 1) / 2)
					{
						this.mat[i, j] = value;
					}
					else
					{
						this.mat[i, j - (i - (Bandwidth + 1) / 2 + 1)] = value;
					}
				}
				get
				{
					if (Math.Abs(i - j) >= (Bandwidth + 1) / 2)
					{
						return 0;
					}
					else if (i < (Bandwidth + 1) / 2)
					{
						return this.mat[i, j];
					}
					else
					{
						return this.mat[i, j - (i - (Bandwidth + 1) / 2 + 1)];
					}
				}
			}

			public int GetLength(int a)
			{
				return (int)size[a];
			}
		}

		public class SymmetricMatrix
		{
			public SymmetricMatrix(int dim)
			{
				int i = 0;
				Dim = dim;
				mat = new double[Dim][];

				for (i = 0; i < dim; i++)
				{
					mat[i] = new double[i + 1];
				}

				size.Add(Dim);
				size.Add(Dim);
			}

			private double[][] mat;
			private int dimension;
			public int Dim { get { return dimension; } private set { this.dimension = value; } }
			private ArrayList size = new ArrayList(2);

			public double this[int i, int j]
			{
				set
				{
					if (i >= j)
					{
						mat[i][j] = value;
					}
					else
					{
						mat[j][i] = value;
					}
				}
				get
				{
					if (i >= j)
					{
						return mat[i][j];
					}
					else
					{
						return mat[j][i];
					}
				}
			}

			public int GetLength(int a)
			{
				return (int)size[a];
			}
		}

		public class SymmetricBandMatrix	//スカイライン法により格納した対称帯行列を通常の行列のインデックスを用いて取り扱うクラス．
		{
			public SymmetricBandMatrix(int dim, int bandwidth)
			{
				int i = 0;
				Bandwidth = bandwidth;
				Dim = dim;
				mat3 = new double[Dim][];
				for (i = 0; i < Dim; i++)
				{
					if (i < (Dim - Bandwidth))
					{
						mat3[i] = new double[Bandwidth];
					}
					else
					{
						mat3[i] = new double[Dim - i];
					}
				}
				size.Add(Dim);
				size.Add(Bandwidth);
			}

			private double[][] mat3;
			private int bandwidth;
			public int Bandwidth { get { return bandwidth; } private set { this.bandwidth = value; } }
			private int dimension;
			public int Dim { get { return dimension; } private set { this.dimension = value; } }
			public int[] band;		//行ごとのバンド幅
			private ArrayList size = new ArrayList(2);

			public double this[int i, int j]
			{
				set
				{
					if ((i < 0) || (i >= Dim) || (j < 0) || (j >= Dim))
					{
						Environment.Exit(0);
					}
					else if (((i - j) >= Bandwidth) || ((j - i) >= Bandwidth))
					{
						//
					}
					else
					{
						if (i <= j)
						{
							this.mat3[i][j - i] = value;
						}
						else
						{
							this.mat3[j][i - j] = value;
						}
					}

				}
				get
				{
					if ((i < 0) || (i >= Dim) || (j < 0) || (j >= Dim))
					{
						Environment.Exit(0);
						return 0;
					}
					else if (((i - j) >= Bandwidth) || ((j - i) >= Bandwidth))
					{
						return 0;
					}
					else
					{
						if (i <= j)
						{
							return this.mat3[i][j - i];
						}
						else
						{
							return this.mat3[j][i - j];
						}
					}
				}
			}

			public int GetLength(int a)
			{
				return (int)size[a];
			}
		}

		public class LowerTriMatrix		//下三角行列クラス
		{
			public LowerTriMatrix(int dim)
			{
				int i = 0;

				mat = new double[dim][];
				for (i = 0; i < dim; i++)
				{
					mat[i] = new double[i + 1];
				}
			}

			private double[][] mat;

			public double this[int i, int j]
			{
				set		//下三角領域以外への書き込みはアウト
				{
					mat[i][j] = value;
				}
				get
				{
					if (i < j)
					{
						return 0;
					}
					else
					{
						return mat[i][j];
					}
				}
			}
		}

		public class UnitLowerTriMatrix		//単位下三角行列クラス
		{
			public UnitLowerTriMatrix(int dim)
			{
				int i = 0;

				mat = new double[dim - 1][];
				for (i = 0; i < (dim - 1); i++)
				{
					mat[i] = new double[i + 1];
				}
			}

			private double[][] mat;

			public double this[int i, int j]
			{
				set		//下三角領域以外への書き込みはアウト
				{
					mat[i - 1][j] = value;
				}
				get
				{
					if (i < j)
					{
						return 0;
					}
					else if (i == j)
					{
						return 1;
					}
					else
					{
						return mat[i - 1][j];
					}
				}
			}
		}
	}

	namespace FastFourierTransform
	{
		public class FastFourierTransform
		{
			/// <summary>
			/// 1次元離散フーリエ変換．定義式通りの実装で低速．
			/// </summary>
			public static void DFT(ref Complex[] x, out Complex[] f)
			{
				int i = 0, j = 0;
				int n = 0;
				Complex temp;

				n = x.Length;
				f = new Complex[n];

				for (i = 0; i < n; i++)
				{
					temp = new Complex(0.0, 0.0);
					for (j = 0; j < n; j++)
					{
						temp += Complex.Multiply(x[j], new Complex(Math.Cos(2 * Math.PI * i * j / n), -Math.Sin(2 * Math.PI * i * j / n)));
					}
					f[i] = temp;
				}
			}

			/// <summary>
			/// 1次元離散フーリエ変換．高速ではない．実数入力に限定．
			/// </summary>
			public static void DFT(ref double[] x, out Complex[] f)
			{
				int i = 0, j = 0;
				int n = 0;
				Complex temp;

				if (x.Length % 2 == 0)
				{
					n = x.Length / 2;
				}
				else
				{
					n = (x.Length + 1) / 2;
				}
				f = new Complex[n];

				for (i = 0; i < n; i++)
				{
					temp = new Complex(0.0, 0.0);
					for (j = 0; j < n; j++)
					{
						temp += new Complex(x[j] * Math.Cos(2 * Math.PI * i * j / n), -x[j] * Math.Sin(2 * Math.PI * i * j / n));
					}
					f[i] = temp;
				}
			}



			/// <summary>
			/// 1次元FFT 
			/// </summary>
			public static void FFT(double[] inputRe, double[] inputIm, out double[] outputRe, out double[] outputIm, int bitSize)
			{
				int dataSize = 1 << bitSize;
				int[] reverseBitArray = BitScrollArray(dataSize);

				outputRe = new double[dataSize];
				outputIm = new double[dataSize];

				// バタフライ演算のための置き換え
				for (int i = 0; i < dataSize; i++)
				{
					outputRe[i] = inputRe[reverseBitArray[i]];
					outputIm[i] = inputIm[reverseBitArray[i]];
				}

				// バタフライ演算
				for (int stage = 1; stage <= bitSize; stage++)
				{
					int butterflyDistance = 1 << stage;
					int numType = butterflyDistance >> 1;
					int butterflySize = butterflyDistance >> 1;

					double wRe = 1.0;
					double wIm = 0.0;
					double uRe = System.Math.Cos(System.Math.PI / butterflySize);
					double uIm = -System.Math.Sin(System.Math.PI / butterflySize);

					for (int type = 0; type < numType; type++)
					{
						for (int j = type; j < dataSize; j += butterflyDistance)
						{
							int jp = j + butterflySize;
							double tempRe = outputRe[jp] * wRe - outputIm[jp] * wIm;
							double tempIm = outputRe[jp] * wIm + outputIm[jp] * wRe;
							outputRe[jp] = outputRe[j] - tempRe;
							outputIm[jp] = outputIm[j] - tempIm;
							outputRe[j] += tempRe;
							outputIm[j] += tempIm;
						}
						double tempWRe = wRe * uRe - wIm * uIm;
						double tempWIm = wRe * uIm + wIm * uRe;
						wRe = tempWRe;
						wIm = tempWIm;
					}
				}
			}

			// ビットを左右反転した配列を返す
			private static int[] BitScrollArray(int arraySize)
			{
				int[] reBitArray = new int[arraySize];
				int arraySizeHarf = arraySize >> 1;

				reBitArray[0] = 0;
				for (int i = 1; i < arraySize; i <<= 1)
				{
					for (int j = 0; j < i; j++)
						reBitArray[j + i] = reBitArray[j] + arraySizeHarf;
					arraySizeHarf >>= 1;
				}
				return reBitArray;
			}
		}
	}

}
