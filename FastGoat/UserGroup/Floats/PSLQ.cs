using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Floats;

/*
The two-level multipair PSLQ algorithm
David H. Bailey
May 2, 2024

Fortran source code
mpfun20-mpfr-v32.tar.gz
tpslqm2.f90
 */
public static class PSLQ
{
    /// <summary>
    /// Step Iteration of one level multipair PSLQ
    /// BigReal arbitrary precision
    /// </summary>
    public static void IterOneLevelMultipair(KMatrix<BigReal> H, KMatrix<BigReal> A, KMatrix<BigReal> B,
        KMatrix<BigReal> T, KMatrix<BigReal> y, (int i, BigReal yi)[] gamma_pow, bool imq = false)
    {
        // 1. Sort the entries of the (n − 1)-long vector {γ^i |Hii |} in decreasing order,
        // producing the sort indices.
        // 2. Beginning at the sort index m1 corresponding to the largest γ^i |Hii |, select
        // pairs of indices (mi , mi + 1), where mi is the sort index. If at any step
        // either mi or mi + 1 has already been selected, pass to the next index in
        // the list. Continue until either βn pairs have been selected, or the list
        // is exhausted. Let p denote the number of pairs actually selected in this
        // manner.
        // 3. For i := 1 to p, exchange the entries of y indexed mi and mi + 1, and the
        // corresponding rows of A, B and H; endfor.
        var m = gamma_pow.Select(e => (e.i, (e.yi * H[e.i, e.i].Absolute)))
            .OrderByDescending(e => e.Item2).ToArray();
        var list = new List<(int, int)>();
        var selected = new HashSet<int>();
        var n = gamma_pow.Length + 1;
        for (int i = 0; i < n - 1; i++)
        {
            var (b, c) = (m[i].i, m[i].i + 1);
            if (selected.Contains(b) || selected.Contains(c))
                continue;

            list.Add((b, c));
            selected.UnionWith([b, c]);
            if (imq || list.Count + 1 > 0.4 * n)
                break;
        }

        foreach (var (mi, _) in list)
        {
            // Console.WriteLine($"Swap {(mi, mi + 1)}");
            (y.Coefs[mi, 0], y.Coefs[mi + 1, 0]) = (y[mi + 1, 0], y[mi, 0]);
            for (int k = 0; k < n; k++)
            {
                (A.Coefs[mi, k], A.Coefs[mi + 1, k]) = (A.Coefs[mi + 1, k], A.Coefs[mi, k]);
                (B.Coefs[mi, k], B.Coefs[mi + 1, k]) = (B.Coefs[mi + 1, k], B.Coefs[mi, k]);
                if (k < n - 1)
                    (H.Coefs[mi, k], H.Coefs[mi + 1, k]) = (H.Coefs[mi + 1, k], H.Coefs[mi, k]);
            }
        }

        // 4. Remove corners on H diagonal: For j := 1 to p: if mj ≤ n − 2 then set
        // t0 := Sqrt(H2mj_mj + H2mj_mj1), t1 := Hmj_mj /t0 and t2 := Hmj_mj1 /t0 ; for 
        // i := mj to n: set t3 := Hi_mj ; t4 := Hi_mj1 ; Hi_mj := t1t3 + t2t4 ; and
        // Hi_mj1 := −t2t3 + t1t4 ; endfor; endif; endfor.
        for (int j = 0; j < list.Count; j++)
        {
            var mj = list[j].Item1;
            if (mj < n - 2)
            {
                var t0 = BigReal.Sqrt(H[mj, mj].Pow(2) + H[mj, mj + 1].Pow(2));
                var t1 = H[mj, mj] / t0;
                var t2 = H[mj, mj + 1] / t0;
                for (int i = mj; i < n; i++)
                {
                    var t3 = H[i, mj];
                    var t4 = H[i, mj + 1];
                    H.Coefs[i, mj] = t1 * t3 + t2 * t4;
                    H.Coefs[i, mj + 1] = -t2 * t3 + t1 * t4;
                }
            }
        }

        for (int i = 0; i < n; i++)
        for (int j = 0; j < n - 1; j++)
            T.Coefs[i, j] = T.KZero;

        // 5. Reduce H: For i := 2 to n: for j := 1 to n − i + 1: set l := i + j − 1; for
        // k := j +1 to l −1: set Hlj := Hlj −Tlk Hkj ; endfor; set Tlj := nint(Hlj /Hjj )
        // and Hlj := Hlj − Tlj Hjj ; endfor; endfor. [Note that the n x (n − 1) integer
        // array T is set before it is used.]
        for (int i = 1; i < n; i++)
        {
            for (int j = 0; j < n - i; j++)
            {
                var l = i + j;
                for (int k = j + 1; k < l; k++)
                    H.Coefs[l, j] -= T[l, k] * H[k, j];

                T.Coefs[l, j] = (H[l, j] / H[j, j]).RoundEven;
                H.Coefs[l, j] -= T[l, j] * H[j, j];
            }
        }

        // 6. Update y: For j := 1 to n − 1: for i := j + 1 to n: set yj := yj + Tij yi ;
        // endfor; endfor.
        // 7. Update A and B: For k := 1 to n: for j := 1 to n − 1: for i := j + 1 to n:
        // set Aik := Aik − Tij Ajk and Bjk := Bjk + Tij Bik ; endfor; endfor; endfor.
        for (int j = 0; j < n - 1; j++)
        {
            for (int i = j + 1; i < n; i++)
                y.Coefs[j, 0] += T[i, j] * y[i, 0];
        }

        for (int k = 0; k < n; k++)
        {
            for (int j = 0; j < n - 1; j++)
            {
                for (int i = j + 1; i < n; i++)
                {
                    A.Coefs[i, k] -= T[i, j] * A[j, k];
                    B.Coefs[j, k] += T[i, j] * B[i, k];
                }
            }
        }
    }

    static void SaveRestore<T>(KMatrix<T> H0, KMatrix<T> A0, KMatrix<T> B0, KMatrix<T> T0, KMatrix<T> y0,
        KMatrix<T> H1, KMatrix<T> A1, KMatrix<T> B1, KMatrix<T> T1, KMatrix<T> y1)
        where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
    {
        var N = H0.M;
        Array.Copy(H0.Coefs, H1.Coefs, N * (N - 1));
        Array.Copy(A0.Coefs, A1.Coefs, N * N);
        Array.Copy(B0.Coefs, B1.Coefs, N * N);
        Array.Copy(T0.Coefs, T1.Coefs, N * (N - 1));
        Array.Copy(y0.Coefs, y1.Coefs, N);
    }

    /// <summary>
    /// One-level multipair PSLQ algorithm
    /// </summary>
    /// <param name="x">Row Vector</param>
    /// <param name="gamma">PSLQ gamma</param>
    /// <returns>Integer Relation coefficients</returns>
    /// <exception cref="Exception">Precision exhausted</exception>
    public static Rational[] OnelevelMultipair(KMatrix<BigReal> x, BigReal gamma)
    {
        // Initialize:
        // 1. For j := 1 to n: for i := 1 to n: if i = j then set Aij := 1 and Bij := 1
        // else set Aij := 0 and Bij := 0; endfor; endfor.
        var n = x.N;
        var z = BigReal.BrZero(x.KOne.O);
        var A = new KMatrix<BigReal>(z, n, n).One;
        var B = new KMatrix<BigReal>(z, n, n).One;
        var T = new KMatrix<BigReal>(z, n, n - 1);

        // 2. For k := 1 to n: set sk :=Sqrt(Sum[j=k to n] xj^2) ; endfor; set t = 1/s1 ; for k := 1 to
        // n: set yk := txk ; sk := tsk ; endfor.
        var s = n.Range().Select(j => (n - j).Range(j).Aggregate(x.KZero, (acc, k) => acc + x[0, k].Pow(2)))
            .Select(e => BigReal.Sqrt(e)).ToArray();
        var t = s[0].Inv();
        var y = x.Select(xi => xi * t).ToKMatrix(n);
        s = s.Select(si => si * t).ToArray();

        // 3. Initial H: For j := 1 to n − 1: for i := 1 to j − 1: set Hij := 0; endfor;
        // set Hjj := sj+1 /sj ; for i := j + 1 to n: set Hij := −yi yj /(sj sj+1 ); endfor;
        // endfor.
        var H = new KMatrix<BigReal>(z, n, n - 1);
        for (int j = 0; j < n - 1; j++)
        {
            for (int i = 0; i < j - 1; i++)
                H.Coefs[i, j] = z;

            H.Coefs[j, j] = s[j + 1] / s[j];
            for (int i = j + 1; i < n; i++)
                H.Coefs[i, j] = (-y[i, 0] * y[j, 0]) / (s[j] * s[j + 1]);
        }

        var step = 0;
        var O1 = z.O;
        var O2 = gamma.O;
        var gamma_pow = (n - 1).Range().Select(i => (i, yi: gamma.Pow(i + 1))).ToArray();

        // Iteration: Repeat the following steps until precision has been exhausted or a
        // relation has been detected.
        while (true)
        {
            ++step;
            IterOneLevelMultipair(H, A, B, T, y, gamma_pow);

            // 8. Norm bound: Compute M := 1/ maxj |Hjj |. Then there can exist no
            // relation vector whose Euclidean norm is less than M .
            // 9. Termination test: If the largest entry of A exceeds the level of numeric
            // precision used, then precision is exhausted. If the smallest entry of the y
            // vector is less than the detection epsilon, a relation has been detected and
            // is given in the corresponding row of B.
            var M = (n - 1).Range().Max(j => H[j, j]).Inv();
            if (A.Any(e => e.V >= O1 + 2))
                throw new("Precision is exhausted");

            if (y.Any(e => e.ToBigReal(O1 - 4).IsZero()))
            {
                var ym = y.Select((yi, i) => (e: yi, i)).OrderBy(c => c.e.Absolute).First();
                if (Logger.Level != LogLevel.Off)
                {
                    Console.WriteLine($"Possible Solution step:{step}");
                    Console.WriteLine(B.GetRow(ym.i).Select(c => c.RoundEven.ToRational).ToKMatrix());
                    Console.WriteLine();
                }

                return B.GetRow(ym.i).Select(c => c.RoundEven.ToRational).ToArray();
            }
        }

        throw new();
    }

    public static (KMatrix<BigReal> Q, KMatrix<BigReal> R) QR(KMatrix<BigReal> A)
    {
        var O = A.KOne.O;
        var (m, n) = A.Dim;
        var Qis = new List<KMatrix<BigReal>>();
        var Ai = new KMatrix<BigReal>(A.Coefs);
        var min = int.Min(m, n);
        for (int i = 0; i < min - 1; i++)
        {
            var x = Ai.GetCol(0);
            var sign = Ai[0, 0].IsZero() ? 1 : Ai[0, 0].Sign;
            var a = sign * BigReal.Sqrt(Ring.SquareNorm2(x));
            var u = x.Select((xk, k) => k == 0 ? xk + a : xk).ToKMatrix(x.M);

            var Ii = new KMatrix<BigReal>(Ring.Diagonal(BigReal.BrOne(O), x.M));
            var Qi = Ii - 2 * u * u.T / Ring.SquareNorm2(u);

            var Qik = i == 0
                ? Qi
                : KMatrix<BigReal>.MergeBlocks(new KMatrix<BigReal>(Ring.Diagonal(BigReal.BrOne(O), m - x.M)), Qi);
            Qis.Add(Qik);
            Ai = (Qi * Ai).Extract(1, Ai.M - 1, 1, Ai.N - 1);
        }

        var Qt = Qis.Aggregate(Qis[0].One, (acc, q) => q * acc);
        return (Qt.T, Qt * A);
    }

    public static (KMatrix<Dble> Q, KMatrix<Dble> R) QR(KMatrix<Dble> A)
    {
        var (m, n) = A.Dim;
        var Qis = new List<KMatrix<Dble>>();
        var Ai = new KMatrix<Dble>(A.Coefs);
        for (int i = 0; i < int.Min(m, n); i++)
        {
            if (i == n - 1)
                break;

            var x = Ai.GetCol(0);
            var sign = Ai[0, 0].IsZero() ? 1 : Ai[0, 0].Sign;
            var a = sign * Dble.Sqrt(Ring.SquareNorm2(x));
            var u = x.Select((xk, k) => k == 0 ? xk + a : xk).ToKMatrix(x.M);

            var Ii = new KMatrix<Dble>(Ring.Diagonal(Dble.DbleOne(), x.M));
            var Qi = Ii - 2 * u * u.T / Ring.SquareNorm2(u);

            var Qik = i == 0
                ? Qi
                : KMatrix<Dble>.MergeBlocks(new KMatrix<Dble>(Ring.Diagonal(Dble.DbleOne(), m - x.M)), Qi);
            Qis.Add(Qik);
            Ai = (Qi * Ai).Extract(1, Ai.M - 1, 1, Ai.N - 1);
        }

        var Qt = Qis.Aggregate(Qis[0].One, (acc, q) => q * acc);
        return (Qt.T, Qt * A);
    }

    public static KMatrix<BigReal> LQpslq(KMatrix<BigReal> A)
    {
        var H = new KMatrix<BigReal>(A.Coefs);
        var (n, m) = A.Dim;
        var min = int.Min(m, n);
        if (n < m)
            throw new("LQ rows >= cols");

        for (int l = 0; l < min; l++)
        {
            if (l == m - 1)
                continue;

            var N2 = (m - l).Range().Select(i => H[l, l + i]).Aggregate(A.KZero, (acc, a) => acc + a * a);
            var N = BigReal.Sqrt(N2);
            if (N.IsZero())
                continue;

            if (!H[l, l].IsZero())
                N *= H[l, l].Sign;

            var Ni = N.Inv();
            for (int i = 0; i < m - l; ++i)
                H.Coefs[l, l + i] *= Ni;

            H.Coefs[l, l] += 1;
            for (int j = l + 1; j < n; ++j)
            {
                var t = -((m - l).Range().Select(i => H[l, l + i] * H[j, l + i])
                            .Aggregate(A.KZero, (acc, a) => acc + a))
                        / H[l, l];

                for (int i = 0; i < m - l; i++)
                    H.Coefs[j, l + i] += t * H[l, l + i];
            }

            H.Coefs[l, l] = -N;
        }

        for (int j = 0; j < m; j++)
        for (int i = 0; i < j; i++)
            H.Coefs[i, j] = H.KZero;

        return new(H.Coefs);
    }

    public static KMatrix<Dble> LQpslq(KMatrix<Dble> A)
    {
        var H = new KMatrix<Dble>(A.Coefs);
        var (n, m) = A.Dim;
        if (n < m)
            throw new("LQ rows >= cols");

        for (int l = 0; l < int.Min(m, n); l++)
        {
            if (l == m - 1)
                continue;

            var N2 = (m - l).Range().Select(i => H[l, l + i]).Aggregate(A.KZero, (acc, a) => acc + a * a);
            var N = Dble.Sqrt(N2);
            if (N.IsZero())
                continue;

            if (!H[l, l].IsZero())
                N *= H[l, l].Sign;

            var Ni = N.Inv();
            for (int i = 0; i < m - l; ++i)
                H.Coefs[l, l + i] *= Ni;

            H.Coefs[l, l] += 1;
            for (int j = l + 1; j < n; ++j)
            {
                var t = -((m - l).Range().Select(i => H[l, l + i] * H[j, l + i])
                            .Aggregate(A.KZero, (acc, a) => acc + a))
                        / H[l, l];

                for (int i = 0; i < m - l; i++)
                    H.Coefs[j, l + i] += t * H[l, l + i];
            }

            H.Coefs[l, l] = -N;
        }

        for (int j = 0; j < m; j++)
        for (int i = 0; i < j; i++)
            H.Coefs[i, j] = H.KZero;

        return new(H.Coefs);
    }

    // PSLQM2 (two-level multipair PSLQ):
    // Input arguments and parameters:
    // N            Int     Length of input vector X and output relation vector R.
    // NDP          Int     Precision level in digits.
    // NDR          Int     log10 of the min acceptable dynamic range of Y vector at detection;
    //                      default = 30. A smaller range is deemed unreliable.
    // NRB          Int     log10 of max size (Euclidean norm) of acceptable relation;
    //                      default = 200.
    // NEP          Int     log10 of full precision epsilon; default = 30 − NDP ; for large
    //                      problems replace 30 by 50 or 60.
    // X            MPR     N-long input multiprecision vector.
    // IQ           Int     Output flag: 0 (unsuccessful) or 1 (successful).
    // R            MPR     N-long output integer relation vector, if successful, else zeroes.
    // NSQ          Int     Second dimension of DSYQ and SYQ; default = 8.
    // DEPS         DP      Double precision epsilon; default = 10^−14 .
    // DREP         DP      Dynamic range epsilon; default = 10^−10 .
    // Integer variables:
    // IT, ITS, IMQ, IZD, IZM
    // Double precision arrays and variables (with dimensions):
    // DA(N, N), DB(N, N), DH(N, N), DYSQ(N, NSQ), DY(N)
    // Multiprecision real arrays and variables (with dimensions):
    // B(N, N), H(N, N), SYQ(N, NSQ), Y(N), EPS, T
    public static KMatrix<BigReal> TwoLevelMultipair(KMatrix<BigReal> X, BigReal gamma)
    {
        // 1. Initialize MPR arrays:
        // a. Set EPS = 10^NEP .
        // b. Set B to an N * N identity matrix.
        // c. Compute initial H matrix as given in multipair PSLQ algorithm above.
        // d. Set SYQ array to zeroes.
        // e. Set IZD, IZM, IMQ, IT and ITS to zero.

        Dble.EpsDouble = 1e-14;
        var N = X.N;
        var o = X.KOne;
        int NDP = gamma.O, NDR = NDP / 20, NRB = NDP / 2, O2 = 16;
        int NEP = -(4 * NDP / 5), NSQ = 8;
        var DEPS = BigReal.FromDouble(1e-14, O2);
        var DREP = BigReal.FromDouble(1e-10, O2);

        var EPS = BigReal.FromDouble(double.Pow(10, NEP), NDP);
        var twoPow72 = BigReal.FromBigInteger(BigInteger.Pow(2, 72), NDP);
        var A = new KMatrix<BigReal>(Ring.Diagonal(o, N));
        var B = new KMatrix<BigReal>(Ring.Diagonal(o, N));
        var HT = new KMatrix<BigReal>(o, N, N - 1);

        var s = new BigReal[N];
        var T = o.Zero;
        for (int i = N - 1; i >= 0; i--)
        {
            T += X[0, i].Pow(2);
            s[i] = BigReal.Sqrt(T);
        }

        T = s[0].Inv();
        var Y = X.Select(xi => xi * T).ToKMatrix(N);
        s = s.Select(si => si * T).ToArray();
        var H = new KMatrix<BigReal>(o.Zero, N, N - 1);
        // Console.WriteLine("X");
        // Console.WriteLine(X);
        // Console.WriteLine();
        // Console.WriteLine("Y");
        // Console.WriteLine(Y);
        // Console.WriteLine();
        for (int j = 0; j < N - 1; j++)
        {
            // for (int i = 0; i < j - 1; i++)
            //     H.Coefs[i, j] = o.Zero;

            H.Coefs[j, j] = s[j + 1] / s[j];
            var tj = Y[j, 0] / (s[j] * s[j + 1]);
            for (int i = j + 1; i < N; i++)
                H.Coefs[i, j] = -Y[i, 0] * tj;
        }

        // Console.WriteLine("H");
        // Console.WriteLine(H);
        // Console.WriteLine();
        var YSQ = NSQ.Range().Select(_ => new KMatrix<BigReal>(o, N, 1)).ToArray();
        int IT = 0, ITS = 0, IMQ = 0, IZD = 0, IZM = 0;

        var gamma_pow = (N - 1).Range().Select(i => (i, yi: gamma.Pow(i + 1))).ToArray();

        var z1 = BigReal.BrZero(O2);
        var DYSQ = NSQ.Range().Select(_ => new KMatrix<BigReal>(z1, N, 1)).ToArray();
        var DH = new KMatrix<BigReal>(z1, N, N - 1);
        var DA = new KMatrix<BigReal>(z1, N, N);
        var DB = new KMatrix<BigReal>(z1, N, N);
        var DY = new KMatrix<BigReal>(z1, N, 1);
        var DT = new KMatrix<BigReal>(z1, N, N - 1);
        var (_DH, _DA, _DB, _DT, _DY) = (DH.Clone, DA.Clone, DB.Clone, DT.Clone, DY.Clone);

        BigReal minY, maxY, maxB;
        int IT0 = IT;

        step2:
        // Console.WriteLine("2. Check if min/max absolute value of Y < DREP");
        // 2. Check if min/max absolute value of Y < DREP . This is often true for
        // the first few tens of iterations. If true, go to Step 7 below.
        minY = Y.Min(e => e.Absolute);
        maxY = Y.Max(e => e.Absolute);
        // Console.WriteLine($"IT:{IT} minY:{minY} maxY:{maxY} dynrange:{minY / maxY} DREP:{DREP}");

        if ((minY / maxY) < DREP)
            goto step7;

        step3:
        // Console.WriteLine("3. Initialize DP arrays");
        // 3. Initialize DP arrays:
        // a. Set T = max absolute value of Y ; set DY = Y/T , rounded to DP.
        // b. Set T = max absolute value of H; set DH = H/T , rounded to DP.
        // c. Set DA and DB to identity matrices.
        // d. Set DYSQ array to zeroes.
        var maxH = (N - 1).Range().Select(k => H[k, k]).Max(e => e.Absolute);
        var mhi = maxH.Inv();
        DY = (Y / maxY).Select(c => c.ToBigReal(O2)).ToKMatrix(N);
        DH = H.Select(c => (c * mhi).ToBigReal(O2)).ToKMatrix(N);

        DA = DA.One;
        DB = DB.One;
        for (int i = 0; i < NSQ; i++)
            Array.Copy(DY.Zero.Coefs, DYSQ[i].Coefs, N);

        // 4. Perform an LQ decomposition on DH using DP arithmetic.
        DH = LQpslq(DH);

        step5:
        // Console.WriteLine("5. Perform one multipair PSLQ iteration using DP arithmetic");
        // 5. Perform one multipair PSLQ iteration using DP arithmetic:
        // a. Increment iteration count: IT := IT + 1; set IZD = 0.
        // b. Save the input DA, DB, DH and DY arrays.
        // c. Follow the steps in multipair PSLQ algorithm above to update DY, DA,
        // DB and DH; however, if flag IMQ = 1 from a previous iteration, only
        // select one pair of entries (i.e., set p = 1 in Steps 2, 3 and 4 in the
        // multipair PSLQ algorithm), then set IMQ = 0. It is not necessary to
        // compute the norm bound in the DP iterations.
        // d. If the min absolute value of DY < DEPS, then set IZD = 1. If the
        // max absolute value of DA or DB exceeds 10^13 , but less than 2^52 , then
        // set IZD = 1. If the max absolute value of DA or DB exceeds 2^52
        // (precision failure), then set IZD = 2 and restore the DA, DB, DH
        // and DY arrays saved above in Step 5b.
        // e. Compare the DY vector with the DY vectors of recent iterations saved
        // in array DYSQ; if a match is found, set flag IMQ = 1, which instructs
        // the next iteration to be performed using only one pair of indices in the
        // multipair scheme (in practice, this occurs only very rarely).
        // f. Save DY vector in row K of DYSQ, where K = 1 + mod(IT, NSQ)
        // (i.e., DY vectors are stored a circular sequence in the DYSQ array).

        IZD = 0;
        // Console.WriteLine($"DP IT:{IT}");
        // DY.Println("DY");
        // Console.WriteLine("DH");
        // Console.WriteLine(DH);
        // Console.WriteLine();
        // Console.WriteLine("DA");
        // Console.WriteLine(DA);
        // Console.WriteLine();
        // Console.WriteLine("DB");
        // Console.WriteLine(DB);
        // Console.WriteLine();
        ++IT;
        
        SaveRestore(DH, DA, DB, DT, DY, _DH, _DA, _DB, _DT, _DY);
        IterOneLevelMultipair(DH, DA, DB, DT, DY, gamma_pow, IMQ == 1);
        IMQ = 0;

        minY = DY.Min(e => e.Absolute);
        IZD = minY < DEPS ? 1 : 0;
        var maxDADB = BigReal.Max(DA.Max(e => e.Absolute), DB.Max(e => e.Absolute));
        var log2MaxDADB = double.Log2(maxDADB.ToDouble);
        var log10MaxDADB = double.Log10(maxDADB.ToDouble);
        if (log10MaxDADB > 13 && log2MaxDADB < 52)
        {
            // Console.WriteLine(
            //     $"Found meth 1# maxDADB:{maxDADB} log10MaxDADB:{log10MaxDADB} log2MaxDADB:{log2MaxDADB}");
            IZD = 1;
        }

        if (log2MaxDADB >= 52)
        {
            IZD = 2;
            SaveRestore(_DH, _DA, _DB, _DT, _DY, DH, DA, DB, DT, DY);
        }

        var testDY = DYSQ.Any(dy => BigReal.Sqrt(Ring.SquareNorm2(DY - dy)) < DEPS);
        if (testDY)
        {
            IMQ = 1;
            // Console.WriteLine("###### IMQ");
        }

        Array.Copy(DY.Coefs, DYSQ[IT % NSQ].Coefs, N);
        
        step6:
        // Console.WriteLine("6. Check flags and, if needed, update MPR arrays from DP arrays");
        // 6. Check flags and, if needed, update MPR arrays from DP arrays:
        // a. If IZD = 0, go to Step 5 (continue DP iterations). If IZD = 2,
        // but IT > ITS + 1 (i.e., if the current iteration is more than one plus
        // the previous iteration when a MPR update was performed), then set
        // IZD = 1, so that after an MPR update, regular DP iterations can
        // continue. But if IT = ITS + 1 (i.e., an MPR update was performed
        // on the previous iteration), then leave IZD = 2.
        // b. Update MPR arrays from DP arrays: Set Y := DB x Y (matrix mul-
        // tiplication), then find the min absolute value of the updated Y . Set
        // B := DB x B (matrix multiplication), then find max absolute value of
        // updated B. Set H := DA x H (matrix multiplication). Set ITS = IT.
        // c. The max norm bound may optionally be computed here, as described
        // in the multipair PSLQ algorithm above, and output, along with the
        // current min and max of Y and other data for informational purposes.
        // d. If the min absolute value of Y is less than EPS x max absolute value
        // of B (tentative detection), then set IZM = 1; else if the min absolute
        // value of Y is less than EPS x 2^72 x max absolute value of B (precision
        // exhausted), then set IZM = 2; else set IZM = 0.
        // e. Test output flag: If IZM = 0, then if IZD = 2, go to Step 7 (start
        // MPR iterations), else go to Step 2 (start DP iterations); else if IZM =
        // 1, go to Step 9 (exit); else if IZM = 2, go to Step 9 (exit).
        if (IZD == 0)
            goto step5;

        if (IZD == 2 && IT > ITS + 1)
            IZD = 1;

        var B0 = DB.Select(c => c.ToBigReal(NDP)).ToKMatrix(N);
        var A0 = DA.Select(c => c.ToBigReal(NDP)).ToKMatrix(N);
        Y = B0 * Y;
        B = B0 * B;
        H = A0 * H;

        ITS = IT;
        maxB = B.Max(e => e.Absolute);
        minY = Y.Min(e => e.Absolute);
        if (minY < EPS * maxB)
            IZM = 1;
        else if (minY < EPS * twoPow72 * maxB)
            IZM = 2;
        else
            IZM = 0;

        // Console.WriteLine(
        //     $"IZM:{IZM} minY:{minY} EPS:{EPS} maxB:{maxB} EPS * maxB:{EPS * maxB} EPS * twoPow72 * maxB:{EPS * twoPow72 * maxB}");
        if (IZM == 0)
        {
            if (IZD != 2)
                goto step2;
        }
        else
            goto step9;

        step7:
        // 7. Perform an LQ decomposition on H using MPR arithmetic.
        H = LQpslq(H);
        // Console.WriteLine("H");
        // Console.WriteLine(H);
        // Console.WriteLine();
        for (int i = 0; i < NSQ; i++)
            Array.Copy(Y.Zero.Coefs, YSQ[i].Coefs, N);

        step8:
        // Console.WriteLine("8. Perform one multipair PSLQ iteration using MPR arithmetic");
        // 8. Perform one multipair PSLQ iteration using MPR arithmetic:
        // a. Increment iteration count: IT := IT + 1; set IZM = 0.
        // b. Perform one iteration of the multipair PSLQ algorithm using MPR
        // arithmetic, as described above (except no need to compute norm bound).
        // c. If the min absolute value of Y is less than EPS x max absolute value
        // of B (tentative detection), then set IZM = 1; else if the min absolute
        // value of Y is less than EPS x 2^72 x max absolute value of B (precision
        // exhausted), then set IZM = 2.
        // d. Test output flag: If IZM = 0, then periodically (every IPM iterations
        // since the last MPR update) check if the min/max dynamic range of
        // Y < DREP ; if true, go to Step 8 (continue MPR iterations); else go
        // to Step 3 (start DP iterations); else if IZM = 1, go to Step 9 (exit);
        // else if IZM = 2, go to Step 9 (exit).
        ++IT;
        // Console.WriteLine($"MP IT:{IT}");

        IZM = 0;
        IterOneLevelMultipair(H, A, B, HT, Y, gamma_pow, IMQ == 1);

        var testY = YSQ.Any(y => Ring.SquareNorm2(Y - y) < EPS);
        if (testY)
            IMQ = 1;

        Array.Copy(Y.Coefs, YSQ[IT % NSQ].Coefs, N);

        maxB = B.Max(e => e.Absolute);
        minY = Y.Min(e => e.Absolute);
        maxY = Y.Max(e => e.Absolute);
        if (minY < EPS * twoPow72 * maxB)
        {
            if (minY < EPS * maxB)
            {
                // Console.WriteLine($"Found meth 2.1 minY:{minY} EPS:{EPS} maxB:{maxB} => {EPS * maxB}");
                IZM = 1;
            }
            else
            {
                // Console.WriteLine(
                //     $"Found meth 2.2 minY:{minY} EPS:{EPS} maxB:{maxB} EPS * twoPow72 * maxB:{EPS * twoPow72 * maxB}");
                IZM = 2;
            }
        }

        if (IZM == 0)
        {
            if ((minY / maxY) < DREP)
                goto step8;
            else
                goto step3;
        }

        step9:
        // 9. Exit:
        // a. If IZM = 1 find the index of Y with the min absolute value; set R =
        // row of B corresponding to that index. If the Euclidean norm of the
        // relation is less than 10^NRB , and the dynamic range of the Y vector
        // is at least 10^NDR , set IQ = 1 (success) and exit; otherwise set R =
        // zeroes and set IQ = 0 (failure), then exit.
        // b. If IZM = 2, then set R = zeroes and set IQ = 0 (failure), then exit.
        if (IZM == 2)
            throw new("Failure IZM");

        var (idx, _) = Y.Select((yi, i) => (i, yi)).MinBy(e => e.yi.Absolute);
        var R = B.GetRow(idx);
        var normR = BigReal.Sqrt(Ring.SquareNorm2(R));
        minY = Y.Min(e => e.Absolute);
        maxY = Y.Max(e => e.Absolute);

        Console.WriteLine($"Possible Solution step:{IT}");
        // Console.WriteLine(R);
        // Console.WriteLine(Y);
        // Console.WriteLine(B);
        if (double.Log10(normR.ToDouble) < NRB && double.Log10((minY / maxY).ToDouble) < NDR)
            return R;

        throw new
            ($"Failure normR or (minY / maxY) IT:{IT}");
    }
}