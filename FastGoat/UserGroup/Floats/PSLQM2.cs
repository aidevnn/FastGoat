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
pslqm2-alg.pdf
https://www.davidhbailey.com/dhbpapers/

Fortran source code
tpslqm2.f90
mpfun20-mpfr-v32.tar.gz
https://www.davidhbailey.com/dhbsoftware/
 */
public static class PSLQ
{
    /// <summary>
    /// Step Iteration of one level multipair PSLQ
    /// </summary>
    public static void IterOneLevelMultipair<K>(KMatrix<K> H, KMatrix<K> A, KMatrix<K> B,
        KMatrix<K> T, KMatrix<K> y, (int i, K yi)[] gamma_pow, bool imq = false)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
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
        var n = H.M;
        var list = new List<(int, int)>(n);
        var selected = new HashSet<int>(n);
        foreach (var (i, j) in T.M.Range().Grid2D(T.N.Range()))
            T.Coefs[i, j] = i == j ? T.KOne : T.KZero;

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
        // i := mj to n: set t3 := Hi_mj ; t4 := Hi_mj1 ; Hi_mj := t*1t3 + t*2t4 ; and
        // Hi_mj1 := −t2t3 + t1t4 ; endfor; endif; endfor.
        for (int j = 0; j < list.Count; j++)
        {
            var mj = list[j].Item1;
            if (mj < n - 2)
            {
                var (b, c) = (H[mj, mj], H[mj, mj + 1]);
                var t0 = K.Sqrt(b.Absolute2 + c.Absolute2);
                var t1 = b / t0;
                var t2 = c / t0;
                var t1c = t1.Conj;
                var t2c = t2.Conj;
                for (int i = mj; i < n; i++)
                {
                    var t3 = H[i, mj];
                    var t4 = H[i, mj + 1];
                    H.Coefs[i, mj] = t1c * t3 + t2c * t4;
                    H.Coefs[i, mj + 1] = -t2 * t3 + t1 * t4;
                }
            }
        }

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

    public static string ToLongString<K>(K e)
    {
        if (e is BigCplx c)
            return c.ToSciForm();
        if (e is BigReal r)
            return r.ToSciForm();
        if (e is Cplx c1)
            return $"{c1.RealPart:E18} + {c1.ImaginaryPart:E18}*I";
        if (e is Dble dp)
            return $"{dp.K:E18}";
        if (e is Dcml dc)
            return $"{dc.K:G27}";

        throw new();
    }

    /// <summary>
    /// One-level multipair PSLQ algorithm
    /// PSLQM1
    /// </summary>
    /// <param name="x">Row Vector</param>
    /// <param name="gamma">PSLQ gamma</param>
    /// <param name="checkPos">Index of control term</param>
    /// <returns>Integer Relation coefficients</returns>
    /// <exception cref="Exception">Precision exhausted</exception>
    public static K[] OnelevelMultipair<K>(KMatrix<K> x, K gamma, int checkPos)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
    {
        // Initialize:
        // 1. For j := 1 to n: for i := 1 to n: if i = j then set Aij := 1 and Bij := 1
        // else set Aij := 0 and Bij := 0; endfor; endfor.
        var n = x.N;
        var z = x.KZero;
        var A = new KMatrix<K>(z, n, n).One;
        var B = new KMatrix<K>(z, n, n).One;
        var T = new KMatrix<K>(z, n, n - 1);

        // 2. For k := 1 to n: set sk :=Sqrt(Sum[j=k to n] xj x*j) ; endfor; set t = 1/s1 ; for k := 1 to
        // n: set yk := txk ; sk := tsk ; endfor.
        var s = n.Range().Select(j => (n - j).Range(j).Aggregate(x.KZero, (acc, k) => acc + x[0, k].Absolute2))
            .Select(e => K.Sqrt(e)).ToArray();
        var t = s[0].Inv();
        var y = x.Select(xi => xi * t).ToKMatrix(n);
        s = s.Select(si => si * t).ToArray();

        // 3. Initial H: For j := 1 to n − 1: for i := 1 to j − 1: set Hij := 0; endfor;
        // set Hjj := sj+1 /sj ; for i := j + 1 to n: set Hij := −y*i yj /(sj sj+1 ); endfor;
        // endfor.
        var H = new KMatrix<K>(z, n, n - 1);
        for (int j = 0; j < n - 1; j++)
        {
            for (int i = 0; i < j - 1; i++)
                H.Coefs[i, j] = z;

            H.Coefs[j, j] = s[j + 1] / s[j];
            for (int i = j + 1; i < n; i++)
                H.Coefs[i, j] = (-y[i, 0].Conj * y[j, 0]) / (s[j] * s[j + 1]);
        }

        var step = 0;
        var O1 = GetO(z);
        var gamma_pow = (n - 1).Range().Select(i => (i, yi: ToNbDigits(gamma.Pow(i + 1), O1))).ToArray();

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
            if (A.Any(e => !e.IsZero() && GetMantissa(e) >= O1 + 2))
                throw new("Precision is exhausted");

            if (y.Select((yk, k) => (yk, k))
                .Any(e => ToNbDigits(e.yk, O1 - n / 2).IsZero() && !B[e.k, checkPos].IsZero()))
            {
                if (Logger.Level != LogLevel.Off)
                    Console.WriteLine($"Possible Solution step:{step}");

                var (_, idx) = y.Select((yk, k) => (yk, k))
                    .First(e => ToNbDigits(e.yk, O1 - n / 2).IsZero() && !B[e.k, checkPos].IsZero());

                Console.WriteLine("End PSLQ algorithm");
                Console.WriteLine($"Step {step}");
                return B.GetRow(idx).Select(c => c.RoundEven).ToArray();
            }
        }

        throw new("#1 Problem");
    }

    public static KMatrix<K> LQ<K>(KMatrix<K> A) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
    {
        var H = new KMatrix<K>(A.Coefs);
        var (n, m) = A.Dim;
        var min = int.Min(m, n);
        if (n < m)
            throw new("LQ rows >= cols");

        for (int l = 0; l < min; l++)
        {
            if (l == m - 1)
                continue;

            var N2 = (m - l).Range().Select(i => H[l, l + i]).Aggregate(A.KZero, (acc, a) => acc + a.Absolute2);
            var N = K.Sqrt(N2);
            if (N.IsZero())
                continue;

            if (!H[l, l].IsZero())
                N *= H[l, l] / H[l, l].Absolute;

            var Ni = N.Inv();
            for (int i = 0; i < m - l; ++i)
                H.Coefs[l, l + i] *= Ni;

            H.Coefs[l, l] += 1;
            for (int j = l + 1; j < n; ++j)
            {
                var t = -((m - l).Range().Select(i => H[l, l + i].Conj * H[j, l + i])
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

    public static void SaveRestore<T>(KMatrix<T> H0, KMatrix<T> A0, KMatrix<T> B0, KMatrix<T> T0, KMatrix<T> y0,
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

    public static K ToNbDigits<K>(K e, int o) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
    {
        if (e is BigReal r)
            return (dynamic)r.ToBigReal(o);
        if (e is BigCplx c)
            return (dynamic)c.ToBigCplx(o);

        return e;
    }

    public static int GetO<K>(K e) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
    {
        if (e is BigReal r)
            return r.O;
        if (e is BigCplx c)
            return c.O;
        if (e is Dcml d)
            return 28;
        else
            return 17;

        throw new();
    }

    public static int GetMantissa<K>(K e) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
    {
        if (e is BigReal r)
            return r.V;
        if (e is BigCplx c)
            return c.NormInf.V;

        throw new();
    }

    public static K FromFixedPrecision<K, F>(F e, int o)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
        where F : struct, IElt<F>, IRingElt<F>, IFieldElt<F>, IFloatElt<F>, IFixedPrecisionElt<F>
    {
        if (e is Dble dp)
            return (dynamic)BigReal.FromFixedPrecision(dp, o);
        if (e is Dcml dc)
            return (dynamic)BigReal.FromFixedPrecision(dc, o);
        if (e is Cplx c)
        {
            return (dynamic)BigCplx.FromBigReal(BigReal.FromDouble(c.RealPart, o),
                BigReal.FromDouble(c.ImaginaryPart, o));
        }

        throw new();
    }
}

public class PSLQM2<F, G>
    where F : struct, IElt<F>, IRingElt<F>, IFieldElt<F>, IFloatElt<F>, IFixedPrecisionElt<F>
    where G : struct, IElt<G>, IRingElt<G>, IFieldElt<G>, IFloatElt<G>
{
    private readonly struct PslqState(int it, int its, int izd, int izm, int imq, int step)
    {
        // return new(IT, ITS, IZD, IZM, IMQ, STEP);
        public int IT { get; } = it;
        public int ITS { get; } = its;
        public int IZD { get; } = izd;
        public int IZM { get; } = izm;
        public int IMQ { get; } = imq;
        public int STEP { get; } = step;

        public override string ToString()
        {
            return $"State IT:{IT,-4} ITS:{ITS,-4} IZD:{IZD,-4} IZM:{IZM,-4} IMQ:{IMQ,-4} STEP:{STEP,-4}";
        }
    }

    private int N { get; }
    private int CheckPos { get; }
    private int NMP { get; }
    private int NDR { get; }
    private int NRB { get; }
    private int NDP { get; }
    private int NEP { get; }
    private int NSQ { get; }
    private G[] Relations { get; set; }
    private KMatrix<G> X { get; }
    private KMatrix<G> H { get; set; }
    private KMatrix<G> A { get; }
    private KMatrix<G> B { get; set; }
    private KMatrix<G> T { get; }
    private KMatrix<G> Y { get; set; }
    private KMatrix<G>[] Yseq { get; }
    private KMatrix<F> DH { get; }
    private KMatrix<F> DA { get; }
    private KMatrix<F> DB { get; }
    private KMatrix<F> DT { get; }
    private KMatrix<F> DY { get; }
    private KMatrix<F>[] DYseq { get; }
    private KMatrix<F> _DH { get; }
    private KMatrix<F> _DA { get; }
    private KMatrix<F> _DB { get; }
    private KMatrix<F> _DT { get; }
    private KMatrix<F> _DY { get; }
    private (int i, G yi)[] gamma_pow { get; }
    private (int i, F yi)[] gamma_pow_dcml { get; }
    private G EPS { get; }
    private F DEPS { get; }
    private G DREP { get; }
    private G TwoPowNMP { get; }

    // PSLQM2test (two-level multipair PSLQ):
    // Input arguments and parameters:
    // N            Int     Length of input vector X and output relation vector R.
    // NDP          Int     Double Precision level in digits.
    // NMP          Int     Multi Precision level in digits.
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
    private PSLQM2(KMatrix<G> X, G gamma, int checkPos)
    {
        N = X.N;
        this.X = new KMatrix<G>(X.Coefs);
        var o = X.KOne;
        CheckPos = checkPos;
        NMP = PSLQ.GetO(o);
        NDP = F.Digits;
        if (NMP < NDP + 3)
            throw new ArgumentException($"{NMP}-digits must be >= {NDP}");

        NDR = F.Digits * 2;
        NRB = NMP;
        NEP = -(9 * NMP / 10);
        NSQ = 8;
        DEPS = F.From((10 * o).Pow(3 - NDP));
        DREP = (10 * o).Pow(7 - NDP);

        EPS = (10 * o).Pow(NEP);
        TwoPowNMP = (2 * o).Pow(NMP / 5);
        A = new KMatrix<G>(Ring.Diagonal(o, N));
        B = new KMatrix<G>(Ring.Diagonal(o, N));
        T = new KMatrix<G>(o, N, N - 1);
        Y = new KMatrix<G>(o, N, 1);
        H = new KMatrix<G>(o.Zero, N, N - 1);
        Yseq = NSQ.Range().Select(_ => new KMatrix<G>(o, N, 1)).ToArray();
        gamma_pow = (N - 1).Range().Select(i => (i, yi: PSLQ.ToNbDigits(gamma.Pow(i + 1), NMP))).ToArray();
        gamma_pow_dcml = gamma_pow.Select(e => (e.i, yi: F.From(e.yi))).ToArray();
        var z1 = F.From(o);
        DYseq = NSQ.Range().Select(_ => new KMatrix<F>(z1, N, 1)).ToArray();
        DH = new KMatrix<F>(z1, N, N - 1);
        DA = new KMatrix<F>(Ring.Diagonal(z1, N));
        DB = new KMatrix<F>(Ring.Diagonal(z1, N));
        DY = new KMatrix<F>(z1, N, 1);
        DT = new KMatrix<F>(z1, N, N - 1);
        (_DH, _DA, _DB, _DT, _DY) = (DH.Clone, DA.Clone, DB.Clone, DT.Clone, DY.Clone);
        Relations = Array.Empty<G>();
    }

    // 1. Initialize MPR arrays:
    // a. Set EPS = 10^NEP .
    // b. Set B to an N * N identity matrix.
    // c. Compute initial H matrix as given in multipair PSLQ algorithm above.
    // d. Set SYQ array to zeroes.
    // e. Set IZD, IZM, IMQ, IT and ITS to zero.
    private PslqState Step1()
    {
        var s = new G[N];
        var t = X.KZero;
        for (int i = N - 1; i >= 0; i--)
        {
            t += X[0, i].Absolute2;
            s[i] = G.Sqrt(t);
        }

        t = s[0].Inv();
        for (int i = 0; i < N; i++)
        {
            Y.Coefs[i, 0] = X[0, i] * t;
            s[i] *= t;
        }

        for (int j = 0; j < N - 1; j++)
        {
            H.Coefs[j, j] = s[j + 1] / s[j];
            var tj = Y[j, 0] / (s[j] * s[j + 1]);
            for (int i = j + 1; i < N; i++)
                H.Coefs[i, j] = -Y[i, 0].Conj * tj;
        }

        return new(0, 0, 0, 0, 0, 2);
    }

    // 2. Check if min/max absolute value of Y < DREP . This is often true for
    // the first few tens of iterations. If true, go to Step 7 below.
    private PslqState Step2(PslqState st)
    {
        var minY = Y.Min(e => e.Absolute);
        var maxY = Y.Max(e => e.Absolute);

        if (minY.CompareTo(DREP * maxY) == -1)
            return new(st.IT, st.ITS, st.IZD, st.IZM, st.IMQ, 7);
        else
            return new(st.IT, st.ITS, st.IZD, st.IZM, st.IMQ, 3);
    }

    // 3. Initialize DP arrays:
    // a. Set T = max absolute value of Y ; set DY = Y/T , rounded to DP.
    // b. Set T = max absolute value of H; set DH = H/T , rounded to DP.
    // c. Set DA and DB to identity matrices.
    // d. Set DYSQ array to zeroes.
    private PslqState Step3(PslqState st)
    {
        var maxY = Y.Max(e => e.Absolute);
        var maxH = (N - 1).Range().Select(k => H[k, k]).Max(e => e.Absolute);
        var (myi, mhi) = (maxY.Inv(), maxH.Inv());
        for (int i = 0; i < N; i++)
        {
            DY.Coefs[i, 0] = F.From(myi * Y.Coefs[i, 0]);
            for (int j = 0; j < N - 1; j++)
                DH.Coefs[i, j] = F.From(mhi * H.Coefs[i, j]);
        }

        foreach (var (i, j) in N.Range().Grid2D())
        {
            if (i == j)
                DA.Coefs[i, j] = DB.Coefs[i, j] = DA.KOne;
            else
                DA.Coefs[i, j] = DB.Coefs[i, j] = DA.KZero;
        }

        var tmp0 = DY.Zero;
        for (int i = 0; i < NSQ; i++)
            Array.Copy(tmp0.Coefs, DYseq[i].Coefs, N);

        return new(st.IT, st.ITS, st.IZD, st.IZM, st.IMQ, 4);
    }

    // 4. Perform an LQ decomposition on DH using DP arithmetic.
    private PslqState Step4(PslqState st)
    {
        var tmp1 = PSLQ.LQ(DH);
        foreach (var (i, j) in N.Range().Grid2D((N - 1).Range()))
            DH.Coefs[i, j] = tmp1[i, j];

        return new(st.IT, st.ITS, st.IZD, st.IZM, st.IMQ, 5);
    }

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
    private PslqState Step5(PslqState st)
    {
        var IT = st.IT + 1;
        PSLQ.SaveRestore(DH, DA, DB, DT, DY, _DH, _DA, _DB, _DT, _DY);
        PSLQ.IterOneLevelMultipair(DH, DA, DB, DT, DY, gamma_pow_dcml, st.IMQ == 1);

        var minDY = DY.Min(e => e.Absolute);
        var IZD = minDY < DEPS ? 1 : 0;
        var maxDADB = DA.Concat(DB).Max(e => F.Abs(e));
        var log2MaxDADB = double.Log2(maxDADB);
        var log10MaxDADB = double.Log10(maxDADB);
        if (log10MaxDADB >= 13 && log2MaxDADB < 52)
            IZD = 1;

        if (log2MaxDADB >= 52)
        {
            IZD = 2;
            PSLQ.SaveRestore(_DH, _DA, _DB, _DT, _DY, DH, DA, DB, DT, DY);
        }

        var testDY = DYseq.Any(dy => (DY - dy).Max(e => e.Absolute) < DEPS);
        var IMQ = testDY ? 1 : 0;

        Array.Copy(DY.Coefs, DYseq[IT % NSQ].Coefs, N);
        return new(IT, st.ITS, IZD, st.IZM, IMQ, 6);
    }

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
    private PslqState Step6(PslqState st)
    {
        var IZD = st.IZD;
        if (IZD == 0)
            return new(st.IT, st.ITS, IZD, st.IZM, st.IMQ, 5);
        else
        {
            if (IZD == 2 && st.IT > st.ITS + 1)
                IZD = 1;

            var B0 = DB.Select(e => double.IsNaN(F.Abs(e)) ? e.One : e)
                .Select(c => PSLQ.FromFixedPrecision<G, F>(c, NMP)).ToKMatrix(N);
            var A0 = DA.Select(e => double.IsNaN(F.Abs(e)) ? e.One : e)
                .Select(c => PSLQ.FromFixedPrecision<G, F>(c, NMP)).ToKMatrix(N);
            var tmpY = B0 * Y;
            var tmpB = B0 * B;
            var tmpH = A0 * H;
            for (int i = 0; i < N; i++)
            {
                Y.Coefs[i, 0] = tmpY[i, 0];
                for (int j = 0; j < N; j++)
                {
                    B.Coefs[i, j] = tmpB[i, j];
                    if (j < N - 1)
                        H.Coefs[i, j] = tmpH[i, j];
                }
            }

            var ITS = st.IT;
            var IZM = 0;
            var maxB = B.Max(e => e.Absolute);
            var minY = Y.Min(e => e.Absolute);
            IZM = 0;
            var lt = Y.Select((yk, k) => (yk, k))
                .Where(e => e.yk.CompareTo(EPS * TwoPowNMP * maxB) == -1 && !B[e.k, CheckPos].IsZero())
                .ToArray();
            if (minY.CompareTo(EPS * TwoPowNMP * maxB) == -1 && lt.Length != 0)
            {
                if (minY.CompareTo(EPS * maxB) == -1)
                    IZM = 1;
                else
                    IZM = 2;
            }

            if (IZM == 0)
            {
                if (IZD == 2)
                    return new(st.IT, ITS, IZD, IZM, st.IMQ, 7);
                else
                    return new(st.IT, ITS, IZD, IZM, st.IMQ, 2);
            }
            else
                return new(st.IT, ITS, IZD, IZM, st.IMQ, 9);
        }
    }

    // 7. Perform an LQ decomposition on H using MPR arithmetic.
    private PslqState Step7(PslqState st)
    {
        var tmpH = PSLQ.LQ(H);
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N - 1; j++)
            H.Coefs[i, j] = tmpH[i, j];

        var tmp0 = Y.Zero;
        for (int i = 0; i < NSQ; i++)
            Array.Copy(tmp0.Coefs, Yseq[i].Coefs, N);

        return new(st.IT, st.ITS, st.IZD, st.IZM, st.IMQ, 8);
    }

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
    private PslqState Step8(PslqState st)
    {
        var IMQ = st.IMQ;
        var IZM = 0;
        
        var IT = st.IT + 1;
        PSLQ.IterOneLevelMultipair(H, A, B, T, Y, gamma_pow, IMQ == 1);
        var testY = Yseq.Any(y => (Y - y).Max(e => e.Absolute).CompareTo(EPS * TwoPowNMP) == -1);
        IMQ = testY ? 1 : 0;

        Array.Copy(Y.Coefs, Yseq[IT % NSQ].Coefs, N);

        var maxB = B.Max(e => e.Absolute);
        var minY = Y.Min(e => e.Absolute);
        var maxY = Y.Max(e => e.Absolute);
        IZM = 0;
        var lt = Y.Select((yk, k) => (yk, k))
            .Where(e => e.yk.CompareTo(EPS * TwoPowNMP * maxB) == -1 && !B[e.k, CheckPos].IsZero())
            .ToArray();
        if (minY.CompareTo(EPS * TwoPowNMP * maxB) == -1 && lt.Length != 0)
        {
            if (minY.CompareTo(EPS * maxB) == -1)
                IZM = 1;
            else
                IZM = 2;
        }

        if (IZM == 0)
        {
            if ((minY / maxY).CompareTo(DREP) == -1)
                return new(IT, st.ITS, st.IZD, IZM, IMQ, 8);
            else
                return new(IT, st.ITS, st.IZD, IZM, IMQ, 3);
        }
        else
            return new(IT, st.ITS, st.IZD, IZM, IMQ, 9);
    }

    // 9. Exit:
    // a. If IZM = 1 find the index of Y with the min absolute value; set R =
    // row of B corresponding to that index. If the Euclidean norm of the
    // relation is less than 10^NRB , and the dynamic range of the Y vector
    // is at least 10^NDR , set IQ = 1 (success) and exit; otherwise set R =
    // zeroes and set IQ = 0 (failure), then exit.
    // b. If IZM = 2, then set R = zeroes and set IQ = 0 (failure), then exit.
    private void Step9(PslqState st)
    {
        if (st.IZM == 2)
        {
            Console.WriteLine("Failure IZM");
            return;
        }

        var maxB = B.Max(e => e.Absolute);
        var lt = Y.Select((yk, k) => (yk, k))
            .Where(e => e.yk.CompareTo(EPS * TwoPowNMP * maxB) == -1 && !B[e.k, CheckPos].IsZero())
            .OrderBy(e =>
                B.N.Range().Where(i => i != CheckPos).Select(i => B[e.k, i].Absolute).Aggregate((a0, a1) => a0 + a1))
            .ThenBy(e => e.yk.Absolute)
            .ToArray();
        
        var idx = lt[0].k;
        var R = B.GetRow(idx);
        var normR = Ring.Norm2f(R);
        var minY = Y.Min(e => e.Absolute);
        var maxY = Y.Max(e => e.Absolute);
        
        if (double.Log10(G.Abs(normR)) < NRB && double.Log10(G.Abs(minY / maxY)) < NDR)
        {
            if (Logger.Level != LogLevel.Off)
            {
                Console.WriteLine("End PSLQ algorithm");
                Console.WriteLine($"Possible Solution step:{st.IT}");
            }
            
            Relations = R.Select(c => c.RoundEven).ToArray();
        }
        else
        {
            Console.WriteLine("B");
            Console.WriteLine(B);
            Console.WriteLine("Y");
            Y.Select(e => G.Abs(e)).Println();
            Console.WriteLine($"Failure normR or (minY / maxY) IT:{st.IT}");
            Console.WriteLine($"log(normR):{double.Log10(G.Abs(normR))} ? {NRB}");
            Console.WriteLine($"log(minY / maxY):{double.Log10(G.Abs(minY / maxY) )} ? {NDR}");
            Console.WriteLine(double.Log10(G.Abs(normR)) < NRB);
            Console.WriteLine(double.Log10(G.Abs(minY / maxY)) - 1 < NDR);
        }
    }
    
    private void Run()
    {
        var st = Step1();
        while (st.STEP != 9)
        {
            // Console.WriteLine(st);
            if (st.STEP == 2)
                st = Step2(st);
            else if (st.STEP == 3)
                st = Step3(st);
            else if (st.STEP == 4)
                st = Step4(st);
            else if (st.STEP == 5)
                st = Step5(st);
            else if (st.STEP == 6)
                st = Step6(st);
            else if (st.STEP == 7)
                st = Step7(st);
            else if (st.STEP == 8)
                st = Step8(st);
            else
                throw new("");
        }

        Step9(st);
    }

    /// <summary>
    /// Two-level multipair PSLQ algorithm
    /// PSLQM2test
    /// </summary>
    /// <param name="X">Row Vector</param>
    /// <param name="gamma">PSLQ gamma</param>
    /// <param name="checkPos">Index not null value</param>
    /// <returns>Integer Relation coefficients</returns>
    /// <exception cref="Exception">Precision exhausted</exception>
    public static G[] TwoLevelMultipair(KMatrix<G> X, G gamma, int checkPos)
    {
        if (F.IsComplex != G.IsComplex)
            throw new();

        var algo = new PSLQM2<F, G>(X, gamma, checkPos);
        algo.Run();
        return algo.Relations.ToArray();
    }
}