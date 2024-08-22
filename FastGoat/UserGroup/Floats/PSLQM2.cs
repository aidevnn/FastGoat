using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Floats;

public class PSLQM2
{
    /// <summary>
    /// Step Iteration of one level multipair PSLQ
    /// </summary>
    static void IterOneLevelMultipair<K>(KMatrix<K> H, KMatrix<K> A, KMatrix<K> B,
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
                var t0 = (H[mj, mj].Pow(2) + H[mj, mj + 1].Pow(2)).Sqrt();
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

    public static KMatrix<K> LQpslq<K>(KMatrix<K> A) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
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

            var N2 = (m - l).Range().Select(i => H[l, l + i]).Aggregate(A.KZero, (acc, a) => acc + a * a);
            var N = N2.Sqrt();
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
    private int NMP { get; }
    private int NDR { get; }
    private int NRB { get; }
    private int NDP { get; }
    private int NEP { get; }
    private int NSQ { get; }
    private Rational[] Relations { get; set; }
    private KMatrix<BigReal> X { get; }
    private KMatrix<BigReal> H { get; set; }
    private KMatrix<BigReal> A { get; }
    private KMatrix<BigReal> B { get; set; }
    private KMatrix<BigReal> T { get; }
    private KMatrix<BigReal> Y { get; set; }
    private KMatrix<BigReal>[] Yseq { get; }
    private KMatrix<Dcml> DH { get; }
    private KMatrix<Dcml> DA { get; }
    private KMatrix<Dcml> DB { get; }
    private KMatrix<Dcml> DT { get; }
    private KMatrix<Dcml> DY { get; }
    private KMatrix<Dcml>[] DYseq { get; }
    private KMatrix<Dcml> _DH { get; }
    private KMatrix<Dcml> _DA { get; }
    private KMatrix<Dcml> _DB { get; }
    private KMatrix<Dcml> _DT { get; }
    private KMatrix<Dcml> _DY { get; }
    private (int i, BigReal yi)[] gamma_pow { get; }
    private (int i, Dcml yi)[] gamma_pow_dble { get; }
    private BigReal EPS { get; }
    private Dcml DEPS { get; }
    private BigReal DREP { get; }
    private BigReal TwoPowNMP { get; }

    private PSLQM2(KMatrix<BigReal> X, BigReal gamma)
    {
        N = X.N;
        this.X = new KMatrix<BigReal>(X.Coefs);
        var o = X.KOne;
        NMP = gamma.O;
        NDR = NMP / 20;
        NRB = NMP / 3;
        NDP = 26;
        NEP = -(9 * NMP / 10);
        NSQ = 8;
        Dcml.EpsDouble = new Dcml(0.1m).Pow(NDP);
        DEPS = new Dcml(0.1m).Pow(NDP - 4);
        DREP = new Dcml(0.1m).Pow(NDP - 8).ToBigReal(NMP);

        EPS = BigReal.FromBigIntegerAndExponent(1, NEP, NMP);
        TwoPowNMP = BigReal.FromBigInteger(BigInteger.Pow(2, NMP / 5), NMP);
        A = new KMatrix<BigReal>(Ring.Diagonal(o, N));
        B = new KMatrix<BigReal>(Ring.Diagonal(o, N));
        T = new KMatrix<BigReal>(o, N, N - 1);
        Y = new KMatrix<BigReal>(o, N, 1);
        H = new KMatrix<BigReal>(o.Zero, N, N - 1);
        Yseq = NSQ.Range().Select(_ => new KMatrix<BigReal>(o, N, 1)).ToArray();
        gamma_pow = (N - 1).Range().Select(i => (i, yi: gamma.Pow(i + 1))).ToArray();
        gamma_pow_dble = (N - 1).Range().Select(i => (i, yi: gamma.Pow(i + 1).ToDcml)).ToArray();
        var z1 = Dcml.DbleZero();
        DYseq = NSQ.Range().Select(_ => new KMatrix<Dcml>(z1, N, 1)).ToArray();
        DH = new KMatrix<Dcml>(z1, N, N - 1);
        DA = new KMatrix<Dcml>(Ring.Diagonal(z1, N));
        DB = new KMatrix<Dcml>(Ring.Diagonal(z1, N));
        DY = new KMatrix<Dcml>(z1, N, 1);
        DT = new KMatrix<Dcml>(z1, N, N - 1);
        (_DH, _DA, _DB, _DT, _DY) = (DH.Clone, DA.Clone, DB.Clone, DT.Clone, DY.Clone);
        Relations = Array.Empty<Rational>();
    }

    // 1. Initialize MPR arrays:
    // a. Set EPS = 10^NEP .
    // b. Set B to an N * N identity matrix.
    // c. Compute initial H matrix as given in multipair PSLQ algorithm above.
    // d. Set SYQ array to zeroes.
    // e. Set IZD, IZM, IMQ, IT and ITS to zero.
    private PslqState Step1()
    {
        var s = new BigReal[N];
        var t = BigReal.BrZero(NMP);
        for (int i = N - 1; i >= 0; i--)
        {
            t += X[0, i].Pow(2);
            s[i] = BigReal.Sqrt(t);
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
                H.Coefs[i, j] = -Y[i, 0] * tj;
        }

        return new(0, 0, 0, 0, 0, 2);
    }

    // 2. Check if min/max absolute value of Y < DREP . This is often true for
    // the first few tens of iterations. If true, go to Step 7 below.
    private PslqState Step2(PslqState st)
    {
        var minY = Y.Min(e => e.Absolute);
        var maxY = Y.Max(e => e.Absolute);

        if ((minY / maxY) < DREP)
            return new(st.IT, st.ITS, st.IZD, st.IZM, st.IMQ, 7);
        else
            return new(st.IT, st.ITS, st.IZD, st.IZM, st.IMQ, 3);
    }

    // 3. Initialize DP arrays:
    // a. Set T = max absolute value of Y ; set DY = Y/T , rounded to DP.
    // b. Set T = max absolute value of H; set DH = H/T , rounded to DP.
    // c. Set DA and DB to identity matrices.
    // d. Set DYSQ array to zeroes.
    // 4. Perform an LQ decomposition on DH using DP arithmetic.
    private PslqState Step3(PslqState st)
    {
        var maxY = Y.Max(e => e.Absolute);
        var maxH = (N - 1).Range().Select(k => H[k, k]).Max(e => e.Absolute);
        var (myi, mhi) = (maxY.Inv(), maxH.Inv());
        for (int i = 0; i < N; i++)
        {
            DY.Coefs[i, 0] = (myi * Y.Coefs[i, 0]).ToDcml;
            for (int j = 0; j < N - 1; j++)
                DH.Coefs[i, j] = (mhi * H.Coefs[i, j]).ToDcml;
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
        var tmp1 = LQpslq(DH);
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
        SaveRestore(DH, DA, DB, DT, DY, _DH, _DA, _DB, _DT, _DY);
        IterOneLevelMultipair(DH, DA, DB, DT, DY, gamma_pow_dble, st.IMQ == 1);

        var minDY = DY.Min(e => e.Absolute);
        var IZD = minDY < DEPS ? 1 : 0;
        var maxDADB = Dcml.Max(DA.Max(e => e.Absolute), DB.Max(e => e.Absolute));
        var log2MaxDADB = double.Log2((double)maxDADB.K);
        var log10MaxDADB = double.Log10((double)maxDADB.K);
        if (log10MaxDADB >= 13 && log2MaxDADB < 52)
            IZD = 1;

        if (log2MaxDADB >= 52)
        {
            IZD = 2;
            SaveRestore(_DH, _DA, _DB, _DT, _DY, DH, DA, DB, DT, DY);
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

            var B0 = DB.Select(c => c.ToBigReal(NMP)).ToKMatrix(N);
            var A0 = DA.Select(c => c.ToBigReal(NMP)).ToKMatrix(N);
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
            if (minY < EPS * TwoPowNMP * maxB)
            {
                if (minY < EPS * maxB)
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
        var tmpH = LQpslq(H);
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
        var IT = st.IT + 1;
        var IZM = 0;
        var IMQ = st.IMQ;
        IterOneLevelMultipair(H, A, B, T, Y, gamma_pow, IMQ == 1);
        var testY = Yseq.Any(y => (Y - y).Max(e => e.Absolute) < EPS * TwoPowNMP);
        IMQ = testY ? 1 : 0;

        Array.Copy(Y.Coefs, Yseq[IT % NSQ].Coefs, N);

        var maxB = B.Max(e => e.Absolute);
        var minY = Y.Min(e => e.Absolute);
        var maxY = Y.Max(e => e.Absolute);
        IZM = 0;
        if (minY < EPS * TwoPowNMP * maxB)
        {
            if (minY < EPS * maxB)
                IZM = 1;
            else
                IZM = 2;
        }

        if (IZM == 0)
        {
            if ((minY / maxY) < DREP)
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

        var (idx, _) = Y.Select((yi, i) => (i, yi)).MinBy(e => e.yi.Absolute);
        var R = B.GetRow(idx);
        var normR = BigReal.Sqrt(Ring.SquareNorm2(R));
        var minY = Y.Min(e => e.Absolute);
        var maxY = Y.Max(e => e.Absolute);

        if (double.Log10(normR.ToDouble) < NRB && double.Log10((minY / maxY).ToDouble) < NDR)
        {
            if (Logger.Level != LogLevel.Off)
                Console.WriteLine($"Possible Solution step:{st.IT}");
            
            Relations = R.Select(c => c.RoundEven.ToRational).ToArray();
        }
        else
        {
            Console.WriteLine($"Failure normR or (minY / maxY) IT:{st.IT}");
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

    public static Rational[] TwoLevelMultipairXP(KMatrix<BigReal> X, BigReal gamma)
    {
        var algo = new PSLQM2(X, gamma);
        algo.Run();
        return algo.Relations;
    }
}
