using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

public struct KMatrix<K> : IVsElt<K, KMatrix<K>>, IElt<KMatrix<K>>, IRingElt<KMatrix<K>>, IFieldElt<KMatrix<K>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public K[,] Coefs { get; }
    public int M { get; }
    public int N { get; }
    public (int m, int n) Dim => (M, N);

    public KMatrix(K scalar, int m, int n)
    {
        Coefs = new K[m, n];
        M = m;
        N = n;
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Coefs[i, j] = scalar;
            }
        }

        P = scalar.P;
        KZero = scalar.Zero;
        Hash = (M, N, P).GetHashCode();
    }

    public KMatrix(K[,] coefs)
    {
        M = coefs.GetLength(0);
        N = coefs.GetLength(1);
        Coefs = new K[M, N];
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Coefs[i, j] = coefs[i, j];
            }
        }

        KZero = Coefs[0, 0].Zero;
        P = Coefs[0, 0].P;
        Hash = (M, N, P).GetHashCode();
    }

    public K KZero { get; }
    public K KOne => KZero.One;

    public KMatrix<K> LeadingCoeff => One;
    public K this[int i, int j] => Coefs[i, j];

    public bool Equals(KMatrix<K> other)
    {
        if (!Dim.Equals(other.Dim))
            return false;

        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (!Coefs[i, j].Equals(other[i, j]))
                    return false;
            }
        }

        return true;
    }

    public int CompareTo(KMatrix<K> other)
    {
        var compDim = Dim.CompareTo(other.Dim);
        if (compDim != 0)
            return compDim;

        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                var comp = Coefs[i, j].CompareTo(other[i, j]);
                if (comp != 0)
                    return comp;
            }
        }

        return 0;
    }

    public KMatrix<K> T => new(Ring.Transpose(Coefs));

    public KMatrix<K> KMul(K k)
    {
        var coefs = new K[M, N];
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                coefs[i, j] = k * Coefs[i, j];
            }
        }

        return new(coefs);
    }

    public int Hash { get; }

    public bool IsZero()
    {
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (!Coefs[i, j].IsZero())
                    return false;
            }
        }

        return true;
    }

    public KMatrix<K> Zero => new(KZero, M, N);
    public KMatrix<K> One => new(Ring.Diagonal(KOne, M));

    public KMatrix<K> Add(KMatrix<K> e)
    {
        if (!Dim.Equals(e.Dim))
            throw new GroupException(GroupExceptionType.GroupDef);

        var coefs = new K[M, N];
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                coefs[i, j] = Coefs[i, j] + e[i, j];
            }
        }

        return new(coefs);
    }

    public KMatrix<K> Sub(KMatrix<K> e)
    {
        if (!Dim.Equals(e.Dim))
            throw new GroupException(GroupExceptionType.GroupDef);

        var coefs = new K[M, N];
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                coefs[i, j] = Coefs[i, j] - e[i, j];
            }
        }

        return new(coefs);
    }

    public KMatrix<K> Opp()
    {
        var coefs = new K[M, N];
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                coefs[i, j] = Coefs[i, j].Opp();
            }
        }

        return new(coefs);
    }

    public KMatrix<K> Mul(KMatrix<K> e) => new(Ring.Dot(Coefs, e.Coefs));

    public (KMatrix<K> quo, KMatrix<K> rem) Div(KMatrix<K> e) => (Mul(e.Inv()), Zero);

    public KMatrix<K> Mul(int k)
    {
        var coefs = new K[M, N];
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                coefs[i, j] = Coefs[i, j] * k;
            }
        }

        return new(coefs);
    }

    public int P { get; }

    public KMatrix<K> Inv()
    {
        var e = Ring.ReducedRowsEchelonForm(Coefs);
        return new(e.P);
    }

    public KMatrix<K> Pow(int k)
    {
        if (M != N)
            throw new GroupException(GroupExceptionType.GroupDef);

        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var pi = this;
        return Enumerable.Repeat(pi, k).Aggregate((a, b) => a.Mul(b));
    }

    public (int nullity, KMatrix<K>) NullSpace()
    {
        var e = Ring.ReducedRowsEchelonForm(T);
        var rgM = M.Range();
        var rgN = N.Range();
        var nullRows = rgN.Where(i => rgM.All(j => e.A0[i, j].IsZero())).ToArray();
        var nullity = nullRows.Length;
        if (nullity == 0)
            return (0, new(KZero, M, 0));

        var coeffs = new K[N, nullity];
        var ep = e.P.T;
        for (int k = 0; k < nullity; k++)
        {
            var j = nullRows[k];
            foreach (var i in rgN)
            {
                coeffs[i, k] = ep[i, j];
            }
        }

        return (nullity, new(coeffs));
    }

    public K Det => Ring.DeterminantByPivot(Coefs);

    public K Trace
    {
        get
        {
            if (M != N)
                throw new ArgumentException();

            var acc = KZero;
            for (int i = 0; i < M; i++)
                acc += Coefs[i, i];

            return acc;
        }
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        return Ring.Matrix2String(Coefs);
    }

    public KMatrix<K> GetRow(int i)
    {
        var r = new K[1, N];
        for (int j = 0; j < M; j++)
            r[0, j] = Coefs[i, j];

        return new(r);
    }

    public KMatrix<K> GetCol(int j)
    {
        var c = new K[M, 1];
        for (int i = 0; i < M; i++)
            c[i, 0] = Coefs[i, j];

        return new(c);
    }

    public KMatrix<K>[] Rows => M.Range().Select(GetRow).ToArray();
    public KMatrix<K>[] Cols => N.Range().Select(GetCol).ToArray();

    public KMatrix<K> DotT(KMatrix<K> m) => Mul(m.T);
    public KMatrix<K> TDot(KMatrix<K> m) => T.Mul(m);

    public static KMatrix<K> MergeSameRows(params KMatrix<K>[] matrices)
    {
        var m = matrices[0].M;
        if (matrices.Any(m0 => m0.M != m))
            throw new ArgumentException();

        var n = matrices.Sum(m0 => m0.N);
        var mat0 = new K[m, n];
        var j = 0;
        foreach (var mat in matrices)
        {
            for (int j0 = 0; j0 < mat.N; j0++, j++)
            {
                for (int i = 0; i < m; i++)
                {
                    mat0[i, j] = mat[i, j0];
                }
            }
        }

        return new(mat0);
    }

    public static KMatrix<K> MergeSameCols(params KMatrix<K>[] matrices)
    {
        var n = matrices[0].N;
        if (matrices.Any(m0 => m0.N != n))
            throw new ArgumentException();

        var m = matrices.Sum(m0 => m0.M);
        var mat0 = new K[m, n];
        var i = 0;
        foreach (var mat in matrices)
        {
            for (int i0 = 0; i0 < mat.M; i0++, i++)
            {
                for (int j = 0; j < m; j++)
                {
                    mat0[i, j] = mat[i0, j];
                }
            }
        }

        return new(mat0);
    }

    public static KMatrix<K> operator +(KMatrix<K> a, KMatrix<K> b) => a.Add(b);

    public static KMatrix<K> operator +(int a, KMatrix<K> b) => b.KOne.Mul(a) + b;

    public static KMatrix<K> operator +(KMatrix<K> a, int b) => a + a.KOne.Mul(b);

    public static KMatrix<K> operator -(KMatrix<K> a) => a.Opp();

    public static KMatrix<K> operator -(KMatrix<K> a, KMatrix<K> b) => a.Sub(b);

    public static KMatrix<K> operator -(int a, KMatrix<K> b) => b.KOne.Mul(a) - b;

    public static KMatrix<K> operator -(KMatrix<K> a, int b) => a - a.KOne.Mul(b);

    public static KMatrix<K> operator *(KMatrix<K> a, KMatrix<K> b) => a.Mul(b);

    public static KMatrix<K> operator *(int a, KMatrix<K> b) => b.Mul(a);

    public static KMatrix<K> operator *(KMatrix<K> a, int b) => a.Mul(b);

    public static KMatrix<K> operator /(KMatrix<K> a, KMatrix<K> b) => a.Div(b).quo;

    public static KMatrix<K> operator /(KMatrix<K> a, int b)
    {
        var coefs = new K[a.M, a.N];
        for (int i = 0; i < a.M; i++)
        {
            for (int j = 0; j < a.N; j++)
            {
                coefs[i, j] = a.Coefs[i, j] / b;
            }
        }

        return new(coefs);
    }

    public KMatrix<EPoly<K>> ToEMatrix(KPoly<K> f)
    {
        var m0 = new EPoly<K>[M, N];
        for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            m0[i, j] = new(f, this[i, j] * f.One);

        return new(m0);
    }

    public static KMatrix<K> operator /(int a, KMatrix<K> b) => b.Inv().Mul(a);

    public static KMatrix<K> operator +(KMatrix<K> a, K b) => a.Add(new(b, a.M, a.N));

    public static KMatrix<K> operator +(K a, KMatrix<K> b) => b.Add(new(a, b.M, b.N));

    public static KMatrix<K> operator -(KMatrix<K> a, K b) => a.Sub(new(b, a.M, a.N));

    public static KMatrix<K> operator -(K a, KMatrix<K> b) => new KMatrix<K>(a, b.M, b.N).Sub(b);

    public static KMatrix<K> operator *(KMatrix<K> a, K b) => a.KMul(b);

    public static KMatrix<K> operator *(K a, KMatrix<K> b) => b.KMul(a);

    public static KMatrix<K> operator /(KMatrix<K> a, K b) => a.KMul(b.Inv());

    public static KMatrix<FracPoly<K>> operator +(KMatrix<K> m, KMatrix<FracPoly<K>> a)
    {
        if (!a.Dim.Equals(m.Dim))
            throw new ArgumentOutOfRangeException();

        var m0 = new FracPoly<K>[a.M, a.N];
        for (int i = 0; i < a.M; i++)
        for (int j = 0; j < a.N; j++)
            m0[i, j] = m[i, j] + a[i, j];

        return new(m0);
    }

    public static KMatrix<FracPoly<K>> operator +(KMatrix<FracPoly<K>> a, KMatrix<K> m) => m + a;
    public static KMatrix<FracPoly<K>> operator -(KMatrix<K> m, KMatrix<FracPoly<K>> a) => m + a.Opp();
    public static KMatrix<FracPoly<K>> operator -(KMatrix<FracPoly<K>> a, KMatrix<K> m) => a + m.Opp();
    public static KMatrix<FracPoly<K>> operator *(KMatrix<K> m, FracPoly<K> a)
    {
        var m0 = new FracPoly<K>[m.M, m.N];
        for (int i = 0; i < m.M; i++)
        for (int j = 0; j < m.N; j++)
            m0[i, j] = m[i, j] * a;

        return new(m0);
    }

    public static KMatrix<FracPoly<K>> operator *(FracPoly<K> a, KMatrix<K> m) => m * a;
}