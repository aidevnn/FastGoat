using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Matrix;

public class GL : IGroup<Mat>
{
    public GL(int n, int p)
    {
        if (n > 5)
            throw new GroupException(GroupExceptionType.GroupDef);

        if (IntExt.Coprimes(p).Count() != p - 1)
            throw new GroupException(GroupExceptionType.GroupDef);

        N = n;
        P = p;
        UnInvertible = IntExt.UnInvertible(p);
        Name = $"GL({n},{p})";
        Fmt = $"{{0,{P.ToString().Length}}}";
        _cache = new int[N * N];
        TableNeutral = new int[N * N];
        for (int i = 0; i < N; ++i)
            TableNeutral[i * (N + 1)] = 1;

        HashNeutral = IntExt.GenHash(P, TableNeutral);
        Hash = (N, P).GetHashCode();
    }

    private int[] _cache;
    public int HashNeutral { get; }
    public int[] TableNeutral { get; }
    private Dictionary<int, int> UnInvertible { get; }
    public int N { get; }
    public int P { get; }
    public string Fmt { get; }

    public IEnumerator<Mat> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<Mat>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public Mat Create(params int[] table)
    {
        var table0 = MatrixExt.ModP(P, table);
        var hash = IntExt.GenHash(P, table0);
        return new(this, hash, table0);
    }

    public Mat At(int[] table, Tuple2Array at, Tuple2Array value)
    {
        if (table.Length != N * N || at.Table.Length != value.Table.Length || at.Table.Any(i => i < 0 || i > N * N))
            throw new GroupException(GroupExceptionType.GroupDef);

        var table0 = MatrixExt.ModP(P, table);
        for (int i = 0; i < at.Table.Length; i++)
        {
            table0[at.Table[i]] = value.Table[i] % P;
        }

        return Create(table0);
    }

    public Mat At(Tuple2Array at, Tuple2Array value) => At(new int[N * N], at, value);

    public Mat At(Tuple2Array at, int value) =>
        At(at, new Tuple2Array(Enumerable.Repeat(value, at.Table.Length).ToArray()));

    public Mat this[params ValueType[] us]
    {
        get
        {
            if (us.Length != N * N)
                throw new GroupException(GroupExceptionType.GroupDef);

            var table = us.Select(i => (int)i).Select(i => i % P).ToArray();
            var hash = IntExt.GenHash(P, table);
            return new(this, hash, table);
        }
    }

    public IEnumerable<Mat> GetElements()
    {
        yield return Neutral();
    }

    public IEnumerable<Mat> GetGenerators()
    {
        yield return Neutral();
    }

    public Mat Neutral() => new(this);

    public Mat Invert(Mat e)
    {
        if (!Equals(e.GL))
            throw new GroupException(GroupExceptionType.BaseGroup);

        if (N == 2)
        {
            var hash = Inv2x2(e.Table, _cache);
            return new(this, hash, _cache);
        }

        if (N == 3)
        {
            var hash = Inv3x3(e.Table, _cache);
            return new(this, hash, _cache);
        }

        if (N == 4)
        {
            var hash = Inv4x4(e.Table, _cache);
            return new(this, hash, _cache);
        }

        if (N == 5)
        {
            var hash = Inv5x5(e.Table, _cache);
            return new(this, hash, _cache);
        }

        throw new GroupException(GroupExceptionType.GroupDef);
    }

    public Mat Op(Mat e1, Mat e2)
    {
        if (!e1.GL.Equals(e2.GL))
            throw new GroupException(GroupExceptionType.BaseGroup);

        var hash = MatrixProduct(P, N, e1.Table, e2.Table, _cache);
        return new(this, hash, _cache);
    }

    public override string ToString() => Name;
    public override int GetHashCode() => Hash;

    int ModP(int e)
    {
        var r = e % P;
        return r < 0 ? r + P : r;
    }

    int MatrixProduct(int p, int aRows, int[] a, int[] b, int[] c)
    {
        var aCols = a.Length / aRows;
        int bRows = aCols;
        var bCols = b.Length / bRows;
        if (a.Length != aRows * aCols || b.Length != bRows * bCols || c.Length != aRows * bCols)
            throw new Exception();

        int hash = 0;
        int pow = 1;
        for (int i = 0; i < aRows; i++)
        {
            for (int j = 0; j < bCols; j++)
            {
                int sum = 0;
                for (int k = 0; k < aCols; k++)
                    sum += a[i * aCols + k] * b[k * bCols + j];

                var v = c[i * aRows + j] = ModP(sum);
                hash += v * pow;
                pow *= p;
            }
        }

        return hash;
    }

    public int Det(Mat e)
    {
        if (N == 2)
            return Det2x2(e.Table);
        if (N == 3)
            return Det3x3(e.Table);
        if (N == 4)
            return Det4x4(e.Table);
        if (N == 5)
            return Det5x5(e.Table);

        throw new();
    }

    int Det2x2(int[] mat)
    {
        var (a, b) = (mat[0], mat[1]);
        var (c, d) = (mat[2], mat[3]);
        var det = (a * d - c * b);
        return ModP(det);
    }

    int Det3x3(int[] mat)
    {
        var (a, b, c) = (mat[0], mat[1], mat[2]);
        var (d, e, f) = (mat[3], mat[4], mat[5]);
        var (g, h, i) = (mat[6], mat[7], mat[8]);
        var det = (a * e * i + d * h * c + g * b * f - a * f * h - b * d * i - c * e * g);
        return ModP(det);
    }

    int Det4x4(int[] mat)
    {
        var (a, b, c, d) = (mat[0], mat[1], mat[2], mat[3]);
        var (e, f, g, h) = (mat[4], mat[5], mat[6], mat[7]);
        var (i, j, k, l) = (mat[8], mat[9], mat[10], mat[11]);
        var (m, n, o, p) = (mat[12], mat[13], mat[14], mat[15]);

        // mathematica.wolframcloud.com
        var det = (d * g * j * m - c * h * j * m - d * f * k * m + b * h * k * m + c * f * l * m - b * g * l * m -
            d * g * i * n + c * h * i * n + d * e * k * n - a * h * k * n - c * e * l * n + a * g * l * n +
            d * f * i * o - b * h * i * o - d * e * j * o + a * h * j * o + b * e * l * o - a * f * l * o -
            c * f * i * p + b * g * i * p + c * e * j * p - a * g * j * p - b * e * k * p + a * f * k * p);
        return ModP(det);
    }

    int Det5x5(int[] mat)
    {
        var (a, b, c, d, e) = (mat[0], mat[1], mat[2], mat[3], mat[4]);
        var (f, g, h, i, j) = (mat[5], mat[6], mat[7], mat[8], mat[9]);
        var (k, l, m, n, o) = (mat[10], mat[11], mat[12], mat[13], mat[14]);
        var (p, q, r, s, t) = (mat[15], mat[16], mat[17], mat[18], mat[19]);
        var (u, v, w, x, y) = (mat[20], mat[21], mat[22], mat[23], mat[24]);

        // mathematica.wolframcloud.com
        var det = e * i * m * q * u - d * j * m * q * u - e * h * n * q * u + c * j * n * q * u + d * h * o * q * u -
            c * i * o * q * u - e * i * l * r * u + d * j * l * r * u + e * g * n * r * u -
            b * j * n * r * u - d * g * o * r * u + b * i * o * r * u + e * h * l * s * u - c * j * l * s * u -
            e * g * m * s * u + b * j * m * s * u + c * g * o * s * u - b * h * o * s * u -
            d * h * l * t * u + c * i * l * t * u + d * g * m * t * u - b * i * m * t * u - c * g * n * t * u +
            b * h * n * t * u - e * i * m * p * v + d * j * m * p * v + e * h * n * p * v -
            c * j * n * p * v - d * h * o * p * v + c * i * o * p * v + e * i * k * r * v - d * j * k * r * v -
            e * f * n * r * v + a * j * n * r * v + d * f * o * r * v - a * i * o * r * v -
            e * h * k * s * v + c * j * k * s * v + e * f * m * s * v - a * j * m * s * v - c * f * o * s * v +
            a * h * o * s * v + d * h * k * t * v - c * i * k * t * v - d * f * m * t * v +
            a * i * m * t * v + c * f * n * t * v - a * h * n * t * v + e * i * l * p * w - d * j * l * p * w -
            e * g * n * p * w + b * j * n * p * w + d * g * o * p * w - b * i * o * p * w -
            e * i * k * q * w + d * j * k * q * w + e * f * n * q * w - a * j * n * q * w - d * f * o * q * w +
            a * i * o * q * w + e * g * k * s * w - b * j * k * s * w - e * f * l * s * w +
            a * j * l * s * w + b * f * o * s * w - a * g * o * s * w - d * g * k * t * w + b * i * k * t * w +
            d * f * l * t * w - a * i * l * t * w - b * f * n * t * w + a * g * n * t * w - e * h * l * p * x +
            c * j * l * p * x + e * g * m * p * x - b * j * m * p * x - c * g * o * p * x + b * h * o * p * x +
            e * h * k * q * x -
            c * j * k * q * x - e * f * m * q * x + a * j * m * q * x + c * f * o * q * x - a * h * o * q * x -
            e * g * k * r * x + b * j * k * r * x + e * f * l * r * x -
            a * j * l * r * x - b * f * o * r * x + a * g * o * r * x + c * g * k * t * x - b * h * k * t * x -
            c * f * l * t * x + a * h * l * t * x + b * f * m * t * x -
            a * g * m * t * x + d * h * l * p * y - c * i * l * p * y - d * g * m * p * y + b * i * m * p * y +
            c * g * n * p * y - b * h * n * p * y - d * h * k * q * y +
            c * i * k * q * y + d * f * m * q * y - a * i * m * q * y - c * f * n * q * y + a * h * n * q * y +
            d * g * k * r * y - b * i * k * r * y - d * f * l * r * y +
            a * i * l * r * y + b * f * n * r * y - a * g * n * r * y - c * g * k * s * y + b * h * k * s * y +
            c * f * l * s * y - a * h * l * s * y - b * f * m * s * y + a * g * m * s * y;
        return ModP(det);
    }

    int Inv2x2(int[] mat, int[] inv)
    {
        var det = Det2x2(mat);
        var idet = UnInvertible[det];
        var (a, b) = (mat[0], mat[1]);
        var (c, d) = (mat[2], mat[3]);

        var a1 = inv[0] = ModP(d * idet);
        var b1 = inv[1] = ModP(-b * idet);
        var c1 = inv[2] = ModP(-c * idet);
        var d1 = inv[3] = ModP(a * idet);

        return a1 + P * (b1 + P * (c1 + P * d1));
    }

    int Inv3x3(int[] mat, int[] inv)
    {
        var det = Det3x3(mat);
        var idet = UnInvertible[det];
        var (a, b, c) = (mat[0], mat[1], mat[2]);
        var (d, e, f) = (mat[3], mat[4], mat[5]);
        var (g, h, i) = (mat[6], mat[7], mat[8]);

        // mathematica.wolframcloud.com
        var a1 = inv[0] = ModP((-f * h + e * i) * idet);
        var b1 = inv[1] = ModP((c * h - b * i) * idet);
        var c1 = inv[2] = ModP((-c * e + b * f) * idet);
        var d1 = inv[3] = ModP((f * g - d * i) * idet);
        var e1 = inv[4] = ModP((-c * g + a * i) * idet);
        var f1 = inv[5] = ModP((c * d - a * f) * idet);
        var g1 = inv[6] = ModP((-e * g + d * h) * idet);
        var h1 = inv[7] = ModP((b * g - a * h) * idet);
        var i1 = inv[8] = ModP((-b * d + a * e) * idet);

        return a1 + P * (b1 + P * (c1 + P * (d1 + P * (e1 + P * (f1 + P * (g1 + P * (h1 + P * i1)))))));
    }

    int Inv4x4(int[] mat, int[] inv)
    {
        var det = Det4x4(mat);
        var idet = UnInvertible[det];
        var (a, b, c, d) = (mat[0], mat[1], mat[2], mat[3]);
        var (e, f, g, h) = (mat[4], mat[5], mat[6], mat[7]);
        var (i, j, k, l) = (mat[8], mat[9], mat[10], mat[11]);
        var (m, n, o, p) = (mat[12], mat[13], mat[14], mat[15]);

        // mathematica.wolframcloud.com
        var a1 = inv[0] = ModP((-h * k * n + g * l * n + h * j * o - f * l * o - g * j * p + f * k * p) * idet);
        var b1 = inv[1] = ModP((d * k * n - c * l * n - d * j * o + b * l * o + c * j * p - b * k * p) * idet);
        var c1 = inv[2] = ModP((-d * g * n + c * h * n + d * f * o - b * h * o - c * f * p + b * g * p) * idet);
        var d1 = inv[3] = ModP((d * g * j - c * h * j - d * f * k + b * h * k + c * f * l - b * g * l) * idet);
        var e1 = inv[4] = ModP((h * k * m - g * l * m - h * i * o + e * l * o + g * i * p - e * k * p) * idet);
        var f1 = inv[5] = ModP((-d * k * m + c * l * m + d * i * o - a * l * o - c * i * p + a * k * p) * idet);
        var g1 = inv[6] = ModP((d * g * m - c * h * m - d * e * o + a * h * o + c * e * p - a * g * p) * idet);
        var h1 = inv[7] = ModP((-d * g * i + c * h * i + d * e * k - a * h * k - c * e * l + a * g * l) * idet);
        var i1 = inv[8] = ModP((-h * j * m + f * l * m + h * i * n - e * l * n - f * i * p + e * j * p) * idet);
        var j1 = inv[9] = ModP((d * j * m - b * l * m - d * i * n + a * l * n + b * i * p - a * j * p) * idet);
        var k1 = inv[10] = ModP((-d * f * m + b * h * m + d * e * n - a * h * n - b * e * p + a * f * p) * idet);
        var l1 = inv[11] = ModP((d * f * i - b * h * i - d * e * j + a * h * j + b * e * l - a * f * l) * idet);
        var m1 = inv[12] = ModP((g * j * m - f * k * m - g * i * n + e * k * n + f * i * o - e * j * o) * idet);
        var n1 = inv[13] = ModP((-c * j * m + b * k * m + c * i * n - a * k * n - b * i * o + a * j * o) * idet);
        var o1 = inv[14] = ModP((c * f * m - b * g * m - c * e * n + a * g * n + b * e * o - a * f * o) * idet);
        var p1 = inv[15] = ModP((-c * f * i + b * g * i + c * e * j - a * g * j - b * e * k + a * f * k) * idet);

        return a1 + P * (b1 + P * (c1 + P * (d1 + P * (e1 + P * (f1 + P * (g1 + P * (h1 + P * (i1 + P * (j1 + P * (k1 +
            P * (l1 + P * (m1 + P * (n1 + P * (o1 + P * p1))))))))))))));
    }

    int Inv5x5(int[] mat, int[] inv)
    {
        var det = Det5x5(mat);
        var idet = UnInvertible[det];
        var (a, b, c, d, e) = (mat[0], mat[1], mat[2], mat[3], mat[4]);
        var (f, g, h, i, j) = (mat[5], mat[6], mat[7], mat[8], mat[9]);
        var (k, l, m, n, o) = (mat[10], mat[11], mat[12], mat[13], mat[14]);
        var (p, q, r, s, t) = (mat[15], mat[16], mat[17], mat[18], mat[19]);
        var (u, v, w, x, y) = (mat[20], mat[21], mat[22], mat[23], mat[24]);

        // mathematica.wolframcloud.com
        var a1 = inv[0] = ModP((j * n * r * v - i * o * r * v - j * m * s * v + h * o * s * v + i * m * t * v -
                                h * n * t * v - j * n * q * w + i * o * q * w + j * l * s * w - g * o * s * w -
                                i * l * t * w +
                                g * n * t * w + j * m * q * x - h * o * q * x - j * l * r * x + g * o * r * x +
                                h * l * t * x -
                                g * m * t * x - i * m * q * y + h * n * q * y + i * l * r * y - g * n * r * y -
                                h * l * s * y +
                                g * m * s * y) * idet);
        var b1 = inv[1] = ModP((-e * n * r * v + d * o * r * v + e * m * s * v - c * o * s * v - d * m * t * v +
                                c * n * t * v + e * n * q * w - d * o * q * w - e * l * s * w + b * o * s * w +
                                d * l * t * w -
                                b * n * t * w - e * m * q * x + c * o * q * x + e * l * r * x - b * o * r * x -
                                c * l * t * x +
                                b * m * t * x + d * m * q * y - c * n * q * y - d * l * r * y + b * n * r * y +
                                c * l * s * y -
                                b * m * s * y) * idet);
        var c1 = inv[2] = ModP((e * i * r * v - d * j * r * v - e * h * s * v + c * j * s * v + d * h * t * v -
                                c * i * t * v - e * i * q * w + d * j * q * w + e * g * s * w - b * j * s * w -
                                d * g * t * w +
                                b * i * t * w + e * h * q * x - c * j * q * x - e * g * r * x + b * j * r * x +
                                c * g * t * x -
                                b * h * t * x - d * h * q * y + c * i * q * y + d * g * r * y - b * i * r * y -
                                c * g * s * y +
                                b * h * s * y) * idet);
        var d1 = inv[3] = ModP((-e * i * m * v + d * j * m * v + e * h * n * v - c * j * n * v - d * h * o * v +
                                c * i * o * v + e * i * l * w - d * j * l * w - e * g * n * w + b * j * n * w +
                                d * g * o * w -
                                b * i * o * w - e * h * l * x + c * j * l * x + e * g * m * x - b * j * m * x -
                                c * g * o * x +
                                b * h * o * x + d * h * l * y - c * i * l * y - d * g * m * y + b * i * m * y +
                                c * g * n * y -
                                b * h * n * y) * idet);
        var e1 = inv[4] = ModP((e * i * m * q - d * j * m * q - e * h * n * q + c * j * n * q + d * h * o * q -
                                c * i * o * q - e * i * l * r + d * j * l * r + e * g * n * r - b * j * n * r -
                                d * g * o * r +
                                b * i * o * r + e * h * l * s - c * j * l * s - e * g * m * s + b * j * m * s +
                                c * g * o * s -
                                b * h * o * s - d * h * l * t + c * i * l * t + d * g * m * t - b * i * m * t -
                                c * g * n * t +
                                b * h * n * t) * idet);
        var f1 = inv[5] = ModP((-j * n * r * u + i * o * r * u + j * m * s * u - h * o * s * u - i * m * t * u +
                                h * n * t * u + j * n * p * w - i * o * p * w - j * k * s * w + f * o * s * w +
                                i * k * t * w -
                                f * n * t * w - j * m * p * x + h * o * p * x + j * k * r * x - f * o * r * x -
                                h * k * t * x +
                                f * m * t * x + i * m * p * y - h * n * p * y - i * k * r * y + f * n * r * y +
                                h * k * s * y -
                                f * m * s * y) * idet);
        var g1 = inv[6] = ModP((e * n * r * u - d * o * r * u - e * m * s * u + c * o * s * u + d * m * t * u -
                                c * n * t * u - e * n * p * w + d * o * p * w + e * k * s * w - a * o * s * w -
                                d * k * t * w +
                                a * n * t * w + e * m * p * x - c * o * p * x - e * k * r * x + a * o * r * x +
                                c * k * t * x -
                                a * m * t * x - d * m * p * y + c * n * p * y + d * k * r * y - a * n * r * y -
                                c * k * s * y +
                                a * m * s * y) * idet);
        var h1 = inv[7] = ModP((-e * i * r * u + d * j * r * u + e * h * s * u - c * j * s * u - d * h * t * u +
                                c * i * t * u + e * i * p * w - d * j * p * w - e * f * s * w + a * j * s * w +
                                d * f * t * w -
                                a * i * t * w - e * h * p * x + c * j * p * x + e * f * r * x - a * j * r * x -
                                c * f * t * x +
                                a * h * t * x + d * h * p * y - c * i * p * y - d * f * r * y + a * i * r * y +
                                c * f * s * y -
                                a * h * s * y) * idet);
        var i1 = inv[8] = ModP((e * i * m * u - d * j * m * u - e * h * n * u + c * j * n * u + d * h * o * u -
                                c * i * o * u - e * i * k * w + d * j * k * w + e * f * n * w - a * j * n * w -
                                d * f * o * w +
                                a * i * o * w + e * h * k * x - c * j * k * x - e * f * m * x + a * j * m * x +
                                c * f * o * x -
                                a * h * o * x - d * h * k * y + c * i * k * y + d * f * m * y - a * i * m * y -
                                c * f * n * y +
                                a * h * n * y) * idet);
        var j1 = inv[9] = ModP((-e * i * m * p + d * j * m * p + e * h * n * p - c * j * n * p - d * h * o * p +
                                c * i * o * p + e * i * k * r - d * j * k * r - e * f * n * r + a * j * n * r +
                                d * f * o * r -
                                a * i * o * r - e * h * k * s + c * j * k * s + e * f * m * s - a * j * m * s -
                                c * f * o * s +
                                a * h * o * s + d * h * k * t - c * i * k * t - d * f * m * t + a * i * m * t +
                                c * f * n * t -
                                a * h * n * t) * idet);
        var k1 = inv[10] = ModP((j * n * q * u - i * o * q * u - j * l * s * u + g * o * s * u + i * l * t * u -
                                 g * n * t * u - j * n * p * v + i * o * p * v + j * k * s * v - f * o * s * v -
                                 i * k * t * v +
                                 f * n * t * v + j * l * p * x - g * o * p * x - j * k * q * x + f * o * q * x +
                                 g * k * t * x -
                                 f * l * t * x - i * l * p * y + g * n * p * y + i * k * q * y - f * n * q * y -
                                 g * k * s * y +
                                 f * l * s * y) * idet);
        var l1 = inv[11] = ModP((-e * n * q * u + d * o * q * u + e * l * s * u - b * o * s * u - d * l * t * u +
                                 b * n * t * u + e * n * p * v - d * o * p * v - e * k * s * v + a * o * s * v +
                                 d * k * t * v -
                                 a * n * t * v - e * l * p * x + b * o * p * x + e * k * q * x - a * o * q * x -
                                 b * k * t * x +
                                 a * l * t * x + d * l * p * y - b * n * p * y - d * k * q * y + a * n * q * y +
                                 b * k * s * y -
                                 a * l * s * y) * idet);
        var m1 = inv[12] = ModP((e * i * q * u - d * j * q * u - e * g * s * u + b * j * s * u + d * g * t * u -
                                 b * i * t * u - e * i * p * v + d * j * p * v + e * f * s * v - a * j * s * v -
                                 d * f * t * v +
                                 a * i * t * v + e * g * p * x - b * j * p * x - e * f * q * x + a * j * q * x +
                                 b * f * t * x -
                                 a * g * t * x - d * g * p * y + b * i * p * y + d * f * q * y - a * i * q * y -
                                 b * f * s * y +
                                 a * g * s * y) * idet);
        var n1 = inv[13] = ModP((-e * i * l * u + d * j * l * u + e * g * n * u - b * j * n * u - d * g * o * u +
                                 b * i * o * u + e * i * k * v - d * j * k * v - e * f * n * v + a * j * n * v +
                                 d * f * o * v -
                                 a * i * o * v - e * g * k * x + b * j * k * x + e * f * l * x - a * j * l * x -
                                 b * f * o * x +
                                 a * g * o * x + d * g * k * y - b * i * k * y - d * f * l * y + a * i * l * y +
                                 b * f * n * y -
                                 a * g * n * y) * idet);
        var o1 = inv[14] = ModP((e * i * l * p - d * j * l * p - e * g * n * p + b * j * n * p + d * g * o * p -
                                 b * i * o * p - e * i * k * q + d * j * k * q + e * f * n * q - a * j * n * q -
                                 d * f * o * q +
                                 a * i * o * q + e * g * k * s - b * j * k * s - e * f * l * s + a * j * l * s +
                                 b * f * o * s -
                                 a * g * o * s - d * g * k * t + b * i * k * t + d * f * l * t - a * i * l * t -
                                 b * f * n * t +
                                 a * g * n * t) * idet);
        var p1 = inv[15] = ModP((-j * m * q * u + h * o * q * u + j * l * r * u - g * o * r * u - h * l * t * u +
                                 g * m * t * u + j * m * p * v - h * o * p * v - j * k * r * v + f * o * r * v +
                                 h * k * t * v -
                                 f * m * t * v - j * l * p * w + g * o * p * w + j * k * q * w - f * o * q * w -
                                 g * k * t * w +
                                 f * l * t * w + h * l * p * y - g * m * p * y - h * k * q * y + f * m * q * y +
                                 g * k * r * y -
                                 f * l * r * y) * idet);
        var q1 = inv[16] = ModP((e * m * q * u - c * o * q * u - e * l * r * u + b * o * r * u + c * l * t * u -
                                 b * m * t * u - e * m * p * v + c * o * p * v + e * k * r * v - a * o * r * v -
                                 c * k * t * v +
                                 a * m * t * v + e * l * p * w - b * o * p * w - e * k * q * w + a * o * q * w +
                                 b * k * t * w -
                                 a * l * t * w - c * l * p * y + b * m * p * y + c * k * q * y - a * m * q * y -
                                 b * k * r * y +
                                 a * l * r * y) * idet);
        var r1 = inv[17] = ModP((-e * h * q * u + c * j * q * u + e * g * r * u - b * j * r * u - c * g * t * u +
                                 b * h * t * u + e * h * p * v - c * j * p * v - e * f * r * v + a * j * r * v +
                                 c * f * t * v -
                                 a * h * t * v - e * g * p * w + b * j * p * w + e * f * q * w - a * j * q * w -
                                 b * f * t * w +
                                 a * g * t * w + c * g * p * y - b * h * p * y - c * f * q * y + a * h * q * y +
                                 b * f * r * y -
                                 a * g * r * y) * idet);
        var s1 = inv[18] = ModP((e * h * l * u - c * j * l * u - e * g * m * u + b * j * m * u + c * g * o * u -
                                 b * h * o * u - e * h * k * v + c * j * k * v + e * f * m * v - a * j * m * v -
                                 c * f * o * v +
                                 a * h * o * v + e * g * k * w - b * j * k * w - e * f * l * w + a * j * l * w +
                                 b * f * o * w -
                                 a * g * o * w - c * g * k * y + b * h * k * y + c * f * l * y - a * h * l * y -
                                 b * f * m * y +
                                 a * g * m * y) * idet);
        var t1 = inv[19] = ModP((-e * h * l * p + c * j * l * p + e * g * m * p - b * j * m * p - c * g * o * p +
                                 b * h * o * p + e * h * k * q - c * j * k * q - e * f * m * q + a * j * m * q +
                                 c * f * o * q -
                                 a * h * o * q - e * g * k * r + b * j * k * r + e * f * l * r - a * j * l * r -
                                 b * f * o * r +
                                 a * g * o * r + c * g * k * t - b * h * k * t - c * f * l * t + a * h * l * t +
                                 b * f * m * t -
                                 a * g * m * t) * idet);
        var u1 = inv[20] = ModP((i * m * q * u - h * n * q * u - i * l * r * u + g * n * r * u + h * l * s * u -
                                 g * m * s * u - i * m * p * v + h * n * p * v + i * k * r * v - f * n * r * v -
                                 h * k * s * v +
                                 f * m * s * v + i * l * p * w - g * n * p * w - i * k * q * w + f * n * q * w +
                                 g * k * s * w -
                                 f * l * s * w - h * l * p * x + g * m * p * x + h * k * q * x - f * m * q * x -
                                 g * k * r * x +
                                 f * l * r * x) * idet);
        var v1 = inv[21] = ModP((-d * m * q * u + c * n * q * u + d * l * r * u - b * n * r * u - c * l * s * u +
                                 b * m * s * u + d * m * p * v - c * n * p * v - d * k * r * v + a * n * r * v +
                                 c * k * s * v -
                                 a * m * s * v - d * l * p * w + b * n * p * w + d * k * q * w - a * n * q * w -
                                 b * k * s * w +
                                 a * l * s * w + c * l * p * x - b * m * p * x - c * k * q * x + a * m * q * x +
                                 b * k * r * x -
                                 a * l * r * x) * idet);
        var w1 = inv[22] = ModP((d * h * q * u - c * i * q * u - d * g * r * u + b * i * r * u + c * g * s * u -
                                 b * h * s * u - d * h * p * v + c * i * p * v + d * f * r * v - a * i * r * v -
                                 c * f * s * v +
                                 a * h * s * v + d * g * p * w - b * i * p * w - d * f * q * w + a * i * q * w +
                                 b * f * s * w -
                                 a * g * s * w - c * g * p * x + b * h * p * x + c * f * q * x - a * h * q * x -
                                 b * f * r * x +
                                 a * g * r * x) * idet);
        var x1 = inv[23] = ModP((-d * h * l * u + c * i * l * u + d * g * m * u - b * i * m * u - c * g * n * u +
                                 b * h * n * u + d * h * k * v - c * i * k * v - d * f * m * v + a * i * m * v +
                                 c * f * n * v -
                                 a * h * n * v - d * g * k * w + b * i * k * w + d * f * l * w - a * i * l * w -
                                 b * f * n * w +
                                 a * g * n * w + c * g * k * x - b * h * k * x - c * f * l * x + a * h * l * x +
                                 b * f * m * x -
                                 a * g * m * x) * idet);
        var y1 = inv[24] = ModP((d * h * l * p - c * i * l * p - d * g * m * p + b * i * m * p + c * g * n * p -
                                 b * h * n * p - d * h * k * q + c * i * k * q + d * f * m * q - a * i * m * q -
                                 c * f * n * q +
                                 a * h * n * q + d * g * k * r - b * i * k * r - d * f * l * r + a * i * l * r +
                                 b * f * n * r -
                                 a * g * n * r - c * g * k * s + b * h * k * s + c * f * l * s - a * h * l * s -
                                 b * f * m * s +
                                 a * g * m * s) * idet);

        return a1 + P * (b1 + P * (c1 + P * (d1 + P * (e1 + P * (f1 + P * (g1 + P * (h1 + P * (i1 + P * (j1 + P * (k1 +
            P * (l1 + P * (m1 + P * (n1 + P * (o1 + P * (p1 + P * (q1 + P * (r1 + P * (s1 + P * (t1 + P * (u1 + P *
                (v1 + P * (w1 + P * (x1 + P * y1)))))))))))))))))))))));
    }
}