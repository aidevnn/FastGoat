using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.Matrix;

public class GLnq : IGroup<MatFq>
{
    public int N { get; }
    public Fq Fq { get; }
    public string Fmt { get; }

    public GLnq(int n, int q)
    {
        if (n < 1 || n > 4)
            throw new GroupException(GroupExceptionType.GroupDef);

        N = n;
        Fq = new(q);
        var zero = Fq.Zero;
        var one = Fq.One;
        var digits = Group.GenerateElements(Fq, Fq['x']).Max(e => $"{e}".Length);
        Fmt = $"{{0,{digits}}}";
        _neutralMat = n.Range().Grid2D().Select(a => a.t1 == a.t2 ? one : zero).ToArray();
        _hashNeutral = _neutralMat.Aggregate(0, (acc, a) => a.GetHashCode() + Fq.Q * acc);

        Hash = (n, q).GetHashCode();
        Name = $"GL({n}, {Fq})";
    }

    private EPoly<ZnInt>[] _neutralMat;
    private int _hashNeutral;
    public IEnumerator<MatFq> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();

    public bool Equals(IGroup<MatFq>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public MatFq this[params ValueType[] us]
    {
        get
        {
            var us0 = Enumerable.Range(0, N * N).Select(i => i < us.Length ? us[i] : 0).ToArray();
            if (us0.Any(u => u is EPoly<ZnInt>))
            {
                var table0 = us0.Select(u =>
                        u is EPoly<ZnInt> u0 ? u0 :
                        u is int u1 ? new EPoly<ZnInt>(Fq.F, Fq.F.One.Mul(u1)) :
                        throw new GroupException(GroupExceptionType.GroupDef))
                    .Select(e => e.Poly.Coefs.Select(ei => new ZnInt(Fq.P, ei.K)).ToArray())
                    .Select(e => new EPoly<ZnInt>(Fq.F, new KPoly<ZnInt>(Fq.F.x, ZnInt.ZnZero(Fq.P), e)))
                    .ToArray();

                var hash = table0.Aggregate(0, (acc, a) => a.GetHashCode() + Fq.Q * acc);
                return new(this, hash, table0);
            }
            else if (us0.All(u => u is int || u is char || u is ValueTuple<char, int>))
            {
                var table = us0.Select(i => Fq[i]).ToArray();
                var hash = table.Aggregate(0, (acc, a) => a.GetHashCode() + Fq.Q * acc);
                return new(this, hash, table);
            }

            throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

    public IEnumerable<MatFq> GetElements()
    {
        yield return Neutral();
    }

    public IEnumerable<MatFq> GetGenerators()
    {
        yield return Neutral();
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;

    public MatFq Create(EPoly<ZnInt>[] arr) => new(this, arr);
    
    public MatFq Neutral() => new(this, _hashNeutral, _neutralMat);

    public MatFq Invert(MatFq e)
    {
        var inv = Invert(e.Table);
        var hash = inv.Aggregate(0, (acc, a0) => a0.GetHashCode() + Fq.Q * acc);
        return new(this, hash, inv);
    }

    public MatFq Op(MatFq e1, MatFq e2)
    {
        var dot = MatrixProduct(N, e1.Table, e2.Table);
        // var dot = Dot(e1.Table, e2.Table); // TO DO 
        var hash = dot.Aggregate(0, (acc, a0) => a0.GetHashCode() + Fq.Q * acc);
        return new(this, hash, dot);
    }

    public EPoly<ZnInt> Determinant(MatFq mat)
    {
        if (N == 2)
            return Det2x2(mat.Table);
        if (N == 3)
            return Det3x3(mat.Table);
        if (N == 4)
            return Det4x4(mat.Table);
        if (N == 5)
            return Det5x5(mat.Table);

        throw new GroupException(GroupExceptionType.GroupDef);
    }

    EPoly<ZnInt>[] Invert(EPoly<ZnInt>[] mat)
    {
        if (N == 2)
            return Inv2x2(mat);
        if (N == 3)
            return Inv3x3(mat);
        if (N == 4)
            return Inv4x4(mat);
        if (N == 5)
            return Inv5x5(mat);

        throw new ArgumentException();
    }

    EPoly<ZnInt>[] Dot(EPoly<ZnInt>[] mat0, EPoly<ZnInt>[] mat1)
    {
        if (N == 2)
            return Dot2x2(mat0, mat1);
        if (N == 3)
            return Dot3x3(mat0, mat1);
        if (N == 4)
            return Dot4x4(mat0, mat1);
        if (N == 5)
            return Dot5x5(mat0, mat1);

        throw new ArgumentException();
    }

    EPoly<ZnInt>[] MatrixProduct(int aRows, EPoly<ZnInt>[] a, EPoly<ZnInt>[] b)
    {
        var aCols = a.Length / aRows;
        int bRows = aCols;
        var bCols = b.Length / bRows;
        if (a.Length != aRows * aCols || b.Length != bRows * bCols)
            throw new Exception();

        var c = new EPoly<ZnInt>[a.Length];
        var maxDegree = a.Max(e => e.Poly.Degree) + b.Max(e => e.Poly.Degree) + 1;
        var sum0 = new KPoly<ZnInt>(Fq.F.x, Fq.F.KZero, Enumerable.Repeat(Fq.F.KZero, maxDegree).ToArray());
        for (int i = 0; i < aRows; i++)
        {
            for (int j = 0; j < bCols; j++)
            {
                var sum = new KPoly<ZnInt>(Fq.F.x, Fq.F.KZero, Enumerable.Repeat(Fq.F.KZero, maxDegree).ToArray());
                for (int k = 0; k < aCols; k++)
                    sum.InPlaceAddProd(a[i * aCols + k].Poly, b[k * bCols + j].Poly);

                c[i * aRows + j] = new(Fq.F,
                    new KPoly<ZnInt>(Fq.F.x, Fq.F.KZero, sum.Coefs.TrimSeq().ToArray()).Div(Fq.F).rem);
            }
        }

        return c;
    }

    EPoly<ZnInt> Det2x2(EPoly<ZnInt>[] mat)
    {
        var (a, b) = (mat[0], mat[1]);
        var (c, d) = (mat[2], mat[3]);

        return (a * d - c * b);
    }

    EPoly<ZnInt>[] Inv2x2(EPoly<ZnInt>[] mat)
    {
        var (a, b) = (mat[0], mat[1]);
        var (c, d) = (mat[2], mat[3]);

        var det = a * d + -b * c;
        var idet = det.Inv();

        var inv = new EPoly<ZnInt>[4];
        inv[0] = (d) * idet;
        inv[1] = (-b) * idet;
        inv[2] = (-c) * idet;
        inv[3] = (a) * idet;

        return inv;
    }

    EPoly<ZnInt> Det3x3(EPoly<ZnInt>[] mat)
    {
        var (a, b, c) = (mat[0], mat[1], mat[2]);
        var (d, e, f) = (mat[3], mat[4], mat[5]);
        var (g, h, i) = (mat[6], mat[7], mat[8]);

        return a * e * i + -a * f * h + -b * d * i + b * f * g + c * d * h + -c * e * g;
    }

    EPoly<ZnInt>[] Inv3x3(EPoly<ZnInt>[] mat)
    {
        var (a, b, c) = (mat[0], mat[1], mat[2]);
        var (d, e, f) = (mat[3], mat[4], mat[5]);
        var (g, h, i) = (mat[6], mat[7], mat[8]);

        var det = a * e * i + -a * f * h + -b * d * i + b * f * g + c * d * h + -c * e * g;
        var idet = det.Inv();

        var inv = new EPoly<ZnInt>[9];
        inv[0] = (e * i + -f * h) * idet;
        inv[1] = (-b * i + c * h) * idet;
        inv[2] = (b * f + -c * e) * idet;
        inv[3] = (-d * i + f * g) * idet;
        inv[4] = (a * i + -c * g) * idet;
        inv[5] = (-a * f + c * d) * idet;
        inv[6] = (d * h + -e * g) * idet;
        inv[7] = (-a * h + b * g) * idet;
        inv[8] = (a * e + -b * d) * idet;

        return inv;
    }

    EPoly<ZnInt> Det4x4(EPoly<ZnInt>[] mat)
    {
        var (a, b, c, d) = (mat[0], mat[1], mat[2], mat[3]);
        var (e, f, g, h) = (mat[4], mat[5], mat[6], mat[7]);
        var (i, j, k, l) = (mat[8], mat[9], mat[10], mat[11]);
        var (m, n, o, p) = (mat[12], mat[13], mat[14], mat[15]);

        return a * f * k * p + -a * f * l * o + -a * g * j * p + a * g * l * n + a * h * j * o + -a * h * k * n +
               -b * e * k * p + b * e * l * o + b * g * i * p + -b * g * l * m + -b * h * i * o + b * h * k * m +
               c * e * j * p + -c * e * l * n + -c * f * i * p + c * f * l * m + c * h * i * n + -c * h * j * m +
               -d * e * j * o + d * e * k * n + d * f * i * o + -d * f * k * m + -d * g * i * n + d * g * j * m;
    }

    EPoly<ZnInt>[] Inv4x4(EPoly<ZnInt>[] mat)
    {
        var (a, b, c, d) = (mat[0], mat[1], mat[2], mat[3]);
        var (e, f, g, h) = (mat[4], mat[5], mat[6], mat[7]);
        var (i, j, k, l) = (mat[8], mat[9], mat[10], mat[11]);
        var (m, n, o, p) = (mat[12], mat[13], mat[14], mat[15]);

        var det = a * f * k * p + -a * f * l * o + -a * g * j * p + a * g * l * n + a * h * j * o + -a * h * k * n +
                  -b * e * k * p + b * e * l * o + b * g * i * p + -b * g * l * m + -b * h * i * o + b * h * k * m +
                  c * e * j * p + -c * e * l * n + -c * f * i * p + c * f * l * m + c * h * i * n + -c * h * j * m +
                  -d * e * j * o + d * e * k * n + d * f * i * o + -d * f * k * m + -d * g * i * n + d * g * j * m;
        var idet = det.Inv();

        var inv = new EPoly<ZnInt>[16];
        inv[0] = (f * k * p + -f * l * o + -g * j * p + g * l * n + h * j * o + -h * k * n) * idet;
        inv[1] = (-b * k * p + b * l * o + c * j * p + -c * l * n + -d * j * o + d * k * n) * idet;
        inv[2] = (b * g * p + -b * h * o + -c * f * p + c * h * n + d * f * o + -d * g * n) * idet;
        inv[3] = (-b * g * l + b * h * k + c * f * l + -c * h * j + -d * f * k + d * g * j) * idet;
        inv[4] = (-e * k * p + e * l * o + g * i * p + -g * l * m + -h * i * o + h * k * m) * idet;
        inv[5] = (a * k * p + -a * l * o + -c * i * p + c * l * m + d * i * o + -d * k * m) * idet;
        inv[6] = (-a * g * p + a * h * o + c * e * p + -c * h * m + -d * e * o + d * g * m) * idet;
        inv[7] = (a * g * l + -a * h * k + -c * e * l + c * h * i + d * e * k + -d * g * i) * idet;
        inv[8] = (e * j * p + -e * l * n + -f * i * p + f * l * m + h * i * n + -h * j * m) * idet;
        inv[9] = (-a * j * p + a * l * n + b * i * p + -b * l * m + -d * i * n + d * j * m) * idet;
        inv[10] = (a * f * p + -a * h * n + -b * e * p + b * h * m + d * e * n + -d * f * m) * idet;
        inv[11] = (-a * f * l + a * h * j + b * e * l + -b * h * i + -d * e * j + d * f * i) * idet;
        inv[12] = (-e * j * o + e * k * n + f * i * o + -f * k * m + -g * i * n + g * j * m) * idet;
        inv[13] = (a * j * o + -a * k * n + -b * i * o + b * k * m + c * i * n + -c * j * m) * idet;
        inv[14] = (-a * f * o + a * g * n + b * e * o + -b * g * m + -c * e * n + c * f * m) * idet;
        inv[15] = (a * f * k + -a * g * j + -b * e * k + b * g * i + c * e * j + -c * f * i) * idet;

        return inv;
    }

    EPoly<ZnInt> Det5x5(EPoly<ZnInt>[] mat)
    {
        var (a, b, c, d, e) = (mat[0], mat[1], mat[2], mat[3], mat[4]);
        var (f, g, h, i, j) = (mat[5], mat[6], mat[7], mat[8], mat[9]);
        var (k, l, m, n, o) = (mat[10], mat[11], mat[12], mat[13], mat[14]);
        var (p, q, r, s, t) = (mat[15], mat[16], mat[17], mat[18], mat[19]);
        var (v, w, x, y, z) = (mat[20], mat[21], mat[22], mat[23], mat[24]);

        return a * g * m * s * z + -a * g * m * t * y + -a * g * n * r * z + a * g * n * t * x + a * g * o * r * y +
               -a * g * o * s * x + -a * h * l * s * z + a * h * l * t * y + a * h * n * q * z + -a * h * n * t * w +
               -a * h * o * q * y + a * h * o * s * w + a * i * l * r * z + -a * i * l * t * x + -a * i * m * q * z +
               a * i * m * t * w + a * i * o * q * x + -a * i * o * r * w + -a * j * l * r * y + a * j * l * s * x +
               a * j * m * q * y + -a * j * m * s * w + -a * j * n * q * x + a * j * n * r * w + -b * f * m * s * z +
               b * f * m * t * y + b * f * n * r * z + -b * f * n * t * x + -b * f * o * r * y + b * f * o * s * x +
               b * h * k * s * z + -b * h * k * t * y + -b * h * n * p * z + b * h * n * t * v + b * h * o * p * y +
               -b * h * o * s * v + -b * i * k * r * z + b * i * k * t * x + b * i * m * p * z + -b * i * m * t * v +
               -b * i * o * p * x + b * i * o * r * v + b * j * k * r * y + -b * j * k * s * x + -b * j * m * p * y +
               b * j * m * s * v + b * j * n * p * x + -b * j * n * r * v + c * f * l * s * z + -c * f * l * t * y +
               -c * f * n * q * z + c * f * n * t * w + c * f * o * q * y + -c * f * o * s * w + -c * g * k * s * z +
               c * g * k * t * y + c * g * n * p * z + -c * g * n * t * v + -c * g * o * p * y + c * g * o * s * v +
               c * i * k * q * z + -c * i * k * t * w + -c * i * l * p * z + c * i * l * t * v + c * i * o * p * w +
               -c * i * o * q * v + -c * j * k * q * y + c * j * k * s * w + c * j * l * p * y + -c * j * l * s * v +
               -c * j * n * p * w + c * j * n * q * v + -d * f * l * r * z + d * f * l * t * x + d * f * m * q * z +
               -d * f * m * t * w + -d * f * o * q * x + d * f * o * r * w + d * g * k * r * z + -d * g * k * t * x +
               -d * g * m * p * z + d * g * m * t * v + d * g * o * p * x + -d * g * o * r * v + -d * h * k * q * z +
               d * h * k * t * w + d * h * l * p * z + -d * h * l * t * v + -d * h * o * p * w + d * h * o * q * v +
               d * j * k * q * x + -d * j * k * r * w + -d * j * l * p * x + d * j * l * r * v + d * j * m * p * w +
               -d * j * m * q * v + e * f * l * r * y + -e * f * l * s * x + -e * f * m * q * y + e * f * m * s * w +
               e * f * n * q * x + -e * f * n * r * w + -e * g * k * r * y + e * g * k * s * x + e * g * m * p * y +
               -e * g * m * s * v + -e * g * n * p * x + e * g * n * r * v + e * h * k * q * y + -e * h * k * s * w +
               -e * h * l * p * y + e * h * l * s * v + e * h * n * p * w + -e * h * n * q * v + -e * i * k * q * x +
               e * i * k * r * w + e * i * l * p * x + -e * i * l * r * v + -e * i * m * p * w + e * i * m * q * v;
    }

    EPoly<ZnInt>[] Inv5x5(EPoly<ZnInt>[] mat)
    {
        var (a, b, c, d, e) = (mat[0], mat[1], mat[2], mat[3], mat[4]);
        var (f, g, h, i, j) = (mat[5], mat[6], mat[7], mat[8], mat[9]);
        var (k, l, m, n, o) = (mat[10], mat[11], mat[12], mat[13], mat[14]);
        var (p, q, r, s, t) = (mat[15], mat[16], mat[17], mat[18], mat[19]);
        var (v, w, x, y, z) = (mat[20], mat[21], mat[22], mat[23], mat[24]);

        var det = a * g * m * s * z + -a * g * m * t * y + -a * g * n * r * z + a * g * n * t * x + a * g * o * r * y +
                  -a * g * o * s * x + -a * h * l * s * z + a * h * l * t * y + a * h * n * q * z + -a * h * n * t * w +
                  -a * h * o * q * y + a * h * o * s * w + a * i * l * r * z + -a * i * l * t * x + -a * i * m * q * z +
                  a * i * m * t * w + a * i * o * q * x + -a * i * o * r * w + -a * j * l * r * y + a * j * l * s * x +
                  a * j * m * q * y + -a * j * m * s * w + -a * j * n * q * x + a * j * n * r * w + -b * f * m * s * z +
                  b * f * m * t * y + b * f * n * r * z + -b * f * n * t * x + -b * f * o * r * y + b * f * o * s * x +
                  b * h * k * s * z + -b * h * k * t * y + -b * h * n * p * z + b * h * n * t * v + b * h * o * p * y +
                  -b * h * o * s * v + -b * i * k * r * z + b * i * k * t * x + b * i * m * p * z + -b * i * m * t * v +
                  -b * i * o * p * x + b * i * o * r * v + b * j * k * r * y + -b * j * k * s * x + -b * j * m * p * y +
                  b * j * m * s * v + b * j * n * p * x + -b * j * n * r * v + c * f * l * s * z + -c * f * l * t * y +
                  -c * f * n * q * z + c * f * n * t * w + c * f * o * q * y + -c * f * o * s * w + -c * g * k * s * z +
                  c * g * k * t * y + c * g * n * p * z + -c * g * n * t * v + -c * g * o * p * y + c * g * o * s * v +
                  c * i * k * q * z + -c * i * k * t * w + -c * i * l * p * z + c * i * l * t * v + c * i * o * p * w +
                  -c * i * o * q * v + -c * j * k * q * y + c * j * k * s * w + c * j * l * p * y + -c * j * l * s * v +
                  -c * j * n * p * w + c * j * n * q * v + -d * f * l * r * z + d * f * l * t * x + d * f * m * q * z +
                  -d * f * m * t * w + -d * f * o * q * x + d * f * o * r * w + d * g * k * r * z + -d * g * k * t * x +
                  -d * g * m * p * z + d * g * m * t * v + d * g * o * p * x + -d * g * o * r * v + -d * h * k * q * z +
                  d * h * k * t * w + d * h * l * p * z + -d * h * l * t * v + -d * h * o * p * w + d * h * o * q * v +
                  d * j * k * q * x + -d * j * k * r * w + -d * j * l * p * x + d * j * l * r * v + d * j * m * p * w +
                  -d * j * m * q * v + e * f * l * r * y + -e * f * l * s * x + -e * f * m * q * y + e * f * m * s * w +
                  e * f * n * q * x + -e * f * n * r * w + -e * g * k * r * y + e * g * k * s * x + e * g * m * p * y +
                  -e * g * m * s * v + -e * g * n * p * x + e * g * n * r * v + e * h * k * q * y + -e * h * k * s * w +
                  -e * h * l * p * y + e * h * l * s * v + e * h * n * p * w + -e * h * n * q * v + -e * i * k * q * x +
                  e * i * k * r * w + e * i * l * p * x + -e * i * l * r * v + -e * i * m * p * w + e * i * m * q * v;
        var idet = det.Inv();

        var inv = new EPoly<ZnInt>[25];
        inv[0] = (g * m * s * z + -g * m * t * y + -g * n * r * z + g * n * t * x + g * o * r * y + -g * o * s * x +
                  -h * l * s * z + h * l * t * y + h * n * q * z + -h * n * t * w + -h * o * q * y + h * o * s * w +
                  i * l * r * z + -i * l * t * x + -i * m * q * z + i * m * t * w + i * o * q * x + -i * o * r * w +
                  -j * l * r * y + j * l * s * x + j * m * q * y + -j * m * s * w + -j * n * q * x + j * n * r * w) *
                 idet;
        inv[1] = (-b * m * s * z + b * m * t * y + b * n * r * z + -b * n * t * x + -b * o * r * y + b * o * s * x +
                  c * l * s * z + -c * l * t * y + -c * n * q * z + c * n * t * w + c * o * q * y + -c * o * s * w +
                  -d * l * r * z + d * l * t * x + d * m * q * z + -d * m * t * w + -d * o * q * x + d * o * r * w +
                  e * l * r * y + -e * l * s * x + -e * m * q * y + e * m * s * w + e * n * q * x + -e * n * r * w) *
                 idet;
        inv[2] = (b * h * s * z + -b * h * t * y + -b * i * r * z + b * i * t * x + b * j * r * y + -b * j * s * x +
                  -c * g * s * z + c * g * t * y + c * i * q * z + -c * i * t * w + -c * j * q * y + c * j * s * w +
                  d * g * r * z + -d * g * t * x + -d * h * q * z + d * h * t * w + d * j * q * x + -d * j * r * w +
                  -e * g * r * y + e * g * s * x + e * h * q * y + -e * h * s * w + -e * i * q * x + e * i * r * w) *
                 idet;
        inv[3] = (-b * h * n * z + b * h * o * y + b * i * m * z + -b * i * o * x + -b * j * m * y + b * j * n * x +
                  c * g * n * z + -c * g * o * y + -c * i * l * z + c * i * o * w + c * j * l * y + -c * j * n * w +
                  -d * g * m * z + d * g * o * x + d * h * l * z + -d * h * o * w + -d * j * l * x + d * j * m * w +
                  e * g * m * y + -e * g * n * x + -e * h * l * y + e * h * n * w + e * i * l * x + -e * i * m * w) *
                 idet;
        inv[4] = (b * h * n * t + -b * h * o * s + -b * i * m * t + b * i * o * r + b * j * m * s + -b * j * n * r +
                  -c * g * n * t + c * g * o * s + c * i * l * t + -c * i * o * q + -c * j * l * s + c * j * n * q +
                  d * g * m * t + -d * g * o * r + -d * h * l * t + d * h * o * q + d * j * l * r + -d * j * m * q +
                  -e * g * m * s + e * g * n * r + e * h * l * s + -e * h * n * q + -e * i * l * r + e * i * m * q) *
                 idet;
        inv[5] = (-f * m * s * z + f * m * t * y + f * n * r * z + -f * n * t * x + -f * o * r * y + f * o * s * x +
                  h * k * s * z + -h * k * t * y + -h * n * p * z + h * n * t * v + h * o * p * y + -h * o * s * v +
                  -i * k * r * z + i * k * t * x + i * m * p * z + -i * m * t * v + -i * o * p * x + i * o * r * v +
                  j * k * r * y + -j * k * s * x + -j * m * p * y + j * m * s * v + j * n * p * x + -j * n * r * v) *
                 idet;
        inv[6] = (a * m * s * z + -a * m * t * y + -a * n * r * z + a * n * t * x + a * o * r * y + -a * o * s * x +
                  -c * k * s * z + c * k * t * y + c * n * p * z + -c * n * t * v + -c * o * p * y + c * o * s * v +
                  d * k * r * z + -d * k * t * x + -d * m * p * z + d * m * t * v + d * o * p * x + -d * o * r * v +
                  -e * k * r * y + e * k * s * x + e * m * p * y + -e * m * s * v + -e * n * p * x + e * n * r * v) *
                 idet;
        inv[7] = (-a * h * s * z + a * h * t * y + a * i * r * z + -a * i * t * x + -a * j * r * y + a * j * s * x +
                  c * f * s * z + -c * f * t * y + -c * i * p * z + c * i * t * v + c * j * p * y + -c * j * s * v +
                  -d * f * r * z + d * f * t * x + d * h * p * z + -d * h * t * v + -d * j * p * x + d * j * r * v +
                  e * f * r * y + -e * f * s * x + -e * h * p * y + e * h * s * v + e * i * p * x + -e * i * r * v) *
                 idet;
        inv[8] = (a * h * n * z + -a * h * o * y + -a * i * m * z + a * i * o * x + a * j * m * y + -a * j * n * x +
                  -c * f * n * z + c * f * o * y + c * i * k * z + -c * i * o * v + -c * j * k * y + c * j * n * v +
                  d * f * m * z + -d * f * o * x + -d * h * k * z + d * h * o * v + d * j * k * x + -d * j * m * v +
                  -e * f * m * y + e * f * n * x + e * h * k * y + -e * h * n * v + -e * i * k * x + e * i * m * v) *
                 idet;
        inv[9] = (-a * h * n * t + a * h * o * s + a * i * m * t + -a * i * o * r + -a * j * m * s + a * j * n * r +
                  c * f * n * t + -c * f * o * s + -c * i * k * t + c * i * o * p + c * j * k * s + -c * j * n * p +
                  -d * f * m * t + d * f * o * r + d * h * k * t + -d * h * o * p + -d * j * k * r + d * j * m * p +
                  e * f * m * s + -e * f * n * r + -e * h * k * s + e * h * n * p + e * i * k * r + -e * i * m * p) *
                 idet;
        inv[10] = (f * l * s * z + -f * l * t * y + -f * n * q * z + f * n * t * w + f * o * q * y + -f * o * s * w +
                   -g * k * s * z + g * k * t * y + g * n * p * z + -g * n * t * v + -g * o * p * y + g * o * s * v +
                   i * k * q * z + -i * k * t * w + -i * l * p * z + i * l * t * v + i * o * p * w + -i * o * q * v +
                   -j * k * q * y + j * k * s * w + j * l * p * y + -j * l * s * v + -j * n * p * w + j * n * q * v) *
                  idet;
        inv[11] = (-a * l * s * z + a * l * t * y + a * n * q * z + -a * n * t * w + -a * o * q * y + a * o * s * w +
                   b * k * s * z + -b * k * t * y + -b * n * p * z + b * n * t * v + b * o * p * y + -b * o * s * v +
                   -d * k * q * z + d * k * t * w + d * l * p * z + -d * l * t * v + -d * o * p * w + d * o * q * v +
                   e * k * q * y + -e * k * s * w + -e * l * p * y + e * l * s * v + e * n * p * w + -e * n * q * v) *
                  idet;
        inv[12] = (a * g * s * z + -a * g * t * y + -a * i * q * z + a * i * t * w + a * j * q * y + -a * j * s * w +
                   -b * f * s * z + b * f * t * y + b * i * p * z + -b * i * t * v + -b * j * p * y + b * j * s * v +
                   d * f * q * z + -d * f * t * w + -d * g * p * z + d * g * t * v + d * j * p * w + -d * j * q * v +
                   -e * f * q * y + e * f * s * w + e * g * p * y + -e * g * s * v + -e * i * p * w + e * i * q * v) *
                  idet;
        inv[13] = (-a * g * n * z + a * g * o * y + a * i * l * z + -a * i * o * w + -a * j * l * y + a * j * n * w +
                   b * f * n * z + -b * f * o * y + -b * i * k * z + b * i * o * v + b * j * k * y + -b * j * n * v +
                   -d * f * l * z + d * f * o * w + d * g * k * z + -d * g * o * v + -d * j * k * w + d * j * l * v +
                   e * f * l * y + -e * f * n * w + -e * g * k * y + e * g * n * v + e * i * k * w + -e * i * l * v) *
                  idet;
        inv[14] = (a * g * n * t + -a * g * o * s + -a * i * l * t + a * i * o * q + a * j * l * s + -a * j * n * q +
                   -b * f * n * t + b * f * o * s + b * i * k * t + -b * i * o * p + -b * j * k * s + b * j * n * p +
                   d * f * l * t + -d * f * o * q + -d * g * k * t + d * g * o * p + d * j * k * q + -d * j * l * p +
                   -e * f * l * s + e * f * n * q + e * g * k * s + -e * g * n * p + -e * i * k * q + e * i * l * p) *
                  idet;
        inv[15] = (-f * l * r * z + f * l * t * x + f * m * q * z + -f * m * t * w + -f * o * q * x + f * o * r * w +
                   g * k * r * z + -g * k * t * x + -g * m * p * z + g * m * t * v + g * o * p * x + -g * o * r * v +
                   -h * k * q * z + h * k * t * w + h * l * p * z + -h * l * t * v + -h * o * p * w + h * o * q * v +
                   j * k * q * x + -j * k * r * w + -j * l * p * x + j * l * r * v + j * m * p * w + -j * m * q * v) *
                  idet;
        inv[16] = (a * l * r * z + -a * l * t * x + -a * m * q * z + a * m * t * w + a * o * q * x + -a * o * r * w +
                   -b * k * r * z + b * k * t * x + b * m * p * z + -b * m * t * v + -b * o * p * x + b * o * r * v +
                   c * k * q * z + -c * k * t * w + -c * l * p * z + c * l * t * v + c * o * p * w + -c * o * q * v +
                   -e * k * q * x + e * k * r * w + e * l * p * x + -e * l * r * v + -e * m * p * w + e * m * q * v) *
                  idet;
        inv[17] = (-a * g * r * z + a * g * t * x + a * h * q * z + -a * h * t * w + -a * j * q * x + a * j * r * w +
                   b * f * r * z + -b * f * t * x + -b * h * p * z + b * h * t * v + b * j * p * x + -b * j * r * v +
                   -c * f * q * z + c * f * t * w + c * g * p * z + -c * g * t * v + -c * j * p * w + c * j * q * v +
                   e * f * q * x + -e * f * r * w + -e * g * p * x + e * g * r * v + e * h * p * w + -e * h * q * v) *
                  idet;
        inv[18] = (a * g * m * z + -a * g * o * x + -a * h * l * z + a * h * o * w + a * j * l * x + -a * j * m * w +
                   -b * f * m * z + b * f * o * x + b * h * k * z + -b * h * o * v + -b * j * k * x + b * j * m * v +
                   c * f * l * z + -c * f * o * w + -c * g * k * z + c * g * o * v + c * j * k * w + -c * j * l * v +
                   -e * f * l * x + e * f * m * w + e * g * k * x + -e * g * m * v + -e * h * k * w + e * h * l * v) *
                  idet;
        inv[19] = (-a * g * m * t + a * g * o * r + a * h * l * t + -a * h * o * q + -a * j * l * r + a * j * m * q +
                   b * f * m * t + -b * f * o * r + -b * h * k * t + b * h * o * p + b * j * k * r + -b * j * m * p +
                   -c * f * l * t + c * f * o * q + c * g * k * t + -c * g * o * p + -c * j * k * q + c * j * l * p +
                   e * f * l * r + -e * f * m * q + -e * g * k * r + e * g * m * p + e * h * k * q + -e * h * l * p) *
                  idet;
        inv[20] = (f * l * r * y + -f * l * s * x + -f * m * q * y + f * m * s * w + f * n * q * x + -f * n * r * w +
                   -g * k * r * y + g * k * s * x + g * m * p * y + -g * m * s * v + -g * n * p * x + g * n * r * v +
                   h * k * q * y + -h * k * s * w + -h * l * p * y + h * l * s * v + h * n * p * w + -h * n * q * v +
                   -i * k * q * x + i * k * r * w + i * l * p * x + -i * l * r * v + -i * m * p * w + i * m * q * v) *
                  idet;
        inv[21] = (-a * l * r * y + a * l * s * x + a * m * q * y + -a * m * s * w + -a * n * q * x + a * n * r * w +
                   b * k * r * y + -b * k * s * x + -b * m * p * y + b * m * s * v + b * n * p * x + -b * n * r * v +
                   -c * k * q * y + c * k * s * w + c * l * p * y + -c * l * s * v + -c * n * p * w + c * n * q * v +
                   d * k * q * x + -d * k * r * w + -d * l * p * x + d * l * r * v + d * m * p * w + -d * m * q * v) *
                  idet;
        inv[22] = (a * g * r * y + -a * g * s * x + -a * h * q * y + a * h * s * w + a * i * q * x + -a * i * r * w +
                   -b * f * r * y + b * f * s * x + b * h * p * y + -b * h * s * v + -b * i * p * x + b * i * r * v +
                   c * f * q * y + -c * f * s * w + -c * g * p * y + c * g * s * v + c * i * p * w + -c * i * q * v +
                   -d * f * q * x + d * f * r * w + d * g * p * x + -d * g * r * v + -d * h * p * w + d * h * q * v) *
                  idet;
        inv[23] = (-a * g * m * y + a * g * n * x + a * h * l * y + -a * h * n * w + -a * i * l * x + a * i * m * w +
                   b * f * m * y + -b * f * n * x + -b * h * k * y + b * h * n * v + b * i * k * x + -b * i * m * v +
                   -c * f * l * y + c * f * n * w + c * g * k * y + -c * g * n * v + -c * i * k * w + c * i * l * v +
                   d * f * l * x + -d * f * m * w + -d * g * k * x + d * g * m * v + d * h * k * w + -d * h * l * v) *
                  idet;
        inv[24] = (a * g * m * s + -a * g * n * r + -a * h * l * s + a * h * n * q + a * i * l * r + -a * i * m * q +
                   -b * f * m * s + b * f * n * r + b * h * k * s + -b * h * n * p + -b * i * k * r + b * i * m * p +
                   c * f * l * s + -c * f * n * q + -c * g * k * s + c * g * n * p + c * i * k * q + -c * i * l * p +
                   -d * f * l * r + d * f * m * q + d * g * k * r + -d * g * m * p + -d * h * k * q + d * h * l * p) *
                  idet;

        return inv;
    }

    public EPoly<ZnInt>[] Dot2x2(EPoly<ZnInt>[] mat0, EPoly<ZnInt>[] mat1)
    {
        var (a, b) = (mat0[0], mat0[1]);
        var (c, d) = (mat0[2], mat0[3]);

        var (A, B) = (mat1[0], mat1[1]);
        var (C, D) = (mat1[2], mat1[3]);

        var mat = new EPoly<ZnInt>[4];
        mat[0] = A * a + C * b;
        mat[1] = B * a + D * b;
        mat[2] = A * c + C * d;
        mat[3] = B * c + D * d;

        return mat;
    }

    EPoly<ZnInt>[] Dot3x3(EPoly<ZnInt>[] mat0, EPoly<ZnInt>[] mat1)
    {
        var (a, b, c) = (mat0[0], mat0[1], mat0[2]);
        var (d, e, f) = (mat0[3], mat0[4], mat0[5]);
        var (g, h, i) = (mat0[6], mat0[7], mat0[8]);

        var (A, B, C) = (mat1[0], mat1[1], mat1[2]);
        var (D, E, F) = (mat1[3], mat1[4], mat1[5]);
        var (G, H, I) = (mat1[6], mat1[7], mat1[8]);

        var mat = new EPoly<ZnInt>[9];
        mat[0] = A * a + D * b + G * c;
        mat[1] = B * a + E * b + H * c;
        mat[2] = C * a + F * b + I * c;
        mat[3] = A * d + D * e + G * f;
        mat[4] = B * d + E * e + H * f;
        mat[5] = C * d + F * e + I * f;
        mat[6] = A * g + D * h + G * i;
        mat[7] = B * g + E * h + H * i;
        mat[8] = C * g + F * h + I * i;

        return mat;
    }

    EPoly<ZnInt>[] Dot4x4(EPoly<ZnInt>[] mat0, EPoly<ZnInt>[] mat1)
    {
        var (a, b, c, d) = (mat0[0], mat0[1], mat0[2], mat0[3]);
        var (e, f, g, h) = (mat0[4], mat0[5], mat0[6], mat0[7]);
        var (i, j, k, l) = (mat0[8], mat0[9], mat0[10], mat0[11]);
        var (m, n, o, p) = (mat0[12], mat0[13], mat0[14], mat0[15]);

        var (A, B, C, D) = (mat1[0], mat1[1], mat1[2], mat1[3]);
        var (E, F, G, H) = (mat1[4], mat1[5], mat1[6], mat1[7]);
        var (I, J, K, L) = (mat1[8], mat1[9], mat1[10], mat1[11]);
        var (M, N, O, P) = (mat1[12], mat1[13], mat1[14], mat1[15]);

        var mat = new EPoly<ZnInt>[16];
        mat[0] = A * a + E * b + I * c + M * d;
        mat[1] = B * a + F * b + J * c + N * d;
        mat[2] = C * a + G * b + K * c + O * d;
        mat[3] = D * a + H * b + L * c + P * d;
        mat[4] = A * e + E * f + I * g + M * h;
        mat[5] = B * e + F * f + J * g + N * h;
        mat[6] = C * e + G * f + K * g + O * h;
        mat[7] = D * e + H * f + L * g + P * h;
        mat[8] = A * i + E * j + I * k + M * l;
        mat[9] = B * i + F * j + J * k + N * l;
        mat[10] = C * i + G * j + K * k + O * l;
        mat[11] = D * i + H * j + L * k + P * l;
        mat[12] = A * m + E * n + I * o + M * p;
        mat[13] = B * m + F * n + J * o + N * p;
        mat[14] = C * m + G * n + K * o + O * p;
        mat[15] = D * m + H * n + L * o + P * p;

        return mat;
    }

    EPoly<ZnInt>[] Dot5x5(EPoly<ZnInt>[] mat0, EPoly<ZnInt>[] mat1)
    {
        var (a, b, c, d, e) = (mat0[0], mat0[1], mat0[2], mat0[3], mat0[4]);
        var (f, g, h, i, j) = (mat0[5], mat0[6], mat0[7], mat0[8], mat0[9]);
        var (k, l, m, n, o) = (mat0[10], mat0[11], mat0[12], mat0[13], mat0[14]);
        var (p, q, r, s, t) = (mat0[15], mat0[16], mat0[17], mat0[18], mat0[19]);
        var (u, v, w, x, y) = (mat0[20], mat0[21], mat0[22], mat0[23], mat0[24]);

        var (A, B, C, D, E) = (mat1[0], mat1[1], mat1[2], mat1[3], mat1[4]);
        var (F, G, H, I, J) = (mat1[5], mat1[6], mat1[7], mat1[8], mat1[9]);
        var (K, L, M, N, O) = (mat1[10], mat1[11], mat1[12], mat1[13], mat1[14]);
        var (P, Q, R, S, T) = (mat1[15], mat1[16], mat1[17], mat1[18], mat1[19]);
        var (U, V, W, X, Y) = (mat1[20], mat1[21], mat1[22], mat1[23], mat1[24]);

        var mat = new EPoly<ZnInt>[25];
        mat[0] = A * a + F * b + K * c + P * d + U * e;
        mat[1] = B * a + G * b + L * c + Q * d + V * e;
        mat[2] = C * a + H * b + M * c + R * d + W * e;
        mat[3] = D * a + I * b + N * c + S * d + X * e;
        mat[4] = E * a + J * b + O * c + T * d + Y * e;
        mat[5] = A * f + F * g + K * h + P * i + U * j;
        mat[6] = B * f + G * g + L * h + Q * i + V * j;
        mat[7] = C * f + H * g + M * h + R * i + W * j;
        mat[8] = D * f + I * g + N * h + S * i + X * j;
        mat[9] = E * f + J * g + O * h + T * i + Y * j;
        mat[10] = A * k + F * l + K * m + P * n + U * o;
        mat[11] = B * k + G * l + L * m + Q * n + V * o;
        mat[12] = C * k + H * l + M * m + R * n + W * o;
        mat[13] = D * k + I * l + N * m + S * n + X * o;
        mat[14] = E * k + J * l + O * m + T * n + Y * o;
        mat[15] = A * p + F * q + K * r + P * s + U * t;
        mat[16] = B * p + G * q + L * r + Q * s + V * t;
        mat[17] = C * p + H * q + M * r + R * s + W * t;
        mat[18] = D * p + I * q + N * r + S * s + X * t;
        mat[19] = E * p + J * q + O * r + T * s + Y * t;
        mat[20] = A * u + F * v + K * w + P * x + U * y;
        mat[21] = B * u + G * v + L * w + Q * x + V * y;
        mat[22] = C * u + H * v + M * w + R * x + W * y;
        mat[23] = D * u + I * v + N * w + S * x + X * y;
        mat[24] = E * u + J * v + O * w + T * x + Y * y;

        return mat;
    }
}