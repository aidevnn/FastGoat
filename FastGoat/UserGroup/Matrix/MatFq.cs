using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Matrix;

public struct MatFq : IElt<MatFq>
{
    public EPoly<ZnInt>[] Table { get; }
    public GLnq GLnq { get; }

    public MatFq(GLnq gl, EPoly<ZnInt>[] table)
    {
        GLnq = gl;
        var n = gl.Fq.F.Degree;
        var hash = table.SelectMany(a => a.Poly.CoefsExtended(n))
            .Aggregate(gl.Hash, (acc, h) => (acc, h.Hash).GetHashCode());
        Table = table.ToArray();
        Hash = hash;
    }

    public EPoly<ZnInt> Det => GLnq.Determinant(this);
    
    public bool IsOrder(int ord)
    {
        return Group.ElementIsOrder(GLnq, this, ord);
    }

    public MatFq T
    {
        get
        {
            var table = Table.ToArray();
            var n = GLnq.N;
            for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                table[j * n + i] = Table[i * n + j];

            return new(GLnq, table);
        }
    }

    public bool IsUT
    {
        get
        {
            var n = GLnq.N;
            var rg = n.Range();
            var t0 = Table;
            return rg.Grid2D(rg).Where(e => e.t1 > e.t2).All(e => t0[e.t1 * n + e.t2].IsZero());
        }
    }

    public bool IsLT
    {
        get
        {
            var n = GLnq.N;
            var rg = n.Range();
            var t0 = Table;
            return rg.Grid2D(rg).Where(e => e.t1 < e.t2).All(e => t0[e.t1 * n + e.t2].IsZero());
        }
    }

    public bool IsDiag
    {
        get
        {
            var n = GLnq.N;
            var rg = n.Range();
            var t0 = Table;
            return rg.Grid2D(rg).Where(e => e.t1 != e.t2).All(e => t0[e.t1 * n + e.t2].IsZero());
        }
    }

    public bool Is2ndDiag
    {
        get
        {
            var n = GLnq.N;
            var rg = n.Range();
            var t0 = Table;
            return rg.Grid2D(rg).Where(e => e.t1 + e.t2 != n - 1).All(e => t0[e.t1 * n + e.t2].IsZero());
        }
    }

    public bool IsSym
    {
        get
        {
            var n = GLnq.N;
            var rg = n.Range();
            var t0 = Table;
            return rg.Grid2D(rg).Where(e => e.t1 < e.t2).All(e => t0[e.t1 * n + e.t2].Equals(t0[e.t2 * n + e.t1]));
        }
    }

    public bool Equals(MatFq other) => Table.SequenceEqual(other.Table);

    public int CompareTo(MatFq other) => Table.SequenceCompareTo(other.Table);

    public int Hash { get; }
    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        return Table.ToKMatrix(GLnq.N).ToString();
    }
}