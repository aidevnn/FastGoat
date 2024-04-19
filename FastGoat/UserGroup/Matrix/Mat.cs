using System.Text;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Matrix;

public struct Mat : IElt<Mat>
{
    public int[] Table { get; }
    public GL GL { get; }

    public Mat(GL gl)
    {
        GL = gl;
        Table = GL.TableNeutral.ToArray();
        Hash = GL.HashNeutral;
    }

    public Mat(GL gl, int hash, int[] table)
    {
        GL = gl;
        Table = table.ToArray();
        Hash = hash;
    }

    public bool IsUT
    {
        get
        {
            var n = GL.N;
            var rg = n.Range();
            var t0 = Table;
            return rg.Grid2D(rg).Where(e => e.t1 > e.t2).All(e => t0[e.t1 * n + e.t2] == 0);
        }
    }

    public bool IsLT
    {
        get
        {
            var n = GL.N;
            var rg = n.Range();
            var t0 = Table;
            return rg.Grid2D(rg).Where(e => e.t1 < e.t2).All(e => t0[e.t1 * n + e.t2] == 0);
        }
    }

    public bool IsDiag
    {
        get
        {
            var n = GL.N;
            var rg = n.Range();
            var t0 = Table;
            return rg.Grid2D(rg).Where(e => e.t1 != e.t2).All(e => t0[e.t1 * n + e.t2] == 0);
        }
    }

    public bool Is2ndDiag
    {
        get
        {
            var n = GL.N;
            var rg = n.Range();
            var t0 = Table;
            return rg.Grid2D(rg).Where(e => e.t1 + e.t2 != n - 1).All(e => t0[e.t1 * n + e.t2] == 0);
        }
    }

    public bool IsSym
    {
        get
        {
            var n = GL.N;
            var rg = n.Range();
            var t0 = Table;
            return rg.Grid2D(rg).Where(e => e.t1 < e.t2).All(e => t0[e.t1 * n + e.t2] == t0[e.t2 * n + e.t1]);
        }
    }

    public int Det => GL.Det(this);
    public Mat T 
    {
        get
        {
            var n = GL.N;
            var rg = n.Range();
            var t0 = Table;
            return GL.Create(rg.Grid2D(rg).Select(e => t0[e.t2 * n + e.t1]).ToArray());
        }
    }

    public Mat At(Tuple2Array at, int value) => GL.At(Table, at, value);
    public bool Equals(Mat other) => Hash == other.Hash && Table.SequenceEqual(other.Table);

    public int CompareTo(Mat other) => Table.SequenceCompareTo(other.Table);

    public int Hash { get; }
    public IGroup<Mat> BaseGroup => GL;

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var mod = GL.P;
        return Table.Select(i => new ZnInt(mod, i)).ToKMatrix(GL.N).ToString();
    }
}