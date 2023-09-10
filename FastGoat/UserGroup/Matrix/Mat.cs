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

    public Mat At(Tuple2Array at, int value) => GL.At(Table, at, value);
    public bool Equals(Mat other) => Hash == other.Hash;

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