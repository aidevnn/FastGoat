using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Matrix;

public struct MatFq : IElt<MatFq>
{
    public EPoly<ZnInt>[] Table { get; }
    public GLnq GLnq { get; }

    public MatFq(GLnq gl, int hash, EPoly<ZnInt>[] table)
    {
        GLnq = gl;
        Table = table.ToArray();
        Hash = hash;
    }

    public bool Equals(MatFq other) => Table.SequenceEqual(other.Table);

    public int CompareTo(MatFq other) => Table.SequenceCompareTo(other.Table);

    public int Hash { get; }
    public override int GetHashCode() => Hash;
    public override string ToString()
    {
        return $"[{Table.Glue("; ", GLnq.Fmt)}]";
    }
}