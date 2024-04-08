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

    public MatFq(GLnq gl, EPoly<ZnInt>[] table)
    {
        GLnq = gl;
        var hash = table.Aggregate(0, (acc, a0) => a0.GetHashCode() + gl.Fq.Q * acc);
        Table = table.ToArray();
        Hash = hash;
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

    public bool Equals(MatFq other) => Table.SequenceEqual(other.Table);

    public int CompareTo(MatFq other) => Table.SequenceCompareTo(other.Table);

    public int Hash { get; }
    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        // return $"[{Table.Glue("; ", GLnq.Fmt)}]";
        return Table.ToKMatrix(GLnq.N).ToString();
    }
}