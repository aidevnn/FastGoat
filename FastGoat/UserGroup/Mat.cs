using System.Text;
namespace FastGoat.UserGroup;

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
        var sb = new StringBuilder();
        for (int i = 0; i < Table.Length; i += GL.N)
            sb.AppendFormat("[{0}]", Table.Skip(i).Take(GL.N).Glue(" ", GL.Fmt));

        return $"[{sb}]";
    }
}
