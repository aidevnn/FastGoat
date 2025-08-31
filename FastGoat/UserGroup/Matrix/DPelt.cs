using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Perms;

namespace FastGoat.UserGroup.Matrix;

public readonly struct DPelt : IElt<DPelt>, IEnumerable<int>
{
    public DPGL Dpgl { get; }
    public int[] Diag { get; }
    public int[] Perm { get; }

    public DPelt(DPGL dpgl, int[] diag, int[] perm)
    {
        if (diag.Length != perm.Length || dpgl.N != perm.Length)
            throw new GroupException(GroupExceptionType.GroupDef);

        Diag = diag;
        Perm = perm;
        Dpgl = dpgl;
        Hash = (IntExt.GenHash(dpgl.P, diag), IntExt.GenHash(dpgl.P, perm)).GetHashCode();
    }

    public bool Equals(DPelt other) => Diag.SequenceEqual(other.Diag) && Perm.SequenceEqual(other.Perm);

    public int CompareTo(DPelt other)
    {
        var compDiag = Diag.SequenceCompareTo(other.Diag);
        if (compDiag != 0)
            return compDiag;

        return Perm.SequenceCompareTo(other.Perm);
    }

    public int Hash { get; }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }
    
    public IEnumerator<int> GetEnumerator()
    {
        for (int i = 0; i < Dpgl.N.Pow(2); i++)
            yield return this[i];
    }

    public override int GetHashCode() => Hash;

    public Mat ToGL() => Dpgl.GL.Create(this.ToArray());

    public void DisplayArrays()
    {
        var len = Dpgl.P.ToString().Length;
        Console.WriteLine($"Diag [{Diag.Glue(", ", $"{{0,{len}}}")}]");
        Console.WriteLine($"Perm ({Perm.Select(i => i + 1).Glue(", ", $"{{0,{len}}}")})");
        Console.WriteLine();
    }

    public int Det
    {
        get
        {
            var p = Dpgl.P;
            return Diag.Aggregate((ai, aj) => ai * aj % p);
        }
    }

    public bool IsOrder(int ord)
    {
        if (IntExt.PowMod(Det, ord, Dpgl.P) != 1)
            return false;

        return Group.ElementIsOrder(Dpgl, this, ord);
    }
    public override string ToString()
    {
        var len = Dpgl.P.ToString().Length;
        var fmt = $"({{0, {len}}},{{1, {len}}})";
        var seq = Diag.Zip(Perm).Select(a => string.Format(fmt, a.First, a.Second + 1)).Glue(" ");
        return $"[{seq}]";
    }

    public int this[(int i, int j) index]
    {
        get
        {
            var pj = Perm[index.i];
            return pj == index.j ? Diag[index.i] : 0;
        }
    }

    public int this[int index]
    {
        get
        {
            var n = Dpgl.N;
            var (i, j) = ((index % (n * n)) / n, index % n);
            return this[(i, j)];
        }
    }
}