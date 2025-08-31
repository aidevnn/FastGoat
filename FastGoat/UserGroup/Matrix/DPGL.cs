using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Perms;

namespace FastGoat.UserGroup.Matrix;

public readonly struct DPGL : IGroup<DPelt>
{
    public GL GL { get; }
    public Sn Sn { get; }
    public int N => GL.N;
    public int P => GL.P;

    public DPGL(int n, int p)
    {
        GL = new GL(n, p);
        Sn = new Sn(n);
        Name = $"DPGL({n},{p})";
        Hash = (n, p, "DPGL").GetHashCode();
    }

    public IEnumerator<DPelt> GetEnumerator()=> GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<DPelt>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public DPelt this[params ValueType[] us]
    {
        get
        {
            var p = P;
            var n = N;
            var diag = us.Select(i => IntExt.AmodP(Convert.ToInt32(i), p))
                .Concat(Enumerable.Repeat(1, n))
                .Take(N)
                .ToArray();
            return new(this, diag, N.Range());
        }
    }

    public DPelt this[int[] diag, params ValueType[] perm]
    {
        get
        {
            var p = P;
            var n = N;
            var diag0 = diag.Select(i => IntExt.AmodP(i, p))
                .Concat(Enumerable.Repeat(1, n))
                .Take(N)
                .ToArray();
            return new(this, diag0, Sn[perm].Table);
        }
    }

    public DPelt this[int[] diag, Perm perm]
    {
        get
        {
            var p = P;
            var n = N;
            var diag0 = diag.Select(i => IntExt.AmodP(i, p))
                .Concat(Enumerable.Repeat(1, n))
                .Take(N)
                .ToArray();
            return new(this, diag0, perm.Table);
        }
    }

    public DPelt this[Perm perm]
    {
        get
        {
            var diag = Enumerable.Repeat(1, N).ToArray();
            return new(this, diag, perm.Table);
        }
    }

    public IEnumerable<DPelt> GetElements()
    {
        yield return Neutral();
    }

    public IEnumerable<DPelt> GetGenerators()
    {
        yield return Neutral();
    }

    public DPelt Neutral() => new(this, Enumerable.Repeat(1, N).ToArray(), N.Range());

    public DPelt Invert(DPelt e)
    {
        var p = P;
        var perm = new int[N];
        var diag = new int[N];
        for (int i = 0; i < N; i++)
        {
            var j = e.Perm[i];
            (perm[j], diag[j]) = (i, IntExt.InvModPbez(e.Diag[i], p));
        }
        
        return new(this, diag, perm);
    }

    public DPelt Op(DPelt e1, DPelt e2)
    {
        var diag = new int[N];
        var perm = new int[N];
        for (int i = 0; i < N; i++)
        {
            var e1pi = e1.Perm[i];
            (diag[i], perm[i]) = (e1.Diag[i] * e2.Diag[e1pi] % P, e2.Perm[e1pi]);
        }
        
        return new(this, diag, perm);
    }
    
    public override string ToString() => Name;
    public override int GetHashCode() => Hash;

}