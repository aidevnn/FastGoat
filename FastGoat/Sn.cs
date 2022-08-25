using FastGoat.Structures.GroupTheory;
using FastGoat.Structures.SetTheory;

namespace FastGoat;

public struct Perm : IElt<Perm>
{
    public int HashCode { get; }
    public int[] Table { get; }
    public IFSet<Perm> FSet { get; }

    public int Sgn { get; }
    public int[] Invariants { get; }
    public int[][] Cycles { get; }

    public Perm(Sn sn)
    {
        FSet = sn;
        Table = sn.N.Range();
        HashCode = Helpers.GenHash(sn.N, Table);
        Invariants = Table.ToArray();
        Sgn = 1;
        Cycles = new int[0][];
    }

    public Perm(Sn sn, int[] arr, int hash)
    {
        FSet = sn;
        Table = arr.ToArray();
        HashCode = hash;
        var orbits = Table.Orbits();
        Invariants = orbits.Where(a => a.Count == 1).SelectMany(a => a).ToArray();
        Cycles = orbits.Where(a => a.Count > 1).Select(a => a.ToArray()).ToArray();
        Sgn = (int)Math.Pow(-1, sn.N - orbits.Count);
    }

    public int CompareTo(Perm other) => Helpers.ArrayCompare(Table, other.Table);
    public bool Equals(Perm other) => HashCode == other.HashCode;
    public override int GetHashCode() => HashCode;
    public override string ToString()
    {
        var inv = Invariants.Length == 0 ? "" : string.Format(" Invariants : [{0}]", Invariants.Select(b => b + 1).Glue(fmt: "{0}"));
        var strCycles = Cycles.Select(a => string.Format("({0})", a.Select(b => b + 1).Glue(fmt: "{0}"))).Glue();
        var cycles = Cycles.Length == 0 ? "" : string.Format(" Cycles : {0}", strCycles);
        var sgn = Sgn == 1 ? "(+)" : "(-)";
        return string.Format("[{0}]{1}{2}{3}", Table.Select(i => i + 1).Glue(sep: ""), sgn, inv, cycles);
    }
}

public class Sn : Group<Perm>
{
    public int N { get; }
    readonly int[] cache;
    public Sn(int n)
    {
        N = n;
        cache = new int[n];
    }

    public override Perm Neutral => new Perm(this);
    public override Perm Invert(Perm a)
    {
        var hash = Helpers.InvertPermutation(a.Table, cache);
        return new Perm(this, cache, hash);
    }

    public override Perm Op(Perm a, Perm b)
    {
        var hash = Helpers.ComposePermutation(a.Table, b.Table, cache);
        return new Perm(this, cache, hash);
    }

    public Perm CreateElement(params int[] vs)
    {
        vs.Add(-1);
        if (!Helpers.CheckTable(N, vs))
            return Neutral;

        var hash = Helpers.GenHash(N, vs);
        var p = new Perm(this, vs, hash);
        AddElement(p);
        return p;
    }

    public Perm Cycle(params int[] vs)
    {
        vs.Add(-1);
        if (!Helpers.CheckCycle(N, vs))
            return Neutral;

        Neutral.Table.CopyTo(cache, 0);
        Helpers.ApplyCycle(cache, vs);
        var hash = Helpers.GenHash(N, cache);
        return new Perm(this, cache, hash);
    }

    public Perm KCycle(int count) => Cycle(count.Range(1));
    public Perm KCycle(int start, int count) => Cycle(count.Range(start));

    public Perm Cycle(Tuple2Array cycle) => Cycle(cycle.Table);
    public Perm ComposesCycles(params Tuple2Array[] cycles)
    {
        Neutral.Table.CopyTo(cache, 0);
        foreach (var e in cycles)
        {
            e.Table.Add(-1);
            if (!Helpers.CheckCycle(N, e.Table)) return Neutral;

            Helpers.ApplyCycle(cache, e.Table);
        }

        var hash = Helpers.GenHash(N, cache);
        return new Perm(this, cache, hash);
    }

    public Perm C(params int[] vs) => Cycle(vs);
    public Perm C(params Tuple2Array[] cycles) => ComposesCycles(cycles);
    public Perm[] AllPerm => Helpers.AllPerms(N, true).Select(CreateElement).ToArray();

    public SubGroup<Perm> GenerateAll()
    {
        var h = AllPerm;
        var g = this.GroupElement();
        g.SetName($"S{N}");
        return g;
    }
}