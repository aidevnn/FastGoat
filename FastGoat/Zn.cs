using FastGoat.Structures.GroupTheory;
using FastGoat.Structures.SetTheory;

namespace FastGoat;

public struct ZnElt : IElt<ZnElt>
{
    public int HashCode { get; }
    public int[] Table { get; }
    public IFSet<ZnElt> FSet { get; }

    public ZnElt(Zn zn)
    {
        FSet = zn;
        Table = new int[zn.Dims.Length];
        HashCode = 0;
    }

    public ZnElt(Zn zn, int[] arr, int hash)
    {
        FSet = zn;
        Table = arr.ToArray();
        HashCode = hash;
    }

    public int CompareTo(ZnElt other) => Helpers.ArrayCompare(Table, other.Table);
    public bool Equals(ZnElt other) => HashCode == other.HashCode;
    public override int GetHashCode() => HashCode;
    public override string ToString() => $"({Table.Glue(sep: ",")})";

}

public class Zn : Group<ZnElt>
{
    public int[] Dims { get; }
    readonly int[] cache;

    public Zn(params int[] dims)
    {
        Dims = dims;
        cache = new int[dims.Length];
    }

    public override ZnElt Neutral => new ZnElt(this);
    public override ZnElt Invert(ZnElt a)
    {
        var hash = Helpers.InvertModulo(Dims, a.Table, cache);
        return new ZnElt(this, cache, hash);
    }

    public override ZnElt Op(ZnElt a, ZnElt b)
    {
        var hash = Helpers.AddModulo(Dims, a.Table, b.Table, cache);
        return new ZnElt(this, cache, hash);
    }

    public ZnElt CreateElement(params int[] vs)
    {
        Helpers.ClearArray(cache);
        int hash = Helpers.AddModulo(Dims, vs, Neutral.Table, cache);
        return new ZnElt(this, cache, hash);
    }

    public ZnElt CE(params int[] vs) => CreateElement(vs);
    public ZnElt[] BaseCanonic => Helpers.BaseCanonic(Dims.Length).Select(CE).ToArray();
    public SubGroup<ZnElt> GenerateAll()
    {
        var g = this.GroupElement(BaseCanonic).Generate();
        g.SetName(Dims.Glue(fmt: "C{0}", sep: " x "));
        return g;
    }
}