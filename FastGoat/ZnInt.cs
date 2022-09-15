// using FiniteGroupsOld;

// var z72 = new Zn(72);
// var g0 = new SubGroup(z72, new ZnInt(z72, 8));
// Console.WriteLine(g0.Glue("; "));

// var g1 = new SubGroup(z72, new ZnInt(z72, 8), new ZnInt(z72, 6));
// Console.WriteLine(g1.Glue("; "));

// var z20 = new Zn(20);
// var z30 = new Zn(30);
// var z20x30 = new Group2(z20, z30);
// var g2 = new SubGroup(z20x30, new Elt2(z20x30, new ZnInt(z20, 5), new ZnInt(z30, 6)));
// Console.WriteLine(g2.Glue("; "));

// var g3 = new SubGroup(g2, new Elt2(z20x30, new ZnInt(z20, 5), new ZnInt(z30, 0)));
// Console.WriteLine(g3.Glue("; "));

// // var z8x18x30 = new ZnTuple(8, 18, 30);
// // var g4 = new SubGroup(z8x18x30, new ZnIntTuple(z8x18x30, 1, 0, 0), new ZnIntTuple(z8x18x30, 0, 1, 0), new ZnIntTuple(z8x18x30, 0, 0, 1));
// // Console.WriteLine(g4.Count());

// var z8 = new Zn(8);
// var z18 = new Zn(18);
// var gz8x18x30 = new Group3(z8, z18, z30);
// var e1 = new Elt3(gz8x18x30, new ZnInt(z8, 1), new ZnInt(z18, 0), new ZnInt(z30, 0));
// var e2 = new Elt3(gz8x18x30, new ZnInt(z8, 0), new ZnInt(z18, 1), new ZnInt(z30, 0));
// var e3 = new Elt3(gz8x18x30, new ZnInt(z8, 0), new ZnInt(z18, 0), new ZnInt(z30, 1));
// var g5 = new SubGroup(gz8x18x30, e1, e2, e3);
// Console.WriteLine(g5.Count());

using FastGoat;

public struct ZnInt : IElt<ZnInt>
{
    public int k { get; }
    public ZnInt()
    {
        throw new Exception();
    }
    public ZnInt(Zn zn)
    {
        k = 0;
        Group = Zn = zn;
        Hash = HashCode.Combine(zn.Hash, k);
    }
    public ZnInt(Zn zn, int k0)
    {
        k = k0 % zn.mod;
        if (k < 0)
            k += zn.mod;

        Group = Zn = zn;
        Hash = HashCode.Combine(zn, k);
    }
    public IGroup<ZnInt> Group { get; }
    public Zn Zn { get; }

    public int Hash { get; }

    public int CompareTo(ZnInt other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return k.CompareTo(other.k);
    }

    public bool Equals(ZnInt other) => Group.Equals(other.Group) && k == other.k;
    public override int GetHashCode() => Hash;
    public override string ToString() => string.Format(Zn.fmt, k);

    public static implicit operator ZnInt((Zn zn, int k) p) => new(p.zn, p.k);
    public static ZnInt operator *(ZnInt a, ZnInt b) => a.Group.Op(a, b);
    public static ZnInt operator ^(ZnInt a, int p) => (a.Zn, a.k * p);
}
