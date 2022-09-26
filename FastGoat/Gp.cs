namespace FastGoat;

public struct Gp<U1> : IGroup<Ep<U1>> where U1 : struct, IElt<U1>
{
    public Gp(IGroup<U1> ig1)
    {
        g1 = ig1;
        Hash = ig1.Hash;
    }
    public IGroup<U1> g1 { get; }
    public int Hash { get; }
    public bool Equals(IGroup<Ep<U1>>? other) => other?.Hash == Hash;
    public Ep<U1> Invert(Ep<U1> a) => new(this, g1.Invert(a.e1));
    public Ep<U1> Neutral() => new(this, g1.Neutral());
    public Ep<U1> Op(Ep<U1> a, Ep<U1> b) => new(this, g1.Op(a.e1, b.e1));
    public Ep<U1> this[int k] => new(this, g1[k]);
    public override string ToString() => $"{g1}";
}
public struct Gp<U1, U2> : IGroup<Ep<U1, U2>> where U1 : struct, IElt<U1> where U2 : struct, IElt<U2>
{
    public Gp(IGroup<U1> ig1, IGroup<U2> ig2)
    {
        g1 = ig1;
        g2 = ig2;
        Hash = HashCode.Combine(g1.Hash, g2.Hash);
    }
    public IGroup<U1> g1 { get; }
    public IGroup<U2> g2 { get; }
    public int Hash { get; }
    public bool Equals(IGroup<Ep<U1, U2>>? other) => other?.Hash == Hash;
    public Ep<U1, U2> Invert(Ep<U1, U2> a) => new(this, g1.Invert(a.e1), g2.Invert(a.e2));
    public Ep<U1, U2> Neutral() => new(this, g1.Neutral(), g2.Neutral());
    public Ep<U1, U2> Op(Ep<U1, U2> a, Ep<U1, U2> b) => new(this, g1.Op(a.e1, b.e1), g2.Op(a.e2, b.e2));
    public Ep<U1, U2> this[int k] => new(this, g1[k], g2[k]);
    public Ep<U1, U2> this[int k1, int k2] => new(this, g1[k1], g2[k2]);
    public override string ToString() => $"({g1} x {g2})";
}
public struct Gp<U1, U2, U3> : IGroup<Ep<U1, U2, U3>>
    where U1 : struct, IElt<U1>
    where U2 : struct, IElt<U2>
    where U3 : struct, IElt<U3>
{
    public Gp(IGroup<U1> ig1, IGroup<U2> ig2, IGroup<U3> ig3)
    {
        g1 = ig1;
        g2 = ig2;
        g3 = ig3;
        Hash = HashCode.Combine(g1.Hash, g2.Hash, g3.Hash);
    }
    public IGroup<U1> g1 { get; }
    public IGroup<U2> g2 { get; }
    public IGroup<U3> g3 { get; }
    public int Hash { get; }
    public bool Equals(IGroup<Ep<U1, U2, U3>>? other) => other?.Hash == Hash;
    public Ep<U1, U2, U3> Invert(Ep<U1, U2, U3> a) => new(this, g1.Invert(a.e1), g2.Invert(a.e2), g3.Invert(a.e3));
    public Ep<U1, U2, U3> Neutral() => new(this, g1.Neutral(), g2.Neutral(), g3.Neutral());
    public Ep<U1, U2, U3> Op(Ep<U1, U2, U3> a, Ep<U1, U2, U3> b) => new(this, g1.Op(a.e1, b.e1), g2.Op(a.e2, b.e2), g3.Op(a.e3, b.e3));
    public Ep<U1, U2, U3> this[int k] => new(this, g1[k], g2[k], g3[k]);
    public Ep<U1, U2, U3> this[int k1, int k2, int k3] => new(this, g1[k1], g2[k2], g3[k3]);
    public override string ToString() => $"({g1} x {g2} x {g3})";
}
public struct Gp<U1, U2, U3, U4> : IGroup<Ep<U1, U2, U3, U4>>
    where U1 : struct, IElt<U1>
    where U2 : struct, IElt<U2>
    where U3 : struct, IElt<U3>
    where U4 : struct, IElt<U4>
{
    public Gp(IGroup<U1> ig1, IGroup<U2> ig2, IGroup<U3> ig3, IGroup<U4> ig4)
    {
        g1 = ig1;
        g2 = ig2;
        g3 = ig3;
        g4 = ig4;
        Hash = HashCode.Combine(g1.Hash, g2.Hash, g3.Hash, g4.Hash);
    }
    public IGroup<U1> g1 { get; }
    public IGroup<U2> g2 { get; }
    public IGroup<U3> g3 { get; }
    public IGroup<U4> g4 { get; }
    public int Hash { get; }
    public bool Equals(IGroup<Ep<U1, U2, U3, U4>>? other) => other?.Hash == Hash;
    public Ep<U1, U2, U3, U4> Invert(Ep<U1, U2, U3, U4> a) => new(this, g1.Invert(a.e1), g2.Invert(a.e2), g3.Invert(a.e3), g4.Invert(a.e4));
    public Ep<U1, U2, U3, U4> Neutral() => new(this, g1.Neutral(), g2.Neutral(), g3.Neutral(), g4.Neutral());
    public Ep<U1, U2, U3, U4> Op(Ep<U1, U2, U3, U4> a, Ep<U1, U2, U3, U4> b) => new(this, g1.Op(a.e1, b.e1), g2.Op(a.e2, b.e2), g3.Op(a.e3, b.e3), g4.Op(a.e4, b.e4));
    public Ep<U1, U2, U3, U4> this[int k] => new(this, g1[k], g2[k], g3[k], g4[k]);
    public Ep<U1, U2, U3, U4> this[int k1, int k2, int k3, int k4] => new(this, g1[k1], g2[k2], g3[k3], g4[k4]);
    public override string ToString() => $"({g1} x {g2} x {g3} x {g4})";
}
public struct Gp<U1, U2, U3, U4, U5> : IGroup<Ep<U1, U2, U3, U4, U5>>
    where U1 : struct, IElt<U1>
    where U2 : struct, IElt<U2>
    where U3 : struct, IElt<U3>
    where U4 : struct, IElt<U4>
    where U5 : struct, IElt<U5>
{
    public Gp(IGroup<U1> ig1, IGroup<U2> ig2, IGroup<U3> ig3, IGroup<U4> ig4, IGroup<U5> ig5)
    {
        g1 = ig1;
        g2 = ig2;
        g3 = ig3;
        g4 = ig4;
        g5 = ig5;
        Hash = HashCode.Combine(g1.Hash, g2.Hash, g3.Hash, g4.Hash, g5.Hash);
    }
    public IGroup<U1> g1 { get; }
    public IGroup<U2> g2 { get; }
    public IGroup<U3> g3 { get; }
    public IGroup<U4> g4 { get; }
    public IGroup<U5> g5 { get; }
    public int Hash { get; }
    public bool Equals(IGroup<Ep<U1, U2, U3, U4, U5>>? other) => other?.Hash == Hash;
    public Ep<U1, U2, U3, U4, U5> Invert(Ep<U1, U2, U3, U4, U5> a) => new(this, g1.Invert(a.e1), g2.Invert(a.e2), g3.Invert(a.e3), g4.Invert(a.e4), g5.Invert(a.e5));
    public Ep<U1, U2, U3, U4, U5> Neutral() => new(this, g1.Neutral(), g2.Neutral(), g3.Neutral(), g4.Neutral(), g5.Neutral());
    public Ep<U1, U2, U3, U4, U5> Op(Ep<U1, U2, U3, U4, U5> a, Ep<U1, U2, U3, U4, U5> b) => new(this, g1.Op(a.e1, b.e1), g2.Op(a.e2, b.e2), g3.Op(a.e3, b.e3), g4.Op(a.e4, b.e4), g5.Op(a.e5, b.e5));
    public Ep<U1, U2, U3, U4, U5> this[int k] => new(this, g1[k], g2[k], g3[k], g4[k], g5[k]);
    public Ep<U1, U2, U3, U4, U5> this[int k1, int k2, int k3, int k4, int k5] => new(this, g1[k1], g2[k2], g3[k3], g4[k4], g5[k5]);
    public override string ToString() => $"({g1} x {g2} x {g3} x {g4} x {g5})";
}
