namespace FastGoat;

public struct Ep<U1> : IElt<Ep<U1>> where U1 : struct, IElt<U1>
{
    public Ep(U1 ie1)
    {
        e1 = ie1;
        Group = new Gp<U1>(e1.Group);
        Hash = e1.Hash;
    }
    public Ep(Gp<U1> gr, U1 ie1)
    {
        e1 = ie1;
        Group = gr;
        Hash = e1.Hash;
    }
    public U1 e1 { get; }
    public IGroup<Ep<U1>> Group { get; }

    public int Hash { get; }

    public int CompareTo(Ep<U1> other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return e1.CompareTo(other.e1);
    }

    public bool Equals(Ep<U1> other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return e1.Equals(other.e1);
    }
    public override int GetHashCode() => Hash;
    public override string ToString() => e1.ToString();

    public static implicit operator Ep<U1>(U1 e) => new(e);
    public static Ep<U1> operator *(Ep<U1> a, Ep<U1> b) => a.Group.Op(a, b);
    public static Ep<U1> operator *(Ep<U1> a, int k)
    {
        var g = a.Group;
        if (k == 0)
            return g.Neutral();

        if (k < 0)
        {
            var ai = g.Invert(a);
            return Enumerable.Repeat(ai, -k).Aggregate((e0, e1) => g.Op(e0, e1));
        }

        return Enumerable.Repeat(a, k).Aggregate((e0, e1) => g.Op(e0, e1));
    }
}
public struct Ep<U1, U2> : IElt<Ep<U1, U2>> where U1 : struct, IElt<U1> where U2 : struct, IElt<U2>
{
    public Ep(U1 ie1, U2 ie2)
    {
        e1 = ie1;
        e2 = ie2;
        Group = new Gp<U1, U2>(e1.Group, e2.Group);
        Hash = HashCode.Combine(e1.Hash, e2.Hash);
    }
    public Ep(Gp<U1, U2> gr, U1 ie1, U2 ie2)
    {
        e1 = ie1;
        e2 = ie2;
        Group = gr;
        Hash = HashCode.Combine(e1.Hash, e2.Hash);
    }
    public U1 e1 { get; }
    public U2 e2 { get; }
    public IGroup<Ep<U1, U2>> Group { get; }

    public int Hash { get; }

    public int CompareTo(Ep<U1, U2> other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return (e1, e2).CompareTo((other.e1, other.e2));
    }

    public bool Equals(Ep<U1, U2> other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return (e1, e2).Equals((other.e1, other.e2));
    }
    public override int GetHashCode() => Hash;
    public override string ToString() => (e1, e2).ToString();

    public static implicit operator Ep<U1, U2>((U1, U2) e) => new(e.Item1, e.Item2);
}
public struct Ep<U1, U2, U3> : IElt<Ep<U1, U2, U3>>
    where U1 : struct, IElt<U1>
    where U2 : struct, IElt<U2>
    where U3 : struct, IElt<U3>
{
    public Ep(U1 ie1, U2 ie2, U3 ie3)
    {
        e1 = ie1;
        e2 = ie2;
        e3 = ie3;
        Group = new Gp<U1, U2, U3>(e1.Group, e2.Group, e3.Group);
        Hash = HashCode.Combine(e1.Hash, e2.Hash, e3.Hash);
    }
    public Ep(Gp<U1, U2, U3> gr, U1 ie1, U2 ie2, U3 ie3)
    {
        e1 = ie1;
        e2 = ie2;
        e3 = ie3;
        Group = gr;
        Hash = HashCode.Combine(e1.Hash, e2.Hash, e3.Hash);
    }
    public U1 e1 { get; }
    public U2 e2 { get; }
    public U3 e3 { get; }
    public IGroup<Ep<U1, U2, U3>> Group { get; }

    public int Hash { get; }

    public int CompareTo(Ep<U1, U2, U3> other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return (e1, e2, e3).CompareTo((other.e1, other.e2, other.e3));
    }

    public bool Equals(Ep<U1, U2, U3> other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return (e1, e2, e3).Equals((other.e1, other.e2, other.e3));
    }
    public override int GetHashCode() => Hash;
    public override string ToString() => (e1, e2, e3).ToString();
    public static implicit operator Ep<U1, U2, U3>((U1, U2, U3) e) => new(e.Item1, e.Item2, e.Item3);
}
public struct Ep<U1, U2, U3, U4> : IElt<Ep<U1, U2, U3, U4>>
    where U1 : struct, IElt<U1>
    where U2 : struct, IElt<U2>
    where U3 : struct, IElt<U3>
    where U4 : struct, IElt<U4>
{
    public Ep(U1 ie1, U2 ie2, U3 ie3, U4 ie4)
    {
        e1 = ie1;
        e2 = ie2;
        e3 = ie3;
        e4 = ie4;
        Group = new Gp<U1, U2, U3, U4>(e1.Group, e2.Group, e3.Group, e4.Group);
        Hash = HashCode.Combine(e1.Hash, e2.Hash, e3.Hash, e4.Hash);
    }
    public Ep(Gp<U1, U2, U3, U4> gr, U1 ie1, U2 ie2, U3 ie3, U4 ie4)
    {
        e1 = ie1;
        e2 = ie2;
        e3 = ie3;
        e4 = ie4;
        Group = gr;
        Hash = HashCode.Combine(e1.Hash, e2.Hash, e3.Hash, e4.Hash);
    }
    public U1 e1 { get; }
    public U2 e2 { get; }
    public U3 e3 { get; }
    public U4 e4 { get; }
    public IGroup<Ep<U1, U2, U3, U4>> Group { get; }
    public int Hash { get; }
    public int CompareTo(Ep<U1, U2, U3, U4> other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return (e1, e2, e3, e4).CompareTo((other.e1, other.e2, other.e3, other.e4));
    }
    public bool Equals(Ep<U1, U2, U3, U4> other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return (e1, e2, e3, e4).Equals((other.e1, other.e2, other.e3, other.e4));
    }
    public override int GetHashCode() => Hash;
    public override string ToString() => (e1, e2, e3, e4).ToString();
    public static implicit operator Ep<U1, U2, U3, U4>((U1, U2, U3, U4) e) => new(e.Item1, e.Item2, e.Item3, e.Item4);
}
public struct Ep<U1, U2, U3, U4, U5> : IElt<Ep<U1, U2, U3, U4, U5>>
    where U1 : struct, IElt<U1>
    where U2 : struct, IElt<U2>
    where U3 : struct, IElt<U3>
    where U4 : struct, IElt<U4>
    where U5 : struct, IElt<U5>
{
    public Ep(U1 ie1, U2 ie2, U3 ie3, U4 ie4, U5 ie5)
    {
        e1 = ie1;
        e2 = ie2;
        e3 = ie3;
        e4 = ie4;
        e5 = ie5;
        Group = new Gp<U1, U2, U3, U4, U5>(e1.Group, e2.Group, e3.Group, e4.Group, e5.Group);
        Hash = HashCode.Combine(e1.Hash, e2.Hash, e3.Hash, e4.Hash, e5.Hash);
    }
    public Ep(Gp<U1, U2, U3, U4, U5> gr, U1 ie1, U2 ie2, U3 ie3, U4 ie4, U5 ie5)
    {
        e1 = ie1;
        e2 = ie2;
        e3 = ie3;
        e4 = ie4;
        e5 = ie5;
        Group = gr;
        Hash = HashCode.Combine(e1.Hash, e2.Hash, e3.Hash, e4.Hash, e5.Hash);
    }
    public U1 e1 { get; }
    public U2 e2 { get; }
    public U3 e3 { get; }
    public U4 e4 { get; }
    public U5 e5 { get; }
    public IGroup<Ep<U1, U2, U3, U4, U5>> Group { get; }
    public int Hash { get; }
    public int CompareTo(Ep<U1, U2, U3, U4, U5> other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return (e1, e2, e3, e4, e5).CompareTo((other.e1, other.e2, other.e3, other.e4, other.e5));
    }
    public bool Equals(Ep<U1, U2, U3, U4, U5> other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return (e1, e2, e3, e4, e5).Equals((other.e1, other.e2, other.e3, other.e4, other.e5));
    }
    public override int GetHashCode() => Hash;
    public override string ToString() => (e1, e2, e3, e4, e5).ToString();
    public static implicit operator Ep<U1, U2, U3, U4, U5>((U1, U2, U3, U4, U5) e) => new(e.Item1, e.Item2, e.Item3, e.Item4, e.Item5);
}
