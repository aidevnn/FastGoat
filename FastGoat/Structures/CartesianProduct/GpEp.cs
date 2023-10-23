using System.Collections;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.CartesianProduct;

public readonly struct Gp2<T1, T2> : IGroup<Ep2<T1, T2>> where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    public string Name { get; }
    public IGroup<T1> G1 { get; }
    public IGroup<T2> G2 { get; }

    public Gp2(IGroup<T1> g1, IGroup<T2> g2)
    {
        G1 = g1;
        G2 = g2;
        Hash = (g1.Hash, g2.Hash).GetHashCode();
        Name = $"{G1} x {G2}";
    }

    public bool Equals(IGroup<Ep2<T1, T2>>? other) => other?.Hash == Hash;
    public int Hash { get; }
    public Ep2<T1, T2> Neutral() => new(G1.Neutral(), G2.Neutral());
    public Ep2<T1, T2> Invert(Ep2<T1, T2> e) => new(G1.Invert(e.E1), G2.Invert(e.E2));
    public Ep2<T1, T2> Op(Ep2<T1, T2> e1, Ep2<T1, T2> e2) => new(G1.Op(e1.E1, e2.E1), G2.Op(e1.E2, e2.E2));

    public Ep2<T1, T2> this[params ValueType[] us]
    {
        get
        {
            dynamic us0 = new ValueType[2];
            if (us.Length == 1)
            {
                dynamic us1 = us[0];
                us0[0] = us1.Item1;
                us0[1] = us1.Item2;
            }
            else if (us.Length == 2)
                us0 = us;
            else
                throw new GroupException(GroupExceptionType.GroupDef);

            var e1 = G1[us0[0]];
            var e2 = G2[us0[1]];
            return new(e1, e2);
        }
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;

    public IEnumerable<Ep2<T1, T2>> GetGenerators()
    {
        foreach (var e in G1.GetGenerators())
            yield return new Ep2<T1, T2>(e, G2.Neutral());

        foreach (var e in G2.GetGenerators())
            yield return new Ep2<T1, T2>(G1.Neutral(), e);
    }

    public IEnumerable<Ep2<T1, T2>> GetElements()
    {
        foreach (var e1 in G1)
        foreach (var e2 in G2)
            yield return new(e1, e2);
    }

    public IEnumerator<Ep2<T1, T2>> GetEnumerator() => GetElements().GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();
}

public readonly struct Ep2<T1, T2> : IElt<Ep2<T1, T2>> where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    public T1 E1 { get; }
    public T2 E2 { get; }

    public Ep2(T1 e1, T2 e2)
    {
        E1 = e1;
        E2 = e2;
        Hash = (e1.Hash, e2.Hash).GetHashCode();
    }

    public bool Equals(Ep2<T1, T2> other) => other.Hash == Hash && (E1, E2).Equals((other.E1, other.E2));

    public int CompareTo(Ep2<T1, T2> other)
    {
        return (E1, E2).CompareTo((other.E1, other.E2));
    }

    public int Hash { get; }
    public override int GetHashCode() => Hash;
    public override string ToString() => $"({E1}, {E2})";

    public void Deconstruct(out T1 e1, out T2 e2)
    {
        (e1, e2) = (E1, E2);
    }
}

public static partial class Product
{
    public static Gp2<T1, T2> Group<T1, T2>(IGroup<T1> g1, IGroup<T2> g2)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return new(g1, g2);
    }

    public static Ep2<T1, T2> Elt<T1, T2>(T1 e1, T2 e2)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return new(e1, e2);
    }

    public static ConcreteGroup<Ep2<T1, T2>> Generate<T1, T2>(IGroup<T1> g1, IGroup<T2> g2)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return new(Group(g1, g2));
    }

    public static ConcreteGroup<Ep2<T1, T2>> Generate<T1, T2>(string name, IGroup<T1> g1, IGroup<T2> g2)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return new(name, Group(g1, g2));
    }
}

public readonly struct Gp3<T1, T2, T3> : IGroup<Ep3<T1, T2, T3>> where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
    where T3 : struct, IElt<T3>
{
    public string Name { get; }
    public IGroup<T1> G1 { get; }
    public IGroup<T2> G2 { get; }
    public IGroup<T3> G3 { get; }

    public Gp3(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3)
    {
        G1 = g1;
        G2 = g2;
        G3 = g3;
        Hash = (g1.Hash, g2.Hash, g3.Hash).GetHashCode();
        Name = $"{G1} x {G2} x {G3}";
    }

    public bool Equals(IGroup<Ep3<T1, T2, T3>>? other) => other?.Hash == Hash;
    public int Hash { get; }
    public Ep3<T1, T2, T3> Neutral() => new(G1.Neutral(), G2.Neutral(), G3.Neutral());
    public Ep3<T1, T2, T3> Invert(Ep3<T1, T2, T3> e) => new(G1.Invert(e.E1), G2.Invert(e.E2), G3.Invert(e.E3));

    public Ep3<T1, T2, T3> Op(Ep3<T1, T2, T3> e1, Ep3<T1, T2, T3> e2) =>
        new(G1.Op(e1.E1, e2.E1), G2.Op(e1.E2, e2.E2), G3.Op(e1.E3, e2.E3));

    public Ep3<T1, T2, T3> this[params ValueType[] us]
    {
        get
        {
            dynamic us0 = new ValueType[3];
            if (us.Length == 1)
            {
                dynamic us1 = us[0];
                us0[0] = us1.Item1;
                us0[1] = us1.Item2;
                us0[2] = us1.Item3;
            }
            else if (us.Length == 3)
                us0 = us;
            else
                throw new GroupException(GroupExceptionType.GroupDef);

            var e1 = G1[us0[0]];
            var e2 = G2[us0[1]];
            var e3 = G3[us0[2]];
            return new(e1, e2, e3);
        }
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;

    public IEnumerable<Ep3<T1, T2, T3>> GetGenerators()
    {
        foreach (var e in G1.GetGenerators())
            yield return new Ep3<T1, T2, T3>(e, G2.Neutral(), G3.Neutral());

        foreach (var e in G2.GetGenerators())
            yield return new Ep3<T1, T2, T3>(G1.Neutral(), e, G3.Neutral());

        foreach (var e in G3.GetGenerators())
            yield return new Ep3<T1, T2, T3>(G1.Neutral(), G2.Neutral(), e);
    }

    public IEnumerable<Ep3<T1, T2, T3>> GetElements()
    {
        foreach (var e1 in G1)
        foreach (var e2 in G2)
        foreach (var e3 in G3)
            yield return new(e1, e2, e3);
    }

    public IEnumerator<Ep3<T1, T2, T3>> GetEnumerator() => GetElements().GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();
}

public readonly struct Ep3<T1, T2, T3> : IElt<Ep3<T1, T2, T3>> where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
    where T3 : struct, IElt<T3>
{
    public T1 E1 { get; }
    public T2 E2 { get; }
    public T3 E3 { get; }

    public Ep3(T1 e1, T2 e2, T3 e3)
    {
        E1 = e1;
        E2 = e2;
        E3 = e3;
        Hash = (e1.Hash, e2.Hash, e3.Hash).GetHashCode();
    }

    public bool Equals(Ep3<T1, T2, T3> other) => other.Hash == Hash && (E1, E2, E3).Equals((other.E1, other.E2, other.E3));

    public int CompareTo(Ep3<T1, T2, T3> other)
    {
        return (E1, E2, E3).CompareTo((other.E1, other.E2, other.E3));
    }

    public int Hash { get; }
    public override int GetHashCode() => Hash;
    public override string ToString() => $"({E1}, {E2}, {E3})";

    public void Deconstruct(out T1 e1, out T2 e2, out T3 e3)
    {
        (e1, e2, e3) = (E1, E2, E3);
    }
}

public static partial class Product
{
    public static Gp3<T1, T2, T3> Group<T1, T2, T3>(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2> where T3 : struct, IElt<T3>
    {
        return new(g1, g2, g3);
    }

    public static Ep3<T1, T2, T3> Elt<T1, T2, T3>(T1 e1, T2 e2, T3 e3)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2> where T3 : struct, IElt<T3>
    {
        return new(e1, e2, e3);
    }

    public static ConcreteGroup<Ep3<T1, T2, T3>> Generate<T1, T2, T3>(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2> where T3 : struct, IElt<T3>
    {
        return new(Group(g1, g2, g3));
    }

    public static ConcreteGroup<Ep3<T1, T2, T3>> Generate<T1, T2, T3>(string name, IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2> where T3 : struct, IElt<T3>
    {
        return new(name, Group(g1, g2, g3));
    }
}

public readonly struct Gp4<T1, T2, T3, T4> : IGroup<Ep4<T1, T2, T3, T4>> where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
    where T3 : struct, IElt<T3>
    where T4 : struct, IElt<T4>
{
    public string Name { get; }
    public IGroup<T1> G1 { get; }
    public IGroup<T2> G2 { get; }
    public IGroup<T3> G3 { get; }
    public IGroup<T4> G4 { get; }

    public Gp4(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3, IGroup<T4> g4)
    {
        G1 = g1;
        G2 = g2;
        G3 = g3;
        G4 = g4;
        Hash = (g1.Hash, g2.Hash, g3.Hash, g4.Hash).GetHashCode();
        Name = $"{G1} x {G2} x {G3} x {G4}";
    }

    public bool Equals(IGroup<Ep4<T1, T2, T3, T4>>? other) => other?.Hash == Hash;
    public int Hash { get; }
    public Ep4<T1, T2, T3, T4> Neutral() => new(G1.Neutral(), G2.Neutral(), G3.Neutral(), G4.Neutral());

    public Ep4<T1, T2, T3, T4> Invert(Ep4<T1, T2, T3, T4> e) =>
        new(G1.Invert(e.E1), G2.Invert(e.E2), G3.Invert(e.E3), G4.Invert(e.E4));

    public Ep4<T1, T2, T3, T4> Op(Ep4<T1, T2, T3, T4> e1, Ep4<T1, T2, T3, T4> e2) =>
        new(G1.Op(e1.E1, e2.E1), G2.Op(e1.E2, e2.E2), G3.Op(e1.E3, e2.E3), G4.Op(e1.E4, e2.E4));

    public Ep4<T1, T2, T3, T4> this[params ValueType[] us]
    {
        get
        {
            dynamic us0 = new ValueType[4];
            if (us.Length == 1)
            {
                dynamic us1 = us[0];
                us0[0] = us1.Item1;
                us0[1] = us1.Item2;
                us0[2] = us1.Item3;
                us0[3] = us1.Item4;
            }
            else if (us.Length == 4)
                us0 = us;
            else
                throw new GroupException(GroupExceptionType.GroupDef);

            var e1 = G1[us0[0]];
            var e2 = G2[us0[1]];
            var e3 = G3[us0[2]];
            var e4 = G4[us0[3]];
            return new(e1, e2, e3, e4);
        }
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;

    public IEnumerable<Ep4<T1, T2, T3, T4>> GetGenerators()
    {
        foreach (var e in G1.GetGenerators())
            yield return new Ep4<T1, T2, T3, T4>(e, G2.Neutral(), G3.Neutral(), G4.Neutral());

        foreach (var e in G2.GetGenerators())
            yield return new Ep4<T1, T2, T3, T4>(G1.Neutral(), e, G3.Neutral(), G4.Neutral());

        foreach (var e in G3.GetGenerators())
            yield return new Ep4<T1, T2, T3, T4>(G1.Neutral(), G2.Neutral(), e, G4.Neutral());

        foreach (var e in G4.GetGenerators())
            yield return new Ep4<T1, T2, T3, T4>(G1.Neutral(), G2.Neutral(), G3.Neutral(), e);
    }

    public IEnumerable<Ep4<T1, T2, T3, T4>> GetElements()
    {
        foreach (var e1 in G1)
        foreach (var e2 in G2)
        foreach (var e3 in G3)
        foreach (var e4 in G4)
            yield return new(e1, e2, e3, e4);
    }

    public IEnumerator<Ep4<T1, T2, T3, T4>> GetEnumerator() => GetElements().GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();
}

public readonly struct Ep4<T1, T2, T3, T4> : IElt<Ep4<T1, T2, T3, T4>> where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
    where T3 : struct, IElt<T3>
    where T4 : struct, IElt<T4>
{
    public T1 E1 { get; }
    public T2 E2 { get; }
    public T3 E3 { get; }
    public T4 E4 { get; }

    public Ep4(T1 e1, T2 e2, T3 e3, T4 e4)
    {
        E1 = e1;
        E2 = e2;
        E3 = e3;
        E4 = e4;
        Hash = (e1.Hash, e2.Hash, e3.Hash, e4.Hash).GetHashCode();
    }

    public bool Equals(Ep4<T1, T2, T3, T4> other) =>
        other.Hash == Hash && (E1, E2, E3, E4).Equals((other.E1, other.E2, other.E3, other.E4));

    public int CompareTo(Ep4<T1, T2, T3, T4> other)
    {
        return (E1, E2, E3, E4).CompareTo((other.E1, other.E2, other.E3, other.E4));
    }

    public int Hash { get; }
    public override int GetHashCode() => Hash;
    public override string ToString() => $"({E1}, {E2}, {E3}, {E4})";

    public void Deconstruct(out T1 e1, out T2 e2, out T3 e3, out T4 e4)
    {
        (e1, e2, e3, e4) = (E1, E2, E3, E4);
    }
}

public static partial class Product
{
    public static Gp4<T1, T2, T3, T4> Group<T1, T2, T3, T4>(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3, IGroup<T4> g4)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2> where T3 : struct, IElt<T3> where T4 : struct, IElt<T4>
    {
        return new(g1, g2, g3, g4);
    }

    public static Ep4<T1, T2, T3, T4> Elt<T1, T2, T3, T4>(T1 e1, T2 e2, T3 e3, T4 e4)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2> where T3 : struct, IElt<T3> where T4 : struct, IElt<T4>
    {
        return new(e1, e2, e3, e4);
    }

    public static ConcreteGroup<Ep4<T1, T2, T3, T4>> Generate<T1, T2, T3, T4>(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3,
        IGroup<T4> g4)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2> where T3 : struct, IElt<T3> where T4 : struct, IElt<T4>
    {
        return new(Group(g1, g2, g3, g4));
    }

    public static ConcreteGroup<Ep4<T1, T2, T3, T4>> Generate<T1, T2, T3, T4>(string name, IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3,
        IGroup<T4> g4)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2> where T3 : struct, IElt<T3> where T4 : struct, IElt<T4>
    {
        return new(name, Group(g1, g2, g3, g4));
    }
}

public readonly struct Gp5<T1, T2, T3, T4, T5> : IGroup<Ep5<T1, T2, T3, T4, T5>> where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
    where T3 : struct, IElt<T3>
    where T4 : struct, IElt<T4>
    where T5 : struct, IElt<T5>
{
    public string Name { get; }
    public IGroup<T1> G1 { get; }
    public IGroup<T2> G2 { get; }
    public IGroup<T3> G3 { get; }
    public IGroup<T4> G4 { get; }
    public IGroup<T5> G5 { get; }

    public Gp5(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3, IGroup<T4> g4, IGroup<T5> g5)
    {
        G1 = g1;
        G2 = g2;
        G3 = g3;
        G4 = g4;
        G5 = g5;
        Hash = (g1.Hash, g2.Hash, g3.Hash, g4.Hash, g5.Hash).GetHashCode();
        Name = $"{G1} x {G2} x {G3} x {G4} x {G5}";
    }

    public bool Equals(IGroup<Ep5<T1, T2, T3, T4, T5>>? other) => other?.Hash == Hash;
    public int Hash { get; }
    public Ep5<T1, T2, T3, T4, T5> Neutral() => new(G1.Neutral(), G2.Neutral(), G3.Neutral(), G4.Neutral(), G5.Neutral());

    public Ep5<T1, T2, T3, T4, T5> Invert(Ep5<T1, T2, T3, T4, T5> e) =>
        new(G1.Invert(e.E1), G2.Invert(e.E2), G3.Invert(e.E3), G4.Invert(e.E4), G5.Invert(e.E5));

    public Ep5<T1, T2, T3, T4, T5> Op(Ep5<T1, T2, T3, T4, T5> e1, Ep5<T1, T2, T3, T4, T5> e2) => new(G1.Op(e1.E1, e2.E1),
        G2.Op(e1.E2, e2.E2), G3.Op(e1.E3, e2.E3), G4.Op(e1.E4, e2.E4), G5.Op(e1.E5, e2.E5));

    public Ep5<T1, T2, T3, T4, T5> this[params ValueType[] us]
    {
        get
        {
            dynamic us0 = new ValueType[5];
            if (us.Length == 1)
            {
                dynamic us1 = us[0];
                us0[0] = us1.Item1;
                us0[1] = us1.Item2;
                us0[2] = us1.Item3;
                us0[3] = us1.Item4;
                us0[4] = us1.Item5;
            }
            else if (us.Length == 5)
                us0 = us;
            else
                throw new GroupException(GroupExceptionType.GroupDef);

            var e1 = G1[us0[0]];
            var e2 = G2[us0[1]];
            var e3 = G3[us0[2]];
            var e4 = G4[us0[3]];
            var e5 = G5[us0[4]];
            return new(e1, e2, e3, e4, e5);
        }
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;

    public IEnumerable<Ep5<T1, T2, T3, T4, T5>> GetGenerators()
    {
        foreach (var e in G1.GetGenerators())
            yield return new Ep5<T1, T2, T3, T4, T5>(e, G2.Neutral(), G3.Neutral(), G4.Neutral(), G5.Neutral());

        foreach (var e in G2.GetGenerators())
            yield return new Ep5<T1, T2, T3, T4, T5>(G1.Neutral(), e, G3.Neutral(), G4.Neutral(), G5.Neutral());

        foreach (var e in G3.GetGenerators())
            yield return new Ep5<T1, T2, T3, T4, T5>(G1.Neutral(), G2.Neutral(), e, G4.Neutral(), G5.Neutral());

        foreach (var e in G4.GetGenerators())
            yield return new Ep5<T1, T2, T3, T4, T5>(G1.Neutral(), G2.Neutral(), G3.Neutral(), e, G5.Neutral());

        foreach (var e in G5.GetGenerators())
            yield return new Ep5<T1, T2, T3, T4, T5>(G1.Neutral(), G2.Neutral(), G3.Neutral(), G4.Neutral(), e);
    }

    public IEnumerable<Ep5<T1, T2, T3, T4, T5>> GetElements()
    {
        foreach (var e1 in G1)
        foreach (var e2 in G2)
        foreach (var e3 in G3)
        foreach (var e4 in G4)
        foreach (var e5 in G5)
            yield return new(e1, e2, e3, e4, e5);
    }

    public IEnumerator<Ep5<T1, T2, T3, T4, T5>> GetEnumerator() => GetElements().GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();
}

public readonly struct Ep5<T1, T2, T3, T4, T5> : IElt<Ep5<T1, T2, T3, T4, T5>> where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
    where T3 : struct, IElt<T3>
    where T4 : struct, IElt<T4>
    where T5 : struct, IElt<T5>
{
    public T1 E1 { get; }
    public T2 E2 { get; }
    public T3 E3 { get; }
    public T4 E4 { get; }
    public T5 E5 { get; }

    public Ep5(T1 e1, T2 e2, T3 e3, T4 e4, T5 e5)
    {
        E1 = e1;
        E2 = e2;
        E3 = e3;
        E4 = e4;
        E5 = e5;
        Hash = (e1.Hash, e2.Hash, e3.Hash, e4.Hash, e5.Hash).GetHashCode();
    }

    public bool Equals(Ep5<T1, T2, T3, T4, T5> other) =>
        other.Hash == Hash && (E1, E2, E3, E4, E5).Equals((other.E1, other.E2, other.E3, other.E4, other.E5));

    public int CompareTo(Ep5<T1, T2, T3, T4, T5> other)
    {
        return (E1, E2, E3, E4, E5).CompareTo((other.E1, other.E2, other.E3, other.E4, other.E5));
    }

    public int Hash { get; }
    public override int GetHashCode() => Hash;
    public override string ToString() => $"({E1}, {E2}, {E3}, {E4}, {E5})";

    public void Deconstruct(out T1 e1, out T2 e2, out T3 e3, out T4 e4, out T5 e5)
    {
        (e1, e2, e3, e4, e5) = (E1, E2, E3, E4, E5);
    }
}

public static partial class Product
{
    public static Gp5<T1, T2, T3, T4, T5> Group<T1, T2, T3, T4, T5>(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3, IGroup<T4> g4,
        IGroup<T5> g5)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
        where T3 : struct, IElt<T3>
        where T4 : struct, IElt<T4>
        where T5 : struct, IElt<T5>
    {
        return new(g1, g2, g3, g4, g5);
    }

    public static Ep5<T1, T2, T3, T4, T5> Elt<T1, T2, T3, T4, T5>(T1 e1, T2 e2, T3 e3, T4 e4, T5 e5)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
        where T3 : struct, IElt<T3>
        where T4 : struct, IElt<T4>
        where T5 : struct, IElt<T5>
    {
        return new(e1, e2, e3, e4, e5);
    }

    public static ConcreteGroup<Ep5<T1, T2, T3, T4, T5>> Generate<T1, T2, T3, T4, T5>(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3,
        IGroup<T4> g4, IGroup<T5> g5)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
        where T3 : struct, IElt<T3>
        where T4 : struct, IElt<T4>
        where T5 : struct, IElt<T5>
    {
        return new(Group(g1, g2, g3, g4, g5));
    }

    public static ConcreteGroup<Ep5<T1, T2, T3, T4, T5>> Generate<T1, T2, T3, T4, T5>(string name, IGroup<T1> g1, IGroup<T2> g2,
        IGroup<T3> g3, IGroup<T4> g4, IGroup<T5> g5)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
        where T3 : struct, IElt<T3>
        where T4 : struct, IElt<T4>
        where T5 : struct, IElt<T5>
    {
        return new(name, Group(g1, g2, g3, g4, g5));
    }
}

public readonly struct Gp6<T1, T2, T3, T4, T5, T6> : IGroup<Ep6<T1, T2, T3, T4, T5, T6>> where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
    where T3 : struct, IElt<T3>
    where T4 : struct, IElt<T4>
    where T5 : struct, IElt<T5>
    where T6 : struct, IElt<T6>
{
    public string Name { get; }
    public IGroup<T1> G1 { get; }
    public IGroup<T2> G2 { get; }
    public IGroup<T3> G3 { get; }
    public IGroup<T4> G4 { get; }
    public IGroup<T5> G5 { get; }
    public IGroup<T6> G6 { get; }

    public Gp6(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3, IGroup<T4> g4, IGroup<T5> g5, IGroup<T6> g6)
    {
        G1 = g1;
        G2 = g2;
        G3 = g3;
        G4 = g4;
        G5 = g5;
        G6 = g6;
        Hash = (g1.Hash, g2.Hash, g3.Hash, g4.Hash, g5.Hash, g6.Hash).GetHashCode();
        Name = $"{G1} x {G2} x {G3} x {G4} x {G5} x {G6}";
    }

    public bool Equals(IGroup<Ep6<T1, T2, T3, T4, T5, T6>>? other) => other?.Hash == Hash;
    public int Hash { get; }

    public Ep6<T1, T2, T3, T4, T5, T6> Neutral() =>
        new(G1.Neutral(), G2.Neutral(), G3.Neutral(), G4.Neutral(), G5.Neutral(), G6.Neutral());

    public Ep6<T1, T2, T3, T4, T5, T6> Invert(Ep6<T1, T2, T3, T4, T5, T6> e) => new(G1.Invert(e.E1), G2.Invert(e.E2), G3.Invert(e.E3),
        G4.Invert(e.E4), G5.Invert(e.E5), G6.Invert(e.E6));

    public Ep6<T1, T2, T3, T4, T5, T6> Op(Ep6<T1, T2, T3, T4, T5, T6> e1, Ep6<T1, T2, T3, T4, T5, T6> e2) => new(G1.Op(e1.E1, e2.E1),
        G2.Op(e1.E2, e2.E2), G3.Op(e1.E3, e2.E3), G4.Op(e1.E4, e2.E4), G5.Op(e1.E5, e2.E5), G6.Op(e1.E6, e2.E6));

    public Ep6<T1, T2, T3, T4, T5, T6> this[params ValueType[] us]
    {
        get
        {
            dynamic us0 = new ValueType[6];
            if (us.Length == 1)
            {
                dynamic us1 = us[0];
                us0[0] = us1.Item1;
                us0[1] = us1.Item2;
                us0[2] = us1.Item3;
                us0[3] = us1.Item4;
                us0[4] = us1.Item5;
                us0[5] = us1.Item6;
            }
            else if (us.Length == 6)
                us0 = us;
            else
                throw new GroupException(GroupExceptionType.GroupDef);

            var e1 = G1[us0[0]];
            var e2 = G2[us0[1]];
            var e3 = G3[us0[2]];
            var e4 = G4[us0[3]];
            var e5 = G5[us0[4]];
            var e6 = G6[us0[5]];
            return new(e1, e2, e3, e4, e5, e6);
        }
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;

    public IEnumerable<Ep6<T1, T2, T3, T4, T5, T6>> GetGenerators()
    {
        foreach (var e in G1.GetGenerators())
            yield return new Ep6<T1, T2, T3, T4, T5, T6>(e, G2.Neutral(), G3.Neutral(), G4.Neutral(), G5.Neutral(), G6.Neutral());

        foreach (var e in G2.GetGenerators())
            yield return new Ep6<T1, T2, T3, T4, T5, T6>(G1.Neutral(), e, G3.Neutral(), G4.Neutral(), G5.Neutral(), G6.Neutral());

        foreach (var e in G3.GetGenerators())
            yield return new Ep6<T1, T2, T3, T4, T5, T6>(G1.Neutral(), G2.Neutral(), e, G4.Neutral(), G5.Neutral(), G6.Neutral());

        foreach (var e in G4.GetGenerators())
            yield return new Ep6<T1, T2, T3, T4, T5, T6>(G1.Neutral(), G2.Neutral(), G3.Neutral(), e, G5.Neutral(), G6.Neutral());

        foreach (var e in G5.GetGenerators())
            yield return new Ep6<T1, T2, T3, T4, T5, T6>(G1.Neutral(), G2.Neutral(), G3.Neutral(), G4.Neutral(), e, G6.Neutral());

        foreach (var e in G6.GetGenerators())
            yield return new Ep6<T1, T2, T3, T4, T5, T6>(G1.Neutral(), G2.Neutral(), G3.Neutral(), G4.Neutral(), G5.Neutral(), e);
    }

    public IEnumerable<Ep6<T1, T2, T3, T4, T5, T6>> GetElements()
    {
        foreach (var e1 in G1)
        foreach (var e2 in G2)
        foreach (var e3 in G3)
        foreach (var e4 in G4)
        foreach (var e5 in G5)
        foreach (var e6 in G6)
            yield return new(e1, e2, e3, e4, e5, e6);
    }

    public IEnumerator<Ep6<T1, T2, T3, T4, T5, T6>> GetEnumerator() => GetElements().GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();
}

public readonly struct Ep6<T1, T2, T3, T4, T5, T6> : IElt<Ep6<T1, T2, T3, T4, T5, T6>> where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
    where T3 : struct, IElt<T3>
    where T4 : struct, IElt<T4>
    where T5 : struct, IElt<T5>
    where T6 : struct, IElt<T6>
{
    public T1 E1 { get; }
    public T2 E2 { get; }
    public T3 E3 { get; }
    public T4 E4 { get; }
    public T5 E5 { get; }
    public T6 E6 { get; }

    public Ep6(T1 e1, T2 e2, T3 e3, T4 e4, T5 e5, T6 e6)
    {
        E1 = e1;
        E2 = e2;
        E3 = e3;
        E4 = e4;
        E5 = e5;
        E6 = e6;
        Hash = (e1.Hash, e2.Hash, e3.Hash, e4.Hash, e5.Hash, e6.Hash).GetHashCode();
    }

    public bool Equals(Ep6<T1, T2, T3, T4, T5, T6> other) => other.Hash == Hash &&
                                                             (E1, E2, E3, E4, E5, E6).Equals((other.E1, other.E2, other.E3, other.E4,
                                                                 other.E5, other.E6));

    public int CompareTo(Ep6<T1, T2, T3, T4, T5, T6> other)
    {
        return (E1, E2, E3, E4, E5, E6).CompareTo((other.E1, other.E2, other.E3, other.E4, other.E5, other.E6));
    }

    public int Hash { get; }
    public override int GetHashCode() => Hash;
    public override string ToString() => $"({E1}, {E2}, {E3}, {E4}, {E5}, {E6})";

    public void Deconstruct(out T1 e1, out T2 e2, out T3 e3, out T4 e4, out T5 e5, out T6 e6)
    {
        (e1, e2, e3, e4, e5, e6) = (E1, E2, E3, E4, E5, E6);
    }
}

public static partial class Product
{
    public static Gp6<T1, T2, T3, T4, T5, T6> Group<T1, T2, T3, T4, T5, T6>(IGroup<T1> g1, IGroup<T2> g2, IGroup<T3> g3, IGroup<T4> g4,
        IGroup<T5> g5, IGroup<T6> g6)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
        where T3 : struct, IElt<T3>
        where T4 : struct, IElt<T4>
        where T5 : struct, IElt<T5>
        where T6 : struct, IElt<T6>
    {
        return new(g1, g2, g3, g4, g5, g6);
    }

    public static Ep6<T1, T2, T3, T4, T5, T6> Elt<T1, T2, T3, T4, T5, T6>(T1 e1, T2 e2, T3 e3, T4 e4, T5 e5, T6 e6)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
        where T3 : struct, IElt<T3>
        where T4 : struct, IElt<T4>
        where T5 : struct, IElt<T5>
        where T6 : struct, IElt<T6>
    {
        return new(e1, e2, e3, e4, e5, e6);
    }

    public static ConcreteGroup<Ep6<T1, T2, T3, T4, T5, T6>> Generate<T1, T2, T3, T4, T5, T6>(IGroup<T1> g1, IGroup<T2> g2,
        IGroup<T3> g3, IGroup<T4> g4, IGroup<T5> g5, IGroup<T6> g6)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
        where T3 : struct, IElt<T3>
        where T4 : struct, IElt<T4>
        where T5 : struct, IElt<T5>
        where T6 : struct, IElt<T6>
    {
        return new(Group(g1, g2, g3, g4, g5, g6));
    }

    public static ConcreteGroup<Ep6<T1, T2, T3, T4, T5, T6>> Generate<T1, T2, T3, T4, T5, T6>(string name, IGroup<T1> g1,
        IGroup<T2> g2, IGroup<T3> g3, IGroup<T4> g4, IGroup<T5> g5, IGroup<T6> g6)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
        where T3 : struct, IElt<T3>
        where T4 : struct, IElt<T4>
        where T5 : struct, IElt<T5>
        where T6 : struct, IElt<T6>
    {
        return new(name, Group(g1, g2, g3, g4, g5, g6));
    }
}