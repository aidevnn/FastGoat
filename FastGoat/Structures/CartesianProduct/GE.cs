using System.Collections;
using FastGoat.Commons;

namespace FastGoat.Structures.CartesianProduct;

public struct Gp<T> : IGroup<Ep<T>> where T : struct, IElt<T>
{
    public string Name { get; }
    public IGroup<T>[] Gi { get; }

    public Gp(IGroup<T>[] gi)
    {
        Gi = gi.ToArray();
        HashCode code = new HashCode();
        foreach (var g in gi)
        {
            code.Add(g);
        }

        Hash = code.ToHashCode();
        Name = gi.Glue(" x ");
    }

    public bool Equals(IGroup<Ep<T>>? other) => other?.Hash == Hash;
    public int Hash { get; }
    public Ep<T> Neutral() => new(Gi.Select(g => g.Neutral()).ToArray());
    public Ep<T> Invert(Ep<T> e) => new(Gi.Select((g, i) => g.Invert(e.Ei[i])).ToArray());
    public Ep<T> Op(Ep<T> e1, Ep<T> e2) => new(Gi.Select((g, i) => g.Op(e1.Ei[i], e2.Ei[i])).ToArray());

    public Ep<T> this[params ValueType[] us]
    {
        get
        {
            dynamic us0 = new ValueType[Gi.Length];
            if (us.Length == 1)
            {
                dynamic us1 = us[0];
                for (int i = 0; i < Gi.Length; ++i)
                    us0[i] = us1[i];
            }
            else if (us.Length == Gi.Length)
                us0 = us;
            else
                throw new GroupException(GroupExceptionType.GroupDef);

            var ei = Gi.Select((g, i) => g[us[i]]).ToArray();
            return new(ei);
        }
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;

    public IEnumerable<Ep<T>> GetGenerators()
    {
        for (int i = 0; i < Gi.Length; ++i)
        {
            var g = Gi[i];
            foreach (var ei in g.GetGenerators())
            {
                var ep = Gi.Select((g0, i0) => i != i0 ? g0.Neutral() : ei).ToArray();
                yield return new Ep<T>(ep);
            }
        }
    }

    public IEnumerable<Ep<T>> GetElements()
    {
        foreach (var ep in Gi.Select(g => g.GetElements()).MultiLoop())
            yield return new(ep.ToArray());
    }

    public IEnumerator<Ep<T>> GetEnumerator() => GetElements().GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();
}

public struct Ep<T> : IElt<Ep<T>> where T : IElt<T>
{
    public T[] Ei { get; }

    public Ep(T[] ei)
    {
        Ei = ei.ToArray();

        HashCode code = new HashCode();
        foreach (var e in ei)
        {
            code.Add(e);
        }

        Hash = code.ToHashCode();
    }

    public bool Equals(Ep<T> other) => other.Hash == Hash;

    public int CompareTo(Ep<T> other)
    {
        if (Ei.Length != other.Ei.Length)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return Ei.SequenceCompareTo(other.Ei);
    }

    public int Hash { get; }
    public override int GetHashCode() => Hash;
    public override string ToString() => $"({Ei.Glue(", ")})";
}

public static partial class Product
{
    public static Gp<T> Gp<T>(params IGroup<T>[] gn) where T : struct, IElt<T>
    {
        return new(gn);
    }

    public static Ep<T> Ep<T>(params T[] en) where T : struct, IElt<T>
    {
        return new(en);
    }

    public static ConcreteGroup<Ep<T>> GpGenerate<T>(params IGroup<T>[] gn) where T : struct, IElt<T>
    {
        return new(Gp(gn));
    }
}