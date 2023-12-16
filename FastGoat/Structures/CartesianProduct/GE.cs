using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.CartesianProduct;

public readonly struct Gp<T> : IGroup<Ep<T>> where T : struct, IElt<T>
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
        Name = gi.Length == 0 ? "C1" : gi.Glue(" x ");
    }

    public bool Equals(IGroup<Ep<T>>? other) => other?.Hash == Hash;
    public int Hash { get; }
    public Ep<T> Neutral() => new(Gi.Select(g => g.Neutral()).ToArray());
    public Ep<T> Invert(Ep<T> e) => new(Gi.Select((g, i) => g.Invert(e.Ei[i])).ToArray());
    public Ep<T> Op(Ep<T> e1, Ep<T> e2) => new(Gi.Select((g, i) => g.Op(e1.Ei[i], e2.Ei[i])).ToArray());
    public Ep<T> Act(T t, Ep<T> e) => new(Gi.Select((g, i) => g.Op(t, e.Ei[i])).ToArray());

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

    public Ep()
    {
        Ei = new T[] { };
        Hash = 0;
    }

    public Ep(params T[] ei)
    {
        Ei = ei;
        Hash = ei.Aggregate(0, (acc, e) => (acc, e).GetHashCode());
    }

    public bool Equals(Ep<T> other) => Ei.SequenceEqual(other.Ei);

    public int CompareTo(Ep<T> other)
    {
        if (Ei.Length != other.Ei.Length)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return Ei.SequenceCompareTo(other.Ei);
    }

    public T this[int index] => Ei[index];
    public Ep<T> SkipAt(int i) => new(Ei.SkipAt(i).ToArray());

    public int Hash { get; }
    public override int GetHashCode() => Hash;
    public override string ToString() => Ei.Length == 1 ? $"{Ei[0]}" : $"({Ei.Glue(", ")})";
}

public static partial class Product
{
    public static Gp<T> Gp<T>(params IGroup<T>[] gn) where T : struct, IElt<T>
    {
        return new(gn);
    }

    public static Gp<T> Gp<T>(IGroup<T> gn, int n) where T : struct, IElt<T>
    {
        return new(Enumerable.Repeat(gn, n).ToArray());
    }

    public static Ep<T> Ep<T>(params T[] en) where T : struct, IElt<T>
    {
        return new(en);
    }

    public static ConcreteGroup<Ep<T>> GpGenerate<T>(string name, params IGroup<T>[] gn) where T : struct, IElt<T>
    {
        return new(name, Gp(gn));
    }

    public static ConcreteGroup<Ep<T>> GpGenerate<T>(params IGroup<T>[] gn) where T : struct, IElt<T>
    {
        return new(Gp(gn));
    }

    public static ConcreteGroup<Ep<T>> GpGenerate<T>(string name, IGroup<T> g, int n) where T : struct, IElt<T>
    {
        return new(name, Gp(Enumerable.Repeat(g, n).ToArray()));
    }

    public static ConcreteGroup<Ep<T>> GpGenerate<T>(IGroup<T> g, int n) where T : struct, IElt<T>
    {
        return new(Gp(Enumerable.Repeat(g, n).ToArray()));
    }
}