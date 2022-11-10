using System.Collections;

namespace FastGoat.Structures.GenericGroup;

public class AutomorphismGroup<T> : IGroup<Automorphism<T>> where T : struct, IElt<T>
{
    public ConcreteGroup<T> G { get; }

    public AutomorphismGroup(ConcreteGroup<T> g)
    {
        G = g;
        Name = $"Aut({g.Name})";
        Hash = (g.Hash, "Aut").GetHashCode();
    }

    public IEnumerable<Automorphism<T>> GetGenerators()
    {
        throw new GroupException(GroupExceptionType.GroupDef);
    }

    public IEnumerable<Automorphism<T>> GetElements()
    {
        yield return Neutral();
    }

    public IEnumerator<Automorphism<T>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();

    public bool Equals(IGroup<Automorphism<T>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public Automorphism<T> Neutral() => new Automorphism<T>(this);

    public Automorphism<T> Invert(Automorphism<T> e)
    {
        if (!Equals(e.AutGroup))
            throw new GroupException(GroupExceptionType.BaseGroup);

        var autMap = e.AutMap.ToDictionary(kp => kp.Value, kp => kp.Key);
        return new Automorphism<T>(this, autMap);
    }

    public Automorphism<T> Op(Automorphism<T> e1, Automorphism<T> e2)
    {
        if (!Equals(e1.AutGroup) || !Equals(e2.AutGroup))
            throw new GroupException(GroupExceptionType.BaseGroup);

        var autMap = e1.AutMap.ToDictionary(kp => kp.Key, kp => e2.AutMap[kp.Value]);
        return new Automorphism<T>(this, autMap);
    }

    public Automorphism<T> Create(IReadOnlyDictionary<T, T> autMap)
    {
        if (autMap.Count != G.Count())
            throw new GroupException(GroupExceptionType.GroupDef);

        return new Automorphism<T>(this, autMap);
    }

    public Automorphism<T> Create(Homomorphism<T, T> hom) => Create(hom.HomMap);

    public Automorphism<T> this[params ValueType[] us]
    {
        get
        {
            var pMap = us.Select(e => (dynamic)e).ToDictionary(kp => (T)kp.Item1, kp => (T)kp.Item2);
            var autMap = Group.AutomorphismMap(G, pMap);
            if (autMap.Count != G.Count())
                throw new GroupException(GroupExceptionType.GroupDef);

            return new Automorphism<T>(this, autMap);
        }
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}