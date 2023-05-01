using System.Collections;
using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

public readonly struct KAut<K> : IGroup<EPoly<K>> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public KPoly<K> F { get; }
    private EPoly<K> X { get; }
    public KAut(KPoly<K> P)
    {
        F = P;
        Hash = P.Hash;
        Name = $"Q[{P.x}]/({P})";
        X = new(F);
    }
    
    public IEnumerator<EPoly<K>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();

    public bool Equals(IGroup<EPoly<K>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public EPoly<K> this[params ValueType[] us]
    {
        get
        {
            var cycles = us.Select(u => (Tuple2Array)u).ToArray();
            if (cycles.Any(c => c.Table.Length == 0))
                throw new GroupException(GroupExceptionType.GroupDef);

            var f = F;
            var one = F.KOne;
            if (cycles.All(c => c.Table.Length == 1) && cycles.Length == F.Degree + 1)
            {
                var coefs = cycles.SelectMany(c => c.Table).Select(c => c * one).ToArray();
                return new(F, new KPoly<K>(F.x, F.KZero, coefs));
            }

            var polys = cycles.Select(c => c.Table.Select(e => e * one).ToArray())
                .Select(coefs => new EPoly<K>(f, new KPoly<K>(f.x, f.KZero, coefs))).ToArray();

            return polys.Aggregate(X, (acc, f0) => acc.Substitute(f0));
        }
    }

    public IEnumerable<EPoly<K>> GetElements()
    {
        yield return Neutral();
    }

    public IEnumerable<EPoly<K>> GetGenerators()
    {
        yield return X;
    }

    public EPoly<K> Neutral() => X.X;

    public EPoly<K> Invert(EPoly<K> e) => e.Inv();

    public EPoly<K> Op(EPoly<K> e1, EPoly<K> e2) => e1.Substitute(e2);
    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}