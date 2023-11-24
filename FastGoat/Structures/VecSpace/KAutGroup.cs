using System.Collections;
using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

public readonly struct KAutGroup<K> : IGroup<KAut<K>> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public KPoly<K> F { get; }
    private EPoly<K> X { get; }

    public KAutGroup(KPoly<K> P)
    {
        F = P;
        Hash = P.Hash;
        Name = $"Q[{P.x}]/({P})";
        X = new(F);
    }

    public IEnumerator<KAut<K>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();

    public bool Equals(IGroup<KAut<K>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public KAut<K> KAut(EPoly<K> e) => new(this, e);

    public KAut<K> this[params ValueType[] us]
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
                var e0 = new EPoly<K>(F, new KPoly<K>(F.x, F.KZero, coefs));
                return new(this, e0);
            }

            var polys = cycles.Select(c => c.Table.Select(e => e * one).ToArray())
                .Select(coefs => new EPoly<K>(f, new KPoly<K>(f.x, f.KZero, coefs))).ToArray();

            var e1 = polys.Aggregate(X, (acc, f0) => acc.Substitute(f0));
            return new(this, e1);
        }
    }

    public IEnumerable<KAut<K>> GetElements()
    {
        yield return Neutral();
    }

    public IEnumerable<KAut<K>> GetGenerators()
    {
        yield return new(this, X);
    }

    public KAut<K> Neutral() => new(this, X);

    public KAut<K> Invert(KAut<K> e)
    {
        var n = e.E.X;
        var tmp0 = e.E.Clone;
        EPoly<K> tmp1;
        do
        {
            tmp1 = tmp0.Clone;
            tmp0 = e.E.Substitute(tmp1);
        } while (!tmp0.Equals(n));

        return new(this, tmp1);
    }

    public KAut<K> Op(KAut<K> e1, KAut<K> e2) => new(this, e1.E.Substitute(e2.E));
    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}