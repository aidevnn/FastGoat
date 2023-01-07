using System.Collections;
using FastGoat.Structures;

namespace FastGoat.Structures.VecSpace;

public readonly struct GLn<K> : IGroup<KMatrix<K>> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    private KMatrix<K>[] Generators { get; }
    private KMatrix<K> idN { get; }
    public int N => idN.N;

    public GLn(KMatrix<K> e0, params KMatrix<K>[] others)
    {
        Generators = others.Prepend(e0).ToArray();
        idN = e0.One;
        var one = idN.KOne;
        if (Generators.Select(m => m.Det).Any(d => !d.Equals(one) || !d.Equals(-one)))
            throw new Exception($"|Det| != 1");

        Hash = Generators.Aggregate(0, (acc, m) => (acc, m.Hash).GetHashCode());
        Name = $"GL({idN.N}, {typeof(K)})";
    }

    public GLn(string name, KMatrix<K> e0, params KMatrix<K>[] others)
    {
        Generators = others.Prepend(e0).ToArray();
        idN = e0.One;
        var one = idN.KOne;
        if (Generators.Select(m => m.Det).Any(d => !d.Equals(one) || !d.Equals(-one)))
            throw new Exception($"|Det| != 1");

        Hash = Generators.Aggregate(0, (acc, m) => (acc, m.Hash).GetHashCode());
        Name = name;
    }

    public GLn(int n, K scalar)
    {
        Generators = Array.Empty<KMatrix<K>>();
        idN = new KMatrix<K>(scalar, n, n).One;
        Name = $"GL({idN.N}, {scalar})";
        Hash = Name.GetHashCode();
    }

    public GLn(string name, int n, K scalar)
    {
        Generators = Array.Empty<KMatrix<K>>();
        idN = new KMatrix<K>(scalar, n, n).One;
        Name = $"GL({idN.N}, {name})";
        Hash = Name.GetHashCode();
    }

    public IEnumerator<KMatrix<K>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<KMatrix<K>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public KMatrix<K> this[params ValueType[] us]
    {
        get
        {
            if (us.Any(e => e is not int && e is not K))
                throw new GroupException(GroupExceptionType.GroupDef);

            var kzero = idN.KZero;
            var kone = idN.KOne;
            var us0 = us.Select(e => e is int e0 ? e0 * kone : e is K e1 ? e1 : kzero).ToArray();
            var m = new KMatrix<K>(Ring.Matrix(N, us0));
            var det = m.Det;
            if (!det.Equals(kone) && !det.Equals(-kone))
                throw new GroupException(GroupExceptionType.GroupDef);

            return m;
        }
    }

    public IEnumerable<KMatrix<K>> GetElements()
    {
        yield return idN;
    }

    public IEnumerable<KMatrix<K>> GetGenerators() => Generators;

    public KMatrix<K> Neutral() => idN;

    public KMatrix<K> Invert(KMatrix<K> e) => e.Inv();

    public KMatrix<K> Op(KMatrix<K> e1, KMatrix<K> e2) => e1 * e2;

    public override string ToString() => Name;
    public override int GetHashCode() => Hash;
}