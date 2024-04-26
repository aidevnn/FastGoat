using System.Collections;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Structures.VecSpace;

public readonly struct GLn<K> : IGroup<KMatrix<K>> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    private KMatrix<K>[] Generators { get; }
    private KMatrix<K> idN { get; }
    public int N => idN.N;
    public K KZero => idN.KZero;

    public GLn(KMatrix<K> e0, params KMatrix<K>[] others)
    {
        Generators = others.Prepend(e0).ToArray();
        idN = e0.One;
        var zero = idN.KZero;
        if (Generators.Select(m => m.Det).Any(d => d.Equals(zero)))
            throw new Exception($"|Det| == 0");

        Hash = Generators.Aggregate(0, (acc, m) => (acc, m.Hash).GetHashCode());
        Name = $"GL({idN.N}, {typeof(K)})";
    }

    public GLn(string name, KMatrix<K> e0, params KMatrix<K>[] others)
    {
        Generators = others.Prepend(e0).ToArray();
        idN = e0.One;
        var zero = idN.KZero;
        if (Generators.Select(m => m.Det).Any(d => d.Equals(zero)))
            throw new Exception($"|Det| == 0");

        Hash = Generators.Aggregate(0, (acc, m) => (acc, m.Hash).GetHashCode());
        Name = name;
    }

    public GLn(int n, K scalar)
    {
        idN = new KMatrix<K>(scalar, n, n).One;
        Generators = [idN];
        Name = $"GL({idN.N}, {scalar})";
        Hash = Name.GetHashCode();
    }

    public GLn(string name, int n, K scalar)
    {
        idN = new KMatrix<K>(scalar, n, n).One;
        Generators = [idN];
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

    public static bool Contains(IGroup<KMatrix<K>> g, KMatrix<K> m)
    {
        var g0 = g is GLn<K> ? (GLn<K>)g : default;
        var kone = g0.idN.KOne;
        var det = m.Det;

        Console.WriteLine($"GLn<K>");
        return det.Equals(kone) || det.Equals(-kone);
    }

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

            if (K.IsValuedField)
            {
                if (Math.Abs(K.Abs(det)) < 1e-7)
                    throw new GroupException(GroupExceptionType.GroupDef);
            }
            else if (det.Equals(kzero))
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