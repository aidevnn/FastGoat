using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public readonly struct CyclotomicGroupBase<K> : IGroup<EPoly<K>> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public int N { get; }
    private EPoly<K> Zeta { get; }
    private KPoly<K> Poly { get; }
    public KPoly<EPoly<K>> X => FG.KPoly('X', Zeta);

    public CyclotomicGroupBase(K scalar, int n, string name)
    {
        N = n;
        var p = scalar.P;
        if (p > 2)
        {
            if (n % p == 0)
                throw new ArgumentException($"In {name}, P={p} divide N={n}");
        }

        var nm = name;
        Name = $"U{n}({nm})";
        var poly = FG.CyclotomicPolynomial(n);
        if (scalar.P == 0)
        {
            Poly = new KPoly<K>('ζ', scalar, poly.Coefs.Select(c => ((int)c.Num) * scalar.One).ToArray());
            Hash = ("CF", p, Poly.Hash).GetHashCode();

            var q = n == 2 ? 1 : IntExt.UnInvertible(n).First(e => e.Key != 1).Key;
            var x = FG.EPoly(Poly);
            Zeta = x.Pow(q);
        }
        else if (scalar is EPoly<ZnInt> s0)
        {
            var poly0 = new KPoly<EPoly<ZnInt>>(poly.x, s0.Zero, poly.Coefs.Select(c => (int)c.Num * s0.One).ToArray());
            var facts = IntFactorisation.Firr(poly0, s0.X).Order().ToArray();
            Poly = new KPoly<K>('ζ', scalar.One, facts.First().Coefs.Cast<K>().ToArray());
            Hash = ("CF", p, Poly.Hash).GetHashCode();

            var q = n == 2 ? 1 : IntExt.UnInvertible(n).First(e => e.Key != 1).Key;
            var x = FG.EPoly(Poly);
            Zeta = x.Pow(q);
        }
        else
            throw new();
    }

    public IEnumerator<EPoly<K>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<EPoly<K>>? other) => other is not null && other.Neutral().F.Equals(Poly);

    public int Hash { get; }
    public string Name { get; }

    public EPoly<K> this[params ValueType[] us]
    {
        get
        {
            var u0 = Convert.ToInt32(us[0]);
            return Zeta.Pow(u0);
        }
    }

    public IEnumerable<EPoly<K>> GetElements()
    {
        for (int i = 0; i < N; i++)
            yield return Zeta.Pow(i);
    }

    public IEnumerable<EPoly<K>> GetGenerators()
    {
        yield return Zeta;
    }

    public EPoly<K> Neutral() => Zeta.One;

    public EPoly<K> Invert(EPoly<K> e) => 1 / e;

    public EPoly<K> Op(EPoly<K> e1, EPoly<K> e2) => e1 * e2;
    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}

public class NthRootQ : ConcreteGroup<EPoly<Rational>>
{
    public NthRootQ(int n) : this(new CyclotomicGroupBase<Rational>(Rational.KZero(), n, "Q"))
    {
    }

    private NthRootQ(CyclotomicGroupBase<Rational> cg) : base(cg)
    {
        CG = cg;
    }

    public CyclotomicGroupBase<Rational> CG { get; }
    public int N => CG.N;
    public KPoly<EPoly<Rational>> X => CG.X;
    public EPoly<Rational>[] PrimitivesRoots() => ElementsOrders.Where(e => e.Value == N).Select(e => e.Key).Order().ToArray();
}

public class NthRootFq : ConcreteGroup<EPoly<EPoly<ZnInt>>>
{
    public NthRootFq(int n, int q) : this(new CyclotomicGroupBase<EPoly<ZnInt>>(FG.FqX(q, 'α'), n, $"F{q}"))
    {
    }

    private NthRootFq(CyclotomicGroupBase<EPoly<ZnInt>> cg) : base(cg)
    {
        CG = cg;
    }

    public CyclotomicGroupBase<EPoly<ZnInt>> CG { get; }
    public int N => CG.N;
    public KPoly<EPoly<EPoly<ZnInt>>> X => CG.X;
    public EPoly<EPoly<ZnInt>>[] PrimitivesRoots() => ElementsOrders.Where(e => e.Value == N).Select(e => e.Key).Order().ToArray();
}