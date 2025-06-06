using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;

namespace FastGoat.UserGroup.EllCurve;

public readonly struct EllPoly<K> : IElt<EllPoly<K>>, IRingElt<EllPoly<K>>, IModuleElt<K, EllPoly<K>>,
    IVsElt<K, EllPoly<K>>, IFieldElt<EllPoly<K>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public SortedList<TriVar, K> Coefs { get; }
    public IndTriVar IndTriVar { get; }

    public EllPoly(K scalar, MonomOrder monomOrder = MonomOrder.Lex)
    {
        IndTriVar = new(monomOrder);
        Coefs = new(IndTriVar.GetComparer()) { [new()] = scalar };
        KZero = scalar.Zero;
        Hash = Coefs.Aggregate(IndTriVar.GetHashCode(), (acc, a) => (acc, (a.Key, a.Value.Hash)).GetHashCode());
    }

    public EllPoly(IndTriVar triVar, K scalar)
    {
        IndTriVar = triVar;
        Coefs = new(IndTriVar.GetComparer()) { [new()] = scalar };
        KZero = scalar.Zero;
        Hash = Coefs.Aggregate(IndTriVar.GetHashCode(), (acc, a) => (acc, (a.Key, a.Value.Hash)).GetHashCode());
    }

    public EllPoly(IndTriVar triVar, K scalar, IDictionary<TriVar, K> coefs)
    {
        IndTriVar = triVar;
        Coefs = new(IndTriVar.GetComparer());
        foreach (var (k, v) in coefs.Where(e => !e.Value.IsZero()))
            Coefs[k] = v;

        if (Coefs.Count == 0)
            Coefs[new()] = scalar.Zero;

        KZero = scalar.Zero;
        Hash = Coefs.Aggregate(IndTriVar.GetHashCode(), (acc, a) => (acc, (a.Key, a.Value.Hash)).GetHashCode());
    }

    public EllPoly(IndTriVar triVar, (TriVar xyz, K c) coef) : this(triVar, coef.c,
        new Dictionary<TriVar, K>() { [coef.xyz] = coef.c })
    {
    }

    public bool Equals(EllPoly<K> other)
    {
        return Coefs.Count == other.Coefs.Count &&
               Coefs.All(e => other.Coefs.ContainsKey(e.Key) && other.Coefs[e.Key].Equals(e.Value));
    }

    public int CompareTo(EllPoly<K> other)
    {
        var coefs = Coefs.OrderByDescending(kp => kp.Key, IndTriVar.GetComparer()).Select(kp => (kp.Key, kp.Value));
        var oCoefs = other.Coefs.OrderByDescending(kp => kp.Key, other.IndTriVar.GetComparer())
            .Select(kp => (kp.Key, kp.Value));
        return coefs.SequenceCompareTo(oCoefs);
    }

    public int Hash { get; }

    public bool IsZero()
    {
        foreach (var (key, value) in Coefs)
        {
            if (!key.IsOne() || !value.IsZero())
                return false;
        }

        return true;
    }

    public int P => KZero.P;
    public K KZero { get; }
    public K KOne => KZero.One;
    public EllPoly<K> Zero => new(IndTriVar, KZero);
    public EllPoly<K> One => new(IndTriVar, KOne);

    public K ConstTerm => this[new()];

    public (K lc, TriVar lm, EllPoly<K> lt) LeadingDetails
    {
        get
        {
            if (IsZero())
                throw new ArgumentException();

            var (lm, lc) = Coefs.Last();
            return (lc, lm, new(IndTriVar, (lm, lc)));
        }
    }

    public EllPoly<K> X3 => new(IndTriVar, ((1, 0, 0), KOne));
    public EllPoly<K> X2 => new(IndTriVar, ((0, 1, 0), KOne));
    public EllPoly<K> X1 => new(IndTriVar, ((0, 0, 1), KOne));

    public int DegreeOfX3 => Coefs.Max(e => e.Key.X3);
    public int DegreeOfX2 => Coefs.Max(e => e.Key.X2);
    public int DegreeOfX1 => Coefs.Max(e => e.Key.X1);

    public KPoly<K> ToKPolyX3()
    {
        var zero = KZero;
        var arr = (DegreeOfX3 + 1).SeqLazy().Select(_ => zero).ToArray();
        foreach (var i in (DegreeOfX3 + 1).SeqLazy())
            arr[i] = this[(i, 0, 0)];

        return arr.ToKPoly(IndTriVar.X3[0]);
    }

    public KPoly<K> ToKPolyX2()
    {
        var zero = KZero;
        var arr = (DegreeOfX2 + 1).SeqLazy().Select(_ => zero).ToArray();
        foreach (var i in (DegreeOfX2 + 1).SeqLazy())
            arr[i] = this[(0, i, 0)];

        return arr.ToKPoly(IndTriVar.X2[0]);
    }

    public KPoly<K> ToKPolyX1()
    {
        var zero = KZero;
        var arr = (DegreeOfX1 + 1).SeqLazy().Select(_ => zero).ToArray();
        foreach (var i in (DegreeOfX1 + 1).SeqLazy())
            arr[i] = this[(0, 0, i)];

        return arr.ToKPoly(IndTriVar.X1[0]);
    }

    public int Degree => Coefs.Last().Key.Degree;

    public EllPoly<K> Add(EllPoly<K> e)
    {
        var coefs = Coefs.ToDictionary(kv => kv.Key, kv => kv.Value);
        foreach (var (xyz, c) in e.Coefs)
        {
            if (coefs.ContainsKey(xyz))
                coefs[xyz] = coefs[xyz] + c;
            else
                coefs[xyz] = c;
        }

        return new(IndTriVar, KZero, coefs);
    }

    public EllPoly<K> Sub(EllPoly<K> e)
    {
        var coefs = Coefs.ToDictionary(kv => kv.Key, kv => kv.Value);
        foreach (var (xyz, c) in e.Coefs)
        {
            if (coefs.ContainsKey(xyz))
                coefs[xyz] = coefs[xyz] - c;
            else
                coefs[xyz] = -c;
        }

        return new(IndTriVar, KZero, coefs);
    }

    public EllPoly<K> Opp()
    {
        var coefs = Coefs.ToDictionary(kv => kv.Key, kv => -kv.Value);
        return new(IndTriVar, KZero, coefs);
    }

    public EllPoly<K> Mul(EllPoly<K> e)
    {
        var coefs = new Dictionary<TriVar, K>(Coefs.Count + e.Coefs.Count);
        foreach (var (a1, c1) in Coefs)
        {
            foreach (var (a2, c2) in e.Coefs)
            {
                var a3 = a1.Mul(a2);
                var c3 = c1 * c2;
                if (coefs.ContainsKey(a3))
                    coefs[a3] = coefs[a3] + c3;
                else
                    coefs[a3] = c3;
            }
        }

        return new(IndTriVar, KZero, coefs);
    }

    private void InPlaceSubMul(SortedList<TriVar, K> A, SortedList<TriVar, K> B, K c, TriVar m)
    {
        foreach (var (k0, v) in B)
        {
            var k1 = k0.Mul(m);
            if (A.ContainsKey(k1))
                A[k1] = A[k1].Sub(v * c);
            else
                A[k1] = -v * c;
        }

        foreach (var (k, _) in A.Where(e => e.Value.IsZero()).ToArray())
            A.Remove(k);

        if (A.Count == 0)
            A[new()] = c.Zero;
    }

    public (EllPoly<K> quo, EllPoly<K> rem) Div(EllPoly<K> e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        var rem = new SortedList<TriVar, K>(Coefs, IndTriVar.GetComparer());
        var quo = new SortedList<TriVar, K>(IndTriVar.GetComparer());
        var (elm, elc) = e.Coefs.Last();
        var k = Coefs.Count - 1;
        var one = new TriVar();
        while (k >= 0)
        {
            var (alm, alc) = (rem.Keys[k], rem.Values[k]);
            var (mnm0, mnm1) = TriVar.Reduce(alm, elm);
            if (!mnm0.IsOne())
            {
                --k;
                continue;
            }

            var (q, r) = alc.Div(elc);
            if (!r.IsZero())
                break;

            InPlaceSubMul(rem, e.Coefs, q, mnm1);
            k = rem.Count - 1;
            quo[mnm1] = q;
            if (rem.Count == 1 && rem.ContainsKey(one) && rem[one].IsZero())
                break;
        }

        return (new(IndTriVar, KZero, quo), new(IndTriVar, KZero, rem));
    }

    public EllPoly<K> Mul(int k)
    {
        var coefs = Coefs.ToDictionary(kv => kv.Key, kv => k * kv.Value);
        return new(IndTriVar, KZero, coefs);
    }

    public EllPoly<K> Pow(int k)
    {
        if (k < 0)
            throw new GroupException(GroupExceptionType.GroupDef);

        if (k == 0)
            return One;

        return this.FastPow(k);
    }

    public EllPoly<K> KMul(K k)
    {
        var coefs = Coefs.ToDictionary(kv => kv.Key, kv => k * kv.Value);
        return new(IndTriVar, KZero, coefs);
    }

    public EllPoly<K> Inv()
    {
        if (Invertible())
            return new(IndTriVar, ConstTerm.Inv());

        throw new DivideByZeroException();
    }

    public bool Invertible() => Degree == 0 && ConstTerm.Invertible();

    public Dictionary<EllPoly<K>, EllPoly<K>> DecomposeX3()
    {
        var (ind, kzero, x3) = (IndTriVar, KZero, X3);
        return Coefs.GroupBy(e => e.Key.X3).ToDictionary(
            e => x3.Pow(e.Key),
            e => new EllPoly<K>(ind, kzero, e.ToDictionary(a => a.Key.GetX2X1(), a => a.Value))
        );
    }
    
    public Dictionary<EllPoly<K>, EllPoly<K>> DecomposeX2()
    {
        var (ind, kzero, x2) = (IndTriVar, KZero, X2);
        return Coefs.GroupBy(e => e.Key.X2).ToDictionary(
            e => x2.Pow(e.Key),
            e => new EllPoly<K>(ind, kzero, e.ToDictionary(a => a.Key.GetX3X1(), a => a.Value))
        );
    }

    public Dictionary<EllPoly<K>, EllPoly<K>> DecomposeX1()
    {
        var (ind, kzero, x1) = (IndTriVar, KZero, X1);
        return Coefs.GroupBy(e => e.Key.X1).ToDictionary(
            e => x1.Pow(e.Key),
            e => new EllPoly<K>(ind, kzero, e.ToDictionary(a => a.Key.GetX3X2(), a => a.Value))
        );
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var (x3, x2, x1) = Ring.Polynomial(KZero, IndTriVar.MonomOrder, IndTriVar.X3, IndTriVar.X2, IndTriVar.X1)
            .Deconstruct();
        var pol = Coefs.Select(e => e.Value * x3.Pow(e.Key.X3) * x2.Pow(e.Key.X2) * x1.Pow(e.Key.X1))
            .Aggregate(x3.Zero, (acc, e) => acc + e);
        return $"{pol}";
    }

    public K this[TriVar index]
    {
        get
        {
            if (Coefs.TryGetValue(index, out K v))
                return v;

            return KZero;
        }
    }

    public static EllPoly<K> operator +(EllPoly<K> a, EllPoly<K> b) => a.Add(b);

    public static EllPoly<K> operator +(int a, EllPoly<K> b) => b.One.Mul(a).Add(b);

    public static EllPoly<K> operator +(EllPoly<K> a, int b) => a.Add(a.One.Mul(b));

    public static EllPoly<K> operator -(EllPoly<K> a) => a.Opp();

    public static EllPoly<K> operator -(EllPoly<K> a, EllPoly<K> b) => a.Sub(b);

    public static EllPoly<K> operator -(int a, EllPoly<K> b) => b.One.Mul(a).Sub(b);

    public static EllPoly<K> operator -(EllPoly<K> a, int b) => a.Sub(a.One.Mul(b));

    public static EllPoly<K> operator *(EllPoly<K> a, EllPoly<K> b) => a.Mul(b);

    public static EllPoly<K> operator *(int a, EllPoly<K> b) => b.Mul(a);

    public static EllPoly<K> operator *(EllPoly<K> a, int b) => a.Mul(b);

    public static EllPoly<K> operator /(EllPoly<K> a, EllPoly<K> b) => a.Div(b).quo;

    public static EllPoly<K> operator /(EllPoly<K> a, int b) => a.Div(a.One.Mul(b)).quo;

    public static EllPoly<K> operator +(EllPoly<K> a, K b) => a.Add(a.One.KMul(b));

    public static EllPoly<K> operator +(K a, EllPoly<K> b) => b.One.KMul(a).Add(b);

    public static EllPoly<K> operator -(EllPoly<K> a, K b) => a.Sub(a.One.KMul(b));

    public static EllPoly<K> operator -(K a, EllPoly<K> b) => b.One.KMul(a).Sub(b);

    public static EllPoly<K> operator *(EllPoly<K> a, K b) => a.KMul(b);

    public static EllPoly<K> operator *(K a, EllPoly<K> b) => b.KMul(a);

    public static EllPoly<K> operator /(EllPoly<K> a, K b) => a.KMul(b.Inv());

    public static EllPoly<K> operator /(int a, EllPoly<K> b)
    {
        throw new NotImplementedException();
    }

    public static double Abs(EllPoly<K> t) => throw new();

    public static bool IsValuedField => false;
}