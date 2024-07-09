using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.Characters;

public readonly struct Character<T> : IElt<Character<T>>, IRingElt<Character<T>>, IFieldElt<Character<T>> where T : struct, IElt<T>
{
    public ConjugacyClasses<T> Classes { get; }
    public ConcreteGroup<T> Gr => Classes.Gr;

    public Character(ConcreteGroup<T> gr)
    {
        Classes = Group.ConjugacyClasses(gr);
        var map = Map = Classes.ToDictionary(e => e, _ => (Cnf?)Cnf.CnfOne);
        Hash = Classes.Aggregate(0, (hash, g) => (hash, map[g]!.Value.GetHashCodeSlow()).GetHashCode());
    }

    public Character(ConjugacyClasses<T> classes, Dictionary<T, Cnf?> map)
    {
        Map = map;
        Classes = classes;
        Hash = Classes.Aggregate(0,
            (hash, g) => (hash, map[g].HasValue ? map[g]!.Value.GetHashCodeSlow() : 0).GetHashCode());
    }

    public Dictionary<T, Cnf?> Map { get; }

    public int Hash { get; }

    public IEnumerable<T> Kernel()
    {
        var chiNeutral = this[Gr.Neutral()]!.Value;
        foreach (var g in Gr)
        {
            if (this[g]!.Value.Equals(chiNeutral))
                yield return g;
        }
    }

    public IEnumerable<T> Centre()
    {
        var c0 = Dim;
        foreach (var g in Gr)
        {
            var abs = this[g]!.Value.Module;
            if (Math.Abs(abs - c0) < 1e-5)
                yield return g;
        }
    }

    public int Dim => (int)Map[Gr.Neutral()]?.Simplify().E[0].Num!;

    public bool HasAllValues => Map.Values.All(c => c.HasValue);
    public bool IsLinear => HasAllValues && (Map[Gr.Neutral()] - 1)!.Value.IsZero();
    public bool IsIdentity => Map.Values.All(c => c.HasValue && c.Value.Equals(Cnf.CnfOne));
    public bool IsFaithfull => HasAllValues && Kernel().Count() == 1;
    public Character<T> Simplify() => new(Classes, Map.ToDictionary(e => e.Key, e => (Cnf?)e.Value!.Value.Simplify()));

    public void SetValue(Cnf c, T e)
    {
        Map[e] = c;
    }

    public Cnf? this[T g] => Map.ContainsKey(g) ? Map[g] : Map[Classes.GetRepresentative(g)];

    public bool Equals(Character<T> other)
    {
        if (!Gr.SetEquals(other.Gr))
            return false;

        foreach (var g in Classes)
        {
            var c0 = Map.TryGetValue(g, out Cnf? e0) && e0.HasValue;
            var c1 = other.Map.TryGetValue(g, out Cnf? e1) && e1.HasValue;
            if (!c0 || !c1 || !e0!.Value.Equals(e1!.Value))
                return false;
        }

        return true;
    }

    public int CompareTo(Character<T> other)
    {
        if (!Gr.SetEquals(other.Gr))
            return 1;

        // Character One priority 
        if (Equals(One))
            return -1;

        if (other.Equals(One))
            return 1;

        // Then by dimension
        var dimThis = this[Gr.Neutral()]?.Module ?? Double.PositiveInfinity;
        var dimOther = other[Gr.Neutral()]?.Module ?? Double.PositiveInfinity;
        var compDim = dimThis.CompareTo(dimOther);
        if (compDim != 0)
            return compDim;

        // Then all values
        var avThis = HasAllValues;
        var avOther = other.HasAllValues;
        var compAV = avThis.CompareTo(avOther);
        if (compAV != 0 || !avThis)
            return compAV;

        // Then rational count
        var t = this;
        var ratThis = Classes.Count(e => (t[e]?.E.Degree ?? 1) == 0);
        var ratOther = Classes.Count(e => (other[e]?.E.Degree ?? 1) == 0);
        var compRat = ratThis.CompareTo(ratOther);
        if (compRat != 0)
            return -compRat;

        // Then by content
        var cls = Classes;
        foreach (var e in Classes.OrderBy(e => cls.GetIndex(e)))
        {
            var (cnfThis, cnfOther) = (this[e], other[e]);
            var cplxThis = cnfThis?.ToComplex ?? Complex.Infinity;
            var cplxOther = cnfOther?.ToComplex ?? Complex.Infinity;
            var compMod = cplxThis.Magnitude.CompareTo(cplxOther.Magnitude);
            if (compMod != 0)
                return compMod;

            var compPhase = cplxThis.Phase.CompareTo(cplxOther.Phase);
            if (compPhase != 0)
                return compPhase;
        }

        return 0;
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var str = new List<string>();
        var cl = Classes;
        foreach (var (e, c0) in Map.OrderBy(e => cl.GetIndex(e.Key)))
        {
            var t = Classes.GetClassName(e);
            if (!c0.HasValue)
            {
                str.Add($"#({t})");
                continue;
            }

            var c = c0.Value;
            if (c.IsZero())
                continue;

            var cstr = $"{FG.PrettyPrintCnf(c).c}";
            if (Gr.Neutral().Equals(e))
            {
                if (cstr.Contains(" + ") || cstr.Contains(" - "))
                    str.Add($"({cstr})");
                else
                    str.Add($"{cstr}");
            }
            else
            {
                if (cstr.Contains(" + ") || cstr.Contains(" - "))
                    str.Add($"({cstr})({t})");
                else
                {
                    if (c.Equals(c.One))
                        str.Add($"({t})");
                    else if (c.Equals(-c.One))
                        str.Add($"-({t})");
                    else
                        str.Add($"{cstr}({t})");
                }
            }
        }

        return str.Count == 0 ? "0" : str.Glue(" + ").Replace("+ -", "- ");
    }

    public bool IsZero() => Map.Values.All(c => c.HasValue && c.Value.IsZero());

    public Character<T> Zero => new(Classes, Map.ToDictionary(e => e.Key, _ => (Cnf?)null));
    public Character<T> One => new(Classes, Map.ToDictionary(e => e.Key, _ => (Cnf?)Cnf.CnfOne));

    public Character<T> Add(Character<T> e) =>
        new(Classes, Map.ToDictionary(k => k.Key, k => (Cnf?)(k.Value!.Value + e[k.Key]!.Value)));

    public Character<T> Sub(Character<T> e) =>
        new(Classes, Map.ToDictionary(k => k.Key, k => (Cnf?)(k.Value!.Value - e[k.Key]!.Value)));

    public Character<T> Opp() => new(Classes, Map.ToDictionary(k => k.Key, k => (Cnf?)-k.Value!.Value));

    public Character<T> Mul(Character<T> e) =>
        new(Classes, Map.ToDictionary(k => k.Key, k => (Cnf?)(k.Value!.Value * e[k.Key]!.Value)));

    public (Character<T> quo, Character<T> rem) Div(Character<T> e)
    {
        throw new NotImplementedException();
    }

    public Character<T> Mul(int k) => new(Classes, Map.ToDictionary(e => e.Key, e => (Cnf?)(e.Value!.Value * k)));

    public Character<T> Pow(int k) => new(Classes, Map.ToDictionary(e => e.Key, e => (Cnf?)(e.Value!.Value.Pow(k))));

    public static Character<T> operator +(Character<T> a, Character<T> b) => a.Add(b);

    public static Character<T> operator +(int a, Character<T> b) => b.One.Mul(a).Add(b);

    public static Character<T> operator +(Character<T> a, int b) => a.Add(a.One.Mul(b));

    public static Character<T> operator -(Character<T> a) => a.Opp();

    public static Character<T> operator -(Character<T> a, Character<T> b) => a.Sub(b);

    public static Character<T> operator -(int a, Character<T> b) => b.One.Mul(a).Sub(b);

    public static Character<T> operator -(Character<T> a, int b) => a.Sub(a.One.Mul(b));

    public static Character<T> operator *(Character<T> a, Character<T> b) => a.Mul(b);

    public static Character<T> operator *(int a, Character<T> b) => b.Mul(a);

    public static Character<T> operator *(Cnf a, Character<T> b) =>
        new(b.Classes, b.Map.ToDictionary(e => e.Key, e => e.Value * a));

    public static Character<T> operator *(Character<T> a, int b) => a.Mul(b);

    public static Character<T> operator /(Character<T> a, Character<T> b) => throw new NotImplementedException();

    public static Character<T> operator /(Character<T> a, int b) =>
        new(a.Classes, a.Map.ToDictionary(e => e.Key, e => (Cnf?)(e.Value!.Value / b)));

    public int P => 0;
    public Character<T> Inv()
    {
        if (!Invertible())
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(Classes, Map.ToDictionary(e => e.Key, e => (Cnf?)e.Value!.Value.Inv()));
    }

    public bool Invertible() => Map.Values.All(c => c.HasValue && !c.Value.IsZero());

    public static Character<T> operator /(int a, Character<T> b)
    {
        return a * b.Inv();
    }

    public static double Abs(Character<T> t) => double.Sqrt(Cnf.Abs(FG.InnerProduct(t, t)));

    public static bool IsValuedField => true;
}