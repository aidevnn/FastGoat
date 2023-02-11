using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.Characters;

public readonly struct Character<T> : IElt<Character<T>>, IRingElt<Character<T>> where T : struct, IElt<T>
{
    public ConjugacyClasses<T> Classes { get; }
    public ConcreteGroup<T> Gr => Classes.Gr;

    public Character(ConcreteGroup<T> gr)
    {
        Classes = Group.ConjugacyClasses(gr);
        Map = Classes.ToDictionary(e => e, _ => (Cnf?)Cnf.CnfOne);
        Hash = Map.Aggregate(0, (hash, kp) => (hash, kp.Key.Hash, kp.Value.GetHashCode()).GetHashCode());
    }

    public Character(ConjugacyClasses<T> classes, Dictionary<T, Cnf?> map)
    {
        Map = map;
        Classes = classes;
        Hash = Map.Aggregate(0, (hash, kp) => (hash, kp.Key.Hash, kp.Value?.GetHashCode() ?? 0).GetHashCode());
    }

    public Dictionary<T, Cnf?> Map { get; }

    public int Hash { get; }

    public IEnumerable<T> Kernel()
    {
        var chiNeutral = this[Gr.Neutral()];
        foreach (var g in Gr)
        {
            if (this[g].Equals(chiNeutral))
                yield return g;
        }
    }

    public IEnumerable<T> Centre()
    {
        var c0 = this[Gr.Neutral()]!.Value.Module;
        foreach (var g in Gr)
        {
            var abs = this[g]!.Value.Module;
            if (Math.Abs(abs - c0) < 1e-5)
                yield return g;
        }
    }

    public bool HasAllValues => Map.Values.All(c => c.HasValue);
    public bool IsIdentity => Map.Values.All(c => c.HasValue && c.Value.Equals(Cnf.CnfOne));
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
            var c0 = Map.ContainsKey(g);
            var c1 = other.Map.ContainsKey(g);
            if (!c0 || !c1 || !Map[g]!.Value.Equals(other.Map[g]!.Value))
                return false;
        }

        return true;
    }

    public int CompareTo(Character<T> other)
    {
        if (!Gr.SetEquals(other.Gr))
            return 1;

        foreach (var e in Classes)
        {
            var x = this[e]?.Module ?? Double.PositiveInfinity;
            var y = other[e]?.Module ?? Double.PositiveInfinity;
            var compMod = x.CompareTo(y);
            if (compMod != 0)
                return compMod;
        }

        return 0;
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var str = new List<string>();
        var cl = Classes;
        foreach (var (e, c0) in Map.OrderBy(e => cl.GetClassName(e.Key)))
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

            var cstr = $"{c}";
            if (Gr.Neutral().Equals(e))
            {
                if (cstr.Contains('+'))
                    str.Add($"({cstr})");
                else
                    str.Add($"{cstr}");
            }
            else
            {
                if (cstr.Contains('+'))
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

        return str.Count == 0 ? "0" : str.Glue(" + ");
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
        var map = Map.ToDictionary(k => k.Key, k => (Cnf?)(k.Value!.Value / e[k.Key]!.Value));
        return (new(Classes, map), Zero);
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
    public static Character<T> operator *(Cnf a, Character<T> b) => new(b.Classes, b.Map.ToDictionary(e => e.Key, e => e.Value * a));

    public static Character<T> operator *(Character<T> a, int b) => a.Mul(b);

    public static Character<T> operator /(Character<T> a, Character<T> b) => a.Div(b).quo;

    public static Character<T> operator /(Character<T> a, int b) =>
        new(a.Classes, a.Map.ToDictionary(e => e.Key, e => (Cnf?)(e.Value!.Value / b)));
}