using System.Collections;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleN;

public readonly struct CrMap<Tn, Tg> : IEnumerable<KeyValuePair<Ep<Tg>, ZNElt<Tn,Tg>>> , IEquatable<CrMap<Tn,Tg>>
    where Tg : struct, IElt<Tg> 
    where Tn : struct, IElt<Tn>
{
    private Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> Map { get; }

    public CrMap()
    {
        Map = new();
    }
    
    public CrMap(Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> map)
    {
        Map = map;
    }
    
    public int R => Map.Count == 0 ? 0 : Map.First().Key.Ei.Length;
    public AbelianDirectSum<Tn> Nab => Map.First().Value.Nab;
    public ConcreteGroup<Tn> N => Map.First().Value.N;
    public ConcreteGroup<Tg> G => Map.First().Value.G;
    public bool IsZero() => Map.Values.All(c => c.IsZero());
    public CrMap<Tn, Tg> Zero => new(Map.ToDictionary(e => e.Key, e => e.Value.ZNZero));
    public ZNElt<Tn, Tg> ZNZero => Map.First().Value.ZNZero;
    public Polynomial<ZnInt, Xi> PZero => Map.First().Value.Zero;
    public CrMap<Tn, Tg> Clone => Recreate();
    public CrMap<Tn, Tg> Recreate(Indeterminates<Xi> ind) => new(Map.ToDictionary(e => e.Key, e => e.Value.Recreate(ind)));
    public CrMap<Tn, Tg> Recreate() => new(Map.ToDictionary(e => e.Key, e => e.Value.Recreate()));

    public CrMap<Tn, Tg> Substitute(Polynomial<ZnInt, Xi> P, Xi xi) =>
        new(Map.ToDictionary(e => e.Key, e => e.Value.Substitute(P, xi).Simplify()));

    public CrMap<Tn, Tg> Substitute((Polynomial<ZnInt, Xi> P, Xi xi)[] subs) =>
        new(Map.ToDictionary(e => e.Key, e => e.Value.Substitute(subs).Simplify()));

    public ZNElt<Tn, Tg> this[Ep<Tg> index] => Map[index];

    public Dictionary<Ep<Tg>, ZNElt<Tn, Tg>>.KeyCollection Keys => Map.Keys;
    public Dictionary<Ep<Tg>, ZNElt<Tn, Tg>>.ValueCollection Values => Map.Values;
    public int Count => Map.Count;

    public CrMap<Tn, Tg>[] Split
    {
        get
        {
            var cr = this;
            var gs = Nab.DecompElementaryMap.Keys.ToArray();
            return gs.Select(t => new CrMap<Tn, Tg>(cr.ToDictionary(e => e.Key, e => e.Value.Get(t)))).ToArray();
        }
    }

    public CrMap<Tn, Tg> Add(CrMap<Tn, Tg> cr) => new(Map.ToDictionary(e => e.Key, e => (e.Value + cr[e.Key]).Simplify()));
    public CrMap<Tn, Tg> Sub(CrMap<Tn, Tg> cr) => new(Map.ToDictionary(e => e.Key, e => (e.Value - cr[e.Key]).Simplify()));
    public CrMap<Tn, Tg> Mul(Polynomial<ZnInt, Xi> c) => new(Map.ToDictionary(e => e.Key, e => (e.Value * c).Simplify()));

    public (int mod, int ord, CrMap<Tn, Tg> map)[] Generators()
    {
        var zero = PZero;
        var xis = Map.SelectMany(e => e.Value.Coefs.Values.SelectMany(p => p.ExtractAllIndeterminates)).ToHashSet();
        var allSubs = xis.OrderBy(e => e.xi).Select(xi => (xi, xis.Except(new[] { xi }).Select(yi => (zero, yi)).ToArray())).ToArray();
        var cr = this;
        var gens = allSubs.Select(subs => (subs.xi, map:cr.Substitute(subs.Item2))).ToArray();
        var nab = Nab;
        var orders = nab.ElemOrders;
        return gens.Select(g =>
                (g.xi, g.map, g.map.SelectMany(e => e.Value.Coefs
                    .Select(c => (mod: nab.DecompElementaryMap[c.Key],
                        ord: c.Value.Coefs.Values.Max(z => orders[nab.DecompElementaryMap[c.Key]][z.K]))))))
            .Select(e => (e.Item3.MaxBy(d => d.ord), map: e.map.Substitute(zero.One, e.xi)))
            .Select(e => (e.Item1.mod, e.Item1.ord, e.map))
            .OrderByDescending(e => e.map.Sum(c => (c.Value.IsZero() ? 10 : 1) + c.Value.Coefs.Count(p => p.Value.IsZero())))
            .ThenByDescending(e => e.ord)
            .ToArray();
    }

    public Xi[] ExtractAllIndeterminates => Map.Values.SelectMany(c => c.ExtractAllIndeterminates).Distinct().Order().ToArray();

    public IEnumerator<KeyValuePair<Ep<Tg>, ZNElt<Tn, Tg>>> GetEnumerator() => Map.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(CrMap<Tn, Tg> other)
    {
        foreach (var ep in Keys)
        {
            if (!this[ep].Equals(other[ep]))
                return false;
        }

        return true;
    }

    public override int GetHashCode() => Map.Count;
}