using System.Collections;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleN;

public readonly struct CrMap<Tn, Tg> : IEnumerable<KeyValuePair<Ep<Tg>, ZNElt<Tn, Tg>>>,
    IEquatable<CrMap<Tn, Tg>>, IComparable<CrMap<Tn, Tg>>
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    private Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> Map { get; }

    public CrMap()
    {
        Map = new();
    }

    public CrMap(ConcreteGroup<Tg> g, Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> map)
    {
        G = g;
        Map = map;
    }

    public CrMap(CrMap<Tn, Tg> cr)
    {
        G = cr.G;
        Map = new(cr.ToDictionary(e => e.Key, e => e.Value));
    }

    public int R => Map.Count == 0 ? 0 : Map.First().Key.Ei.Length;
    public Gp<Tg> Gr => Product.Gp(G, R);
    public AbelianDirectSum<Tn> Nab => Map.First().Value.Nab;
    public ConcreteGroup<Tn> N => Map.First().Value.N;
    public ConcreteGroup<Tg> G { get; }
    public MapElt<Tg, Automorphism<Tn>> L => Map.First().Value.L;
    public bool IsZero() => Map.Values.All(c => c.IsZero());
    public CrMap<Tn, Tg> Zero => new(G, Map.ToDictionary(e => e.Key, e => e.Value.ZNZero));
    public ZNElt<Tn, Tg> ZNZero => Map.First().Value.ZNZero;
    public Polynomial<ZnInt, Xi> ZZero => Map.First().Value.Zero;
    public CrMap<Tn, Tg> Clone() => new(this);

    public CrMap<Tn, Tg> Clone(ConcreteGroup<Tg> g)
    {
        var gp = Product.GpGenerate(g, R);
        var map0 = Map;
        var zero = ZNZero;
        var map = gp.ToDictionary(e => e, e => map0.ContainsKey(e) ? map0[e].Clone() : zero);
        return new(g, map);
    }
    public CrMap<Tn, Tg> Recreate(Indeterminates<Xi> ind) => new(G, Map.ToDictionary(e => e.Key, e => e.Value.Recreate(ind)));
    public Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> getMap() => Map.ToDictionary(e => e.Key, e => e.Value.Clone());

    public CrMap<Tn, Tg> Substitute(Polynomial<ZnInt, Xi> P, Xi xi) =>
        new(G, Map.ToDictionary(e => e.Key, e => e.Value.Substitute(P, xi).Simplify()));

    public CrMap<Tn, Tg> Substitute((Polynomial<ZnInt, Xi> P, Xi xi)[] subs) =>
        new(G, Map.ToDictionary(e => e.Key, e => e.Value.Substitute(subs).Simplify()));

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
            var g0 = G;
            return gs.Select(t => new CrMap<Tn, Tg>(g0, cr.ToDictionary(e => e.Key, e => e.Value.Get(t)))).ToArray();
        }
    }

    public CrMap<Tn, Tg> Add(CrMap<Tn, Tg> cr) => new(G, Map.ToDictionary(e => e.Key, e => (e.Value + cr[e.Key]).Simplify()));
    public CrMap<Tn, Tg> Sub(CrMap<Tn, Tg> cr) => new(G, Map.ToDictionary(e => e.Key, e => (e.Value - cr[e.Key]).Simplify()));
    public CrMap<Tn, Tg> Mul(Polynomial<ZnInt, Xi> c) => new(G, Map.ToDictionary(e => e.Key, e => (e.Value * c).Simplify()));

    public CrMap<Tn, Tg> Mul(CrMap<Tn, Tg> cr)
    {
        var map = Map.Grid2D(cr.Map).ToDictionary(e => Product.Ep(e.t1.Key.Ei.Concat(e.t2.Key.Ei).ToArray()),
            e => e.t1.Value.Mul(e.t2.Value).Simplify());
        return new(G, map);
    }
    public CrMap<Tn, Tg> Mul(int c) => new(G, Map.ToDictionary(e => e.Key, e => (c * e.Value).Simplify()));

    public CrMap<Tn, Tg> ConsTerm()
    {
        var xis = Map.SelectMany(e => e.Value.Coefs.Values.SelectMany(p => p.ExtractAllIndeterminates)).ToHashSet();
        var zero = ZZero;
        var subs = xis.Select(xi => (zero, xi)).ToArray();
        return Substitute(subs);
    }

    public CrMap<Tn, Tg>[] Cycle()
    {
        if (!Equals(ConsTerm()))
        {
            ZNSolver.DisplayCrMap("Error Cycle", this, ConsTerm());
            throw new("Only constant map");
        }

        var cycle = new List<CrMap<Tn, Tg>>() { Zero };
        while (true)
        {
            var m0 = Add(cycle.Last());
            if (m0.IsZero())
                break;
            cycle.Add(m0);
        }

        return cycle.ToArray();
    }

    public (int ord, CrMap<Tn, Tg> map)[] Generators()
    {
        var zero = ZZero;
        var xis = Map.SelectMany(e => e.Value.Coefs.Values.SelectMany(p => p.ExtractAllIndeterminates)).ToHashSet();
        if (xis.Count == 0)
            return new[] { (1, Zero) };

        var allSubs = xis.OrderBy(e => e.xi).Select(xi => (xi, xis.Except(new[] { xi }).Select(yi => (zero, yi)).ToArray())).ToArray();
        var cr = this;
        return allSubs.Select(e => cr.Substitute(e.Item2).Substitute(zero.One, e.xi)).Select(e => (e.Cycle().Length, e))
            .OrderBy(e => e.Length).ToArray();
    }

    public CrMap<Tn, Tg> Act(Tg g)
    {
        var gr = Gr;
        var map = Map.ToDictionary(e => gr.Act(g, e.Key), e => e.Value.Act(g).Simplify());
        return new(G, map);
    }

    public MapElt<Ep<Tg>, Tn> ToMapElt
    {
        get
        {
            var ct = ConsTerm();
            if (!Equals(ct))
                throw new();

            var gr = Product.GpGenerate(G, R);
            return new MapElt<Ep<Tg>, Tn>(gr, N, ct.ToDictionary(e => e.Key, e => e.Value.Expand));
        }
    }

    public CrMap<Tn, Tg>[] Decomp
    {
        get
        {
            var map = Map;
            var gr = Gr;
            var g0 = G;
            var maps = gr.ToDictionary(g => g, g => map[g].Decomp);
            var nb = maps.First().Value.Length;
            return nb.Range().Select(k => new CrMap<Tn, Tg>(g0, gr.ToDictionary(g => g, g => maps[g][k]))).ToArray();
        }
    }

    public ConcreteGroup<MapElt<Ep<Tg>, Tn>> ToGroupMapElt(string name = "")
    {
        var gens = Generators().Select(e => e.map.ToMapElt).ToArray();
        var mapGr = new MapGroupBase<Ep<Tg>, Tn>(gens.First().Domain, N);
        var g = Group.Generate(mapGr, gens);
        if (!string.IsNullOrEmpty(name))
            g.Name = name;

        return g;
    }

    public CrMap<Tn, Tg> ToHomogenous()
    {
        Ep<Tg> Chg(Ep<Tg> e, ConcreteGroup<Tg> g)
        {
            var r = e.Ei.Length + 1;
            var e0 = r.Range().Select(k => g.OpSeq(e.Ei.Take(k))).ToArray();
            return new(e0);
        }

        var G0 = G;
        var map0 = Map.ToDictionary(e => Chg(e.Key, G0), e => e.Value);
        var Gr1 = Product.Gp(G, R + 1);
        var map1 = G.SelectMany(g => map0.ToDictionary(e => Gr1.Act(g, e.Key), e => e.Value.Act(g).Simplify()))
            .ToDictionary(e => e.Key, e => e.Value);

        return new(G, map1);
    }

    public CrMap<Tn, Tg> ToInhomogenous()
    {
        Ep<Tg> Chg(Ep<Tg> e, ConcreteGroup<Tg> g)
        {
            var r = e.Ei.Length + 1;
            var e0 = r.Range().Select(k => g.OpSeq(e.Ei.Take(k))).ToArray();
            return new(e0);
        }

        var Gr1 = Product.Gp(G, R - 1);
        var G0 = G;
        var map = Map;
        var map0 = Gr1.ToDictionary(e => e, e => map[Chg(e, G0)]);
        return new(G, map0);
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

    public int CompareTo(CrMap<Tn, Tg> other)
    {
        return Map.OrderKeys(G).Select(e => e.Value).SequenceCompareTo(other.OrderKeys(G).Select(e => e.Value));
    }

    public override int GetHashCode() => Map.Count;
}