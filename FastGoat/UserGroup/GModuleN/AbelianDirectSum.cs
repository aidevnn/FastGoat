using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleN;

public readonly struct AbelianDirectSum<T> where T : struct, IElt<T>
{
    public ConcreteGroup<T> Ab { get; }

    public List<(T g, int o)> Decomp { get; }
    public Dictionary<T, int> DecompMap { get; }
    public ConcreteGroup<Ep<ZnInt>> AbCanonic { get; }
    private Dictionary<T, Ep<ZnInt>> ToCanonic { get; }
    private Dictionary<Ep<ZnInt>, T> FromCanonic { get; }

    public List<(T g, int o)> DecompElementary { get; }
    public Dictionary<T, int> DecompElementaryMap { get; }
    public ConcreteGroup<Ep<ZnInt>> AbElementaries { get; }
    private Dictionary<T, Ep<ZnInt>> ToElementaries { get; }
    private Dictionary<Ep<ZnInt>, T> FromElementaries { get; }

    public AbelianDirectSum(ConcreteGroup<T> G)
    {
        if (G.GroupType != GroupType.AbelianGroup)
            throw new GroupException(GroupExceptionType.OnlyAbelianGroups);

        Ab = G;
        Decomp = Group.AbelianInvariants(Ab).OrderBy(e => e.o).ThenBy(e => e.g).ToList();
        DecompMap = Decomp.ToDictionary(e => e.g, e => e.o);
        AbCanonic = FG.Abelian(Decomp.Select(e => e.o).ToArray());
        var l1 = Decomp.Count;
        var ords1 = Decomp.Select(e => e.o).ToArray();
        var rg1 = l1.Range();
        var pMap1 = Decomp.Select((e, i) => (e.g, e.o, i)).ToDictionary(
            e => e.g,
            e => new Ep<ZnInt>(rg1.Select(k => new ZnInt(ords1[k], k == e.i ? 1 : 0)).ToArray()));
        var isoMap1 = Group.IsomorphismMap(Ab, AbCanonic, pMap1);
        ToCanonic = new(isoMap1);
        FromCanonic = isoMap1.ToDictionary(e => e.Value, e => e.Key);

        DecompElementary = Decomp.SelectMany(e => (e.o == 1 ? new Dictionary<int, int>() { [1] = 1 } : IntExt.PrimesDec(e.o))
                .Select(kv => kv.Key.Pow(kv.Value))
                .Select(pk => (G.Times(e.g, e.o / pk), pk)))
            .OrderBy(e => e.pk).ToList();
        DecompElementaryMap = DecompElementary.ToDictionary(e => e.g, e => e.o);
        AbElementaries = FG.Abelian(DecompElementary.Select(e => e.o).ToArray());
        var l2 = DecompElementary.Count;
        var ords2 = DecompElementary.Select(e => e.o).ToArray();
        var rg2 = l2.Range();
        var pMap2 = DecompElementary.Select((e, i) => (e.g, e.o, i)).ToDictionary(
            e => e.g,
            e => new Ep<ZnInt>(rg2.Select(k => new ZnInt(ords2[k], k == e.i ? 1 : 0)).ToArray()));
        var isoMap2 = Group.IsomorphismMap(Ab, AbElementaries, pMap2);
        ToElementaries = new(isoMap2);
        FromElementaries = isoMap2.ToDictionary(e => e.Value, e => e.Key);

        var mods = DecompElementary.Select(e => e.o).Distinct().Order().ToDictionary(k => k, k => k.Range());
        ElemOrders = mods.Keys.ToDictionary(e => e, e => (new Cn(e)).ElementsOrders.ToDictionary(a => a.Key.K, a => a.Value));
        ElemInvertible = mods.Keys.ToDictionary(e => e, e => IntExt.UnInvertible(e));
        ElemSolve = mods.ToDictionary(
            m => m.Key,
            m => m.Value.Grid2D(m.Value).ToDictionary(
                e => e,
                e => m.Value.MinBy(k => IntExt.AmodP(e.t1 * k + e.t2, m.Key))));
    }

    public Dictionary<int, Dictionary<int, int>> ElemInvertible { get; }
    public Dictionary<int, Dictionary<int, int>> ElemOrders { get; }
    public Dictionary<int, Dictionary<(int, int), int>> ElemSolve { get; }

    public T CanToGElt(Ep<ZnInt> z) => FromCanonic[z];
    public Ep<ZnInt> GEltToCan(T g) => ToCanonic[g];

    public Dictionary<T, ZnInt> GEltToCanMap(T g)
    {
        var decomp = Decomp.ToArray();
        return ToCanonic[g].Ei.Select((ei, k) => (ei, k)).ToDictionary(a => decomp[a.k].g, a => a.ei);
    }

    public T ElemToGElt(Ep<ZnInt> z) => FromElementaries[z];
    public Ep<ZnInt> GEltToElem(T g) => ToElementaries[g];

    public Dictionary<T, ZnInt> GEltToElemMap(T g)
    {
        var decompElem = DecompElementary.ToArray();
        return ToElementaries[g].Ei.Select((ei, k) => (ei, k)).ToDictionary(a => decompElem[a.k].g, a => a.ei);
    }
}