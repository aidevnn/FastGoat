using FastGoat.Commons;
using FastGoat.Structures.CartesianProduct;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Structures.GenericGroup;

public readonly struct AbelianDirectSum<T> where T : struct, IElt<T>
{
    public ConcreteGroup<T> Ab { get; }
    public ConcreteGroup<Ep<ZnInt>> AbCanonic { get; }
    public ConcreteGroup<T>[] AbSubGroupsCanonic { get; }
    private Dictionary<T, Ep<ZnInt>> ToCanonic { get; }
    private Dictionary<Ep<ZnInt>, T> FromCanonic { get; }
    public List<(T g, int o)> Decomp { get; }

    public AbelianDirectSum(ConcreteGroup<T> G)
    {
        if (G.GroupType != GroupType.AbelianGroup)
            throw new GroupException(GroupExceptionType.OnlyAbelianGroups);

        Ab = G;
        Decomp = Group.AbelianInvariants(Ab).OrderBy(e => e.o).ThenBy(e => e.g).ToList();
        var l = Decomp.Count;
        var ords = Decomp.Select(e => e.o).ToArray();
        var rg = l.Range();
        AbCanonic = FG.Abelian(Decomp.Select(e => e.o).ToArray());
        var pMap = Decomp.Select((e, i) => (e.g, e.o, i)).ToDictionary(
            e => e.g,
            e => new Ep<ZnInt>(rg.Select(k => new ZnInt(ords[k], k == e.i ? 1 : 0)).ToArray()));
        var isoMap = Group.IsomorphismMap(Ab, AbCanonic, pMap);
        ToCanonic = new(isoMap);
        FromCanonic = isoMap.ToDictionary(e => e.Value, e => e.Key);

        AbSubGroupsCanonic = Decomp.Select(e => Group.Generate($"C{e.o}", G, e.g)).ToArray();
    }

    public Ep<ZnInt> GEltToCan(T g) => ToCanonic[g];

    public Dictionary<T, ZnInt> GEltToCanMap(T g)
    {
        var decomp = Decomp.ToArray();
        return ToCanonic[g].Ei.Select((ei, k) => (ei, k)).ToDictionary(a => decomp[a.k].g, a => a.ei);
    }

    public T CanToGElt(Ep<ZnInt> z) => FromCanonic[z];
}