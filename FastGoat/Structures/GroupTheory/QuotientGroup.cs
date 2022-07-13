using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;

public class QuotientGroup<U> : SubGroup<U> where U : struct, IElt<U>
{
    public QuotientGroup(SubGroup<U> g, SubGroup<U> h) : base(g.UpperGroup)
    {
        G = new GroupElement<U>(g.UpperGroup, g.AllElements()) { SortBy = SortBy.Value };
        H = h;
        representatives = new Dictionary<U, U>(g.Count, new EltEquality<U>());
        classOf = new Dictionary<U, List<U>>(g.Count / h.Count, new EltEquality<U>());

        Init();
        SetName($"{g.Infos.Name}/{h.Infos.Name}");
    }

    readonly Dictionary<U, U> representatives;
    readonly Dictionary<U, List<U>> classOf;
    SubGroup<U> G { get; }
    SubGroup<U> H { get; }

    void Init()
    {
        if (!G.UpperGroup.Equals(H.UpperGroup))
            return;

        if (!G.IsGroup() || !H.IsGroup())
            return;

        var h = H.AllElements().ToHashSet(new EltEquality<U>());
        if (!h.IsSubsetOf(G.AllElements()))
            return;

        HashSet<SubSet<U>> GH = new HashSet<SubSet<U>>(new EqSubSet<U>());
        HashSet<SubSet<U>> HG = new HashSet<SubSet<U>>(new EqSubSet<U>());
        var listG = G.AllElements().ToList();
        listG.Sort(G.CompareElt);

        foreach (var x in listG)
        {
            GH.Add(x.Gop(H));
            HG.Add(H.Gop(x));
        }

        if (!GH.SetEquals(HG))
            return;

        foreach (var lh in GH)
        {
            var lu = new List<U>(lh.AllElements());
            lu.Sort((a, b) => a.CompareTo(b));
            var r = lu.First();
            classOf[r] = lu;
            AddElement(r);
            lu.ForEach(e => representatives[e] = r);
        }
    }

    public override U Neutral => UpperGroup.Neutral;
    public override U Invert(U a)
    {
        var ia = G.Invert(a);
        return representatives[ia];
    }

    public override U Op(U a, U b)
    {
        var c = G.Op(a, b);
        return representatives[c];
    }

    public U ClassOf(U e) => representatives[e];

    public void DisplayClasses()
    {
        foreach (var kp in classOf)
        {
            Console.WriteLine("Class of : {0}", kp.Key);
            foreach (var e in kp.Value)
                Console.WriteLine("\t{0}", e);
        }

        Console.WriteLine();
    }
}

public static partial class GroupExt
{
    public static QuotientGroup<U> Over<U>(this SubGroup<U> g, SubGroup<U> h) where U : struct, IElt<U>
    {
        return new QuotientGroup<U>(g, h);
    }
}