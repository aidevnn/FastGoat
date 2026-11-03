using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Subgroups;

namespace Craft.Craft;

public static partial class GroupCraft
{
    public static Comparer<T> ByOrder<T>(this ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        return Comparer<T>.Create((a, b) => g.ElementsOrders[a].CompareTo(g.ElementsOrders[b]));
    }

    public static Comparer<T> ByOrderDescending<T>(this ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        return Comparer<T>.Create((a, b) => -g.ElementsOrders[a].CompareTo(g.ElementsOrders[b]));
    }

    public static bool AreConjugate<T>(ConcreteGroup<T> g, GroupSubset<T> sg1, GroupSubset<T> sg2)
        where T : struct, IElt<T>
    {
        return sg1.Equals(sg2) || g.Any(x => sg1.SetEquals(sg2.Select(s0 => g.OpSeq(x, s0, g.Invert(x)))));
    }

    public static List<GroupSubset<T>> SubGroupsConjugates<T>(ConcreteGroup<T> g, GroupSubset<T> h)
        where T : struct, IElt<T>
    {
        if (!h.SubSetOf(g))
            throw new GroupException(GroupExceptionType.NotSubGroup);

        if (g.GroupType == GroupType.AbelianGroup)
            return [h];

        var all = new List<GroupSubset<T>>();
        var gens = g.GetGenerators().ToHashSet();
        var conjsH = Group.Orbits(gens, Group.ByConjugateSet(g), h);
        foreach (var cj in conjsH)
        {
            var gensConjs = cj.Generators.ToArray();
            var elts = Group.GenerateElements(g, gensConjs);
            all.Add(new([..gensConjs], elts));
        }

        return all;
    }

    public static bool IsInnerAutomorphism<T>(this Automorphism<T> aut) where T : struct, IElt<T>
    {
        var G = aut.Domain;
        var act = Group.ByConjugate(G);
        return G.Any(g => aut.AutMap.All(f => aut[f.Key].Equals(act(g, f.Key))));
    }

    public static T[] RecreateGenerators<T>(IGroup<T> g, T[] generators, Comparer<T> comp)
        where T : struct, IElt<T>
    {
        if (generators.Length == 1)
            return [g.Neutral()];

        HashSet<T> tmpElements = new() { g.Neutral() };
        List<T> newGens = new();
        var set = generators.ToHashSet();
        set.ExceptWith(tmpElements);
        var setG = Group.GenerateElements(g, generators);
        do
        {
            newGens.Add(set.OrderBy(e => e, comp).First());
            tmpElements = Group.GenerateElements(g, tmpElements, newGens);
            set.ExceptWith(tmpElements);
        } while (setG.Count != tmpElements.Count && set.Any());

        return newGens.ToArray();
    }

    public static T[] RecreateGenerators<T>(IGroup<T> g, T[] generators) where T : struct, IElt<T>
    {
        var idx = generators.Index().ToDictionary(e => e.Item, e => e.Index);
        var comp = Comparer<T>.Create((a, b) => idx[a].CompareTo(idx[b]));
        return RecreateGenerators(g, generators, comp);
    }

    public static T[] RecreateGenerators<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var classes = AllOrbits(g);
        var gens = new List<T>();
        foreach (var (_, conjs) in classes.GroupBy(e => g.ElementsOrders[e.Key])
                     .ToDictionary(e => e.Key, e => e.OrderBy(f => f.Value.Count).ToArray())
                     .OrderByDescending(e => e.Key))
        {
            gens.AddRange(conjs.Select(e => e.Key));
            gens.AddRange(conjs.SelectMany(f => f.Value).Where(e => !gens.Contains(e)));
        }

        return RecreateGenerators(g, gens.ToArray(), g.ByOrderDescending());
    }

    public static T[] RecreateGenerators<T>(ConcreteGroup<T> g, Comparer<T> comp) where T : struct, IElt<T>
    {
        return RecreateGenerators(g, g.ToArray(), comp);
    }

    public static HashSet<T> GenerateElementsLimited<T>(IGroup<T> bg, T[] generators, int limits)
        where T : struct, IElt<T>
    {
        var elements = new HashSet<T>([bg.Neutral()]);
        var q = new Queue<T>(elements);
        HashSet<T> generatedElements = new HashSet<T>(elements);
        while (q.Count != 0 && limits >= generatedElements.Count)
        {
            var e1 = q.Dequeue();
            foreach (var e2 in generators)
            {
                var e3 = bg.Op(e2, e1);
                if (generatedElements.Add(e3))
                    q.Enqueue(e3);

                if (limits < generatedElements.Count)
                    return [];
            }
        }

        if (limits < generatedElements.Count)
            return [];

        return generatedElements;
    }

    public static HashSet<T> GenerateElementsLimited<T>(IGroup<T> bg, HashSet<T> elements, IEnumerable<T> generators,
        int limits)
        where T : struct, IElt<T>
    {
        if (!elements.Contains(bg.Neutral()))
            throw new GroupException(GroupExceptionType.BaseGroup);

        var q = new Queue<T>(elements);
        HashSet<T> generatedElements = new HashSet<T>(elements);
        var gens = generators.ToArray();
        while (q.Count != 0)
        {
            var e1 = q.Dequeue();
            foreach (var e2 in gens)
            {
                var e3 = bg.Op(e2, e1);
                if (generatedElements.Add(e3))
                    q.Enqueue(e3);

                if (limits < generatedElements.Count)
                    return [];
            }
        }

        if (limits < generatedElements.Count)
            return [];

        return generatedElements;
    }

    public static HashSet<T> Orbits<T>(HashSet<T> gens, GroupAction<T, T> act, T x) where T : struct, IElt<T>
    {
        HashSet<T> set = new() { x };
        Queue<T> q = new Queue<T>();
        q.Enqueue(x);
        while (q.Count != 0)
        {
            var x0 = q.Dequeue();
            foreach (var g in gens)
            {
                var x1 = act(g, x0);
                if (set.Add(x1))
                    q.Enqueue(x1);
            }
        }

        return set;
    }

    public static Dictionary<T, HashSet<T>> AllOrbits<T>(IGroup<T> G, HashSet<T> gens, T[] set, GroupAction<T, T> act)
        where T : struct, IElt<T>
    {
        var setG = set.ToHashSet();
        var allOrbx = new Dictionary<T, HashSet<T>>();
        while (setG.Count != 0)
        {
            var x = setG.First();
            var orbx = Orbits(gens, act, x);
            setG.ExceptWith(orbx);
            allOrbx[x] = orbx;
        }

        return allOrbx;
    }

    public static Dictionary<T, HashSet<T>> AllOrbits<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        var bs = G.BaseGroup;
        return AllOrbits(bs, G.GetGenerators().ToHashSet(), G.ToArray(), Group.ByConjugate(bs));
    }

    public static bool IsPGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        return IntExt.PrimesDec(G.Order).Count == 1;
    }

    public static bool IsSylow<T>(ConcreteGroup<T> G, ConcreteGroup<T> H) where T : struct, IElt<T>
    {
        return IsPGroup(H) && IntExt.Gcd(G.Order / H.Order, H.Order) == 1;
    }

    public static bool HasComplement<T>(ConcreteGroup<T> G, ConcreteGroup<T> H) where T : struct, IElt<T>
    {
        if (!IsSylow(G, H))
            return false;

        var (C, N) = (Group.Centralize(G, H), Group.Normalize(G, H));
        return C.SetEquals(N);
    }

    public static void DisplayIsNormal<T>(ConcreteGroup<T> G, ConcreteGroup<T> H) where T : struct, IElt<T>
    {
        Console.WriteLine($"{H} IsNormalSubgroup of {G}:{Group.IsNormalSubgroup(G, H)}");
        Console.WriteLine();
    }

    public static GroupSubset<T> MonogenicSet<T>(ConcreteGroup<T> G, T e) where T : struct, IElt<T>
    {
        return new([e], Group.GenerateElements(G, e));
    }

    public static ConcreteGroup<T> MonogenicGroup<T>(ConcreteGroup<T> G, T e) where T : struct, IElt<T>
    {
        var ord = G.ElementsOrders[e];
        return Group.Generate($"C{ord}", G, e);
    }
    
    public static ConcreteGroup<T> MonogenicGroup<T>(string name, ConcreteGroup<T> G, T e) where T : struct, IElt<T>
    {
        return Group.Generate(name, G, e);
    }
    
    public static ConcreteGroup<T> DirectProduct<T>(ConcreteGroup<T> g1, ConcreteGroup<T> g2) where T : struct, IElt<T>
    {
        var (ng1, ng2) = (g1.NameParenthesis(), g2.NameParenthesis());
        var prod = Group.DirectProduct($"{ng1} * {ng2}", g1, g2);
        var norm1 = Group.IsNormalSubgroup(prod, g1);
        var norm2 = Group.IsNormalSubgroup(prod, g2);
        if (prod.Order != g1.Order * g2.Order)
            return prod;

        if (norm1 && norm2)
        {
            prod.Name = $"{ng1} x {ng2}";
            return prod;
        }

        if (norm1)
        {
            prod.Name = $"{ng1} x: {ng2}";
            return prod;
        }

        prod.Name = $"{ng2} x: {ng1}";
        return prod;
    }

}