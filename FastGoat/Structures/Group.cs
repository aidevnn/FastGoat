using System.Collections.ObjectModel;
using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures;

public static partial class Group
{
    public static string WithParenthesis(this string name) => (name.First() == '(' && name.Last() == ')') || !name.Contains(' ')
        ? name
        : $"({name})";
    public static string NameParenthesis<T>(this IGroup<T> g) where T : struct, IElt<T> => g.Name.WithParenthesis();

    public static T Times<T>(this IGroup<T> g, T e, int p) where T : struct, IElt<T>
    {
        var acc = g.Neutral();
        if (p == 0)
            return acc;

        var e0 = p > 0 ? e : g.Invert(e);
        var p0 = p > 0 ? p : -p;
        for (var k = 0; k < p0; ++k)
            acc = g.Op(e0, acc);

        return acc;
    }

    public static T OpSeq<T>(this IGroup<T> g, IEnumerable<T> seq) where T : struct, IElt<T>
    {
        return seq.Aggregate(g.Neutral(), (acc, e) => g.Op(acc, e));
    }

    public static HashSet<T> GenerateElements<T>(IGroup<T> bg, HashSet<T> elements, List<T> generators)
        where T : struct, IElt<T>
    {
        if (!elements.Contains(bg.Neutral()))
            throw new GroupException(GroupExceptionType.BaseGroup);

        var q = new Queue<T>(elements);
        HashSet<T> generatedElements = new HashSet<T>(elements);
        while (q.Count != 0)
        {
            var e1 = q.Dequeue();
            foreach (var e2 in generators)
            {
                var e3 = bg.Op(e2, e1);
                if (generatedElements.Add(e3))
                    q.Enqueue(e3);
            }
        }

        return generatedElements;
    }

    private static HashSet<T> GenerateElements<T>(IGroup<T> bg, List<T> generators) where T : struct, IElt<T>
    {
        return GenerateElements(bg, new HashSet<T>() { bg.Neutral() }, generators);
    }

    public static HashSet<T> GenerateElements<T>(IGroup<T> bg, params T[] elements) where T : struct, IElt<T>
    {
        return GenerateElements(bg, elements.ToList());
    }

    public static ReadOnlyDictionary<T, int> Cycle<T>(IGroup<T> g, T e) where T : struct, IElt<T>
    {
        var e0 = e;
        var n = g.Neutral();
        var p = 1;
        var elements = new Dictionary<T, int> { [e] = p };
        while (!n.Equals(e0))
        {
            e0 = g.Op(e, e0);
            elements[e0] = ++p;
        }

        return new ReadOnlyDictionary<T, int>(elements);
    }
    
    public static IEnumerable<T> CycleExceptNeutral<T>(ConcreteGroup<T> g, T e) where T : struct, IElt<T>
    {
        var a = g.Neutral();
        var o = g.ElementsOrders[e];
        for (int i = 0; i < o - 1; i++)
        {
            a = g.Op(a, e);
            yield return a;
        }
    }

    public static ReadOnlyDictionary<T, int> ElementsOrders<T>(IGroup<T> g, IEnumerable<T> elements)
        where T : struct, IElt<T>
    {
        var orders = new Dictionary<T, int>();
        foreach (var e in elements)
        {
            var cycle = Cycle(g, e);
            orders[e] = cycle.Count;
        }

        return new ReadOnlyDictionary<T, int>(orders);
    }

    public static bool IsCommutative<T>(IGroup<T> g, IEnumerable<T> ts) where T : struct, IElt<T>
    {
        var ts0 = ts.ToArray();
        foreach (var e1 in ts0)
        foreach (var e2 in ts0)
        {
            var e12 = g.Op(e1, e2);
            var e21 = g.Op(e2, e1);
            if (!e12.Equals(e21))
                return false;
        }

        return true;
    }

    public static (HashSet<T> elements, List<T> uniqueGenerators) UniqueGenerators<T>(IGroup<T> g, T[] generators)
        where T : struct, IElt<T>
    {
        HashSet<T> tmpElements = new() { g.Neutral() };
        List<T> uniqueGenerators = new();
        foreach (var elt in generators)
        {
            if (tmpElements.Contains(elt))
                continue;

            uniqueGenerators.Add(elt);
            tmpElements = GenerateElements(g, tmpElements, uniqueGenerators);
        }

        return (tmpElements, uniqueGenerators);
    }

    public static Dictionary<T, Coset<T>> Cosets<T>(ConcreteGroup<T> grG, ConcreteGroup<T> grH, Comparer<T> comparer,
        CosetType cosetType = CosetType.Both) where T : struct, IElt<T>
    {
        var setH = grH.ToHashSet();
        var setG = grG.ToHashSet();
        var cosets = new Dictionary<T, Coset<T>>(setG.Count / setH.Count);
        if (!setH.IsSubsetOf(setG) || !grH.BaseGroup.Equals(grG.BaseGroup))
            throw new GroupException(GroupExceptionType.NotSubGroup);

        var ngH = new Coset<T>(grG, grH);
        foreach (var h in setH)
            cosets[h] = ngH;

        foreach (var x in grG.OrderBy(a => a, comparer))
        {
            if (cosets.ContainsKey(x))
                continue;

            var xi = grG.Invert(x);
            var xH = new Coset<T>(grG, grH, x);
            foreach (var xh in xH)
            {
                cosets[xh] = xH;
                if (cosetType == CosetType.Both)
                {
                    if (!setH.Contains(grG.Op(xh, xi)))
                        throw new GroupException(GroupExceptionType.NotNormal);
                }
            }
        }

        return cosets;
    }

    public static Dictionary<T, Coset<T>> Cosets<T>(ConcreteGroup<T> grG, ConcreteGroup<T> grH,
        CosetType cosetType = CosetType.Both) where T : struct, IElt<T>
    {
        return Cosets(grG, grH, Comparer<T>.Default, cosetType);
    }

    public static T[,] CayleyTable<T>(IGroup<T> g, T[] elements) where T : struct, IElt<T>
    {
        var cayleyTable = new T[elements.Length, elements.Length];
        var lt = Enumerable.Range(0, elements.Length).ToArray();
        foreach (var i in lt)
        {
            foreach (var j in lt)
            {
                cayleyTable[i, j] = g.Op(elements[i], elements[j]);
            }
        }

        return cayleyTable;
    }

    public static bool IsGroup<T>(T[,] table) where T : struct, IElt<T>
    {
        var n = table.GetLength(0);
        if (n != table.GetLength(1))
            return false;

        var lt = Enumerable.Range(0, n).ToArray();
        var group = lt.Select(i => table[0, i]).ToHashSet();
        if (group.Count != n)
            return false;

        var checkRows = lt.Select(i0 => lt.Select(j0 => table[i0, j0]).ToHashSet()).All(group.SetEquals);
        var checkCols = lt.Select(i0 => lt.Select(j0 => table[j0, i0]).ToHashSet()).All(group.SetEquals);
        return checkRows && checkCols;
    }

    public static bool IsGroup<T>(IGroup<T> g, IEnumerable<T> elements) where T : struct, IElt<T>
    {
        var group = elements.ToHashSet();
        var n = g.Neutral();
        if (!group.Contains(n))
            return false;

        if (group.Select(e => (e, ei: g.Invert(e))).Any(a => !group.Contains(a.ei) || !g.Op(a.e, a.ei).Equals(n)))
            return false;

        var cayleyTable = CayleyTable(g, group.ToArray());
        var isLatinSquare = IsGroup(cayleyTable);
        if (!isLatinSquare)
            return false;

        return group.Grid3D(group, group).All(e => g.Op(e.t1, g.Op(e.t2, e.t3)).Equals(g.Op(g.Op(e.t1, e.t2), e.t3)));
    }

    public static bool IsGroup<T>(IGroup<T> g) where T : struct, IElt<T>
    {
        return IsGroup(g, g.GetElements().Ascending());
    }

    public static ConcreteGroup<T> Create<T>(string name, IGroup<T> g) where T : struct, IElt<T>
    {
        return new ConcreteGroup<T>(name, g);
    }

    public static ConcreteGroup<T> Create<T>(IGroup<T> g) where T : struct, IElt<T>
    {
        return new ConcreteGroup<T>(g);
    }

    public static ConcreteGroup<T> Generate<T>(IGroup<T> g, params T[] generators) where T : struct, IElt<T>
    {
        generators = generators.Length == 0 ? g.GetGenerators().ToArray() : generators;
        return new ConcreteGroup<T>(g, generators);
    }

    public static ConcreteGroup<T> Generate<T>(string name, IGroup<T> g, params T[] generators)
        where T : struct, IElt<T>
    {
        generators = generators.Length == 0 ? g.GetGenerators().ToArray() : generators;
        return new ConcreteGroup<T>(name, g, generators);
    }

    public static ConcreteGroup<T> DirectProduct<T>(string name, ConcreteGroup<T> g1, ConcreteGroup<T> g2) where T : struct, IElt<T>
    {
        if (!g1.BaseGroup.Equals(g2.BaseGroup))
            throw new GroupException(GroupExceptionType.GroupDef);

        var generators = g1.GetGenerators().Union(g2.GetGenerators()).ToArray();
        if (g1.SuperGroup is not null && g1.SuperGroup.SuperSetOf(generators))
            return new ConcreteGroup<T>(name, g1.SuperGroup, generators);

        return new ConcreteGroup<T>(name, g1.BaseGroup, generators);
    }

    public static ConcreteGroup<T> DirectProduct<T>(ConcreteGroup<T> g1, ConcreteGroup<T> g2) where T : struct, IElt<T>
    {
        return DirectProduct($"{g1.NameParenthesis()} x {g2.NameParenthesis()}", g1, g2);
    }

    public static ConcreteGroup<Coset<T>> Over<T>(this ConcreteGroup<T> g, ConcreteGroup<T> h)
        where T : struct, IElt<T>
    {
        var quo = new Quotient<T>(g, h);
        return new ConcreteGroup<Coset<T>>(quo);
    }

    public static ConcreteGroup<Coset<T>> Over<T>(this ConcreteGroup<T> g, ConcreteGroup<T> h, string name)
        where T : struct, IElt<T>
    {
        var quo = new Quotient<T>(g, h);
        return new ConcreteGroup<Coset<T>>(name, quo);
    }

    public static ConcreteGroup<K> AddGroup<K>(string name, params K[] gens)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var fg = new OpGroup<K>(name, gens, FGroupOp.Additive);
        return new ConcreteGroup<K>(fg);
    }

    public static ConcreteGroup<K> MulGroup<K>(string name, params K[] gens)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var fg = new OpGroup<K>(name, gens, FGroupOp.Multiplicative);
        return new ConcreteGroup<K>(fg);
    }
}