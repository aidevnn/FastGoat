using System.Collections.ObjectModel;
using FastGoat.Commons;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;

namespace FastGoat.Structures;

public static partial class Group
{
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

    public static IEnumerable<T> GenerateElements<T>(IGroup<T> bg, params T[] elements) where T : struct, IElt<T>
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

    public static ReadOnlyDictionary<T, ReadOnlyDictionary<T, int>> LongestCycles<T>(IGroup<T> g,
        IEnumerable<T> elements)
        where T : struct, IElt<T>
    {
        var set = elements.Ascending().ToHashSet();
        var allCycles = new Dictionary<T, ReadOnlyDictionary<T, int>>(set.Count);

        while (set.Count != 0)
        {
            var e0 = set.First();
            var cycle0 = Cycle(g, e0);
            set.ExceptWith(cycle0.Keys);
            if (e0.Equals(g.Neutral()))
                continue;

            allCycles[e0] = cycle0;
        }

        return new ReadOnlyDictionary<T, ReadOnlyDictionary<T, int>>(allCycles);
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
        var cosets = new Dictionary<T, Coset<T>>();
        var setH = grH.ToHashSet();
        var setG = grG.ToHashSet();
        if (!setH.IsSubsetOf(setG))
            throw new GroupException(GroupExceptionType.NotSubGroup);

        var ng = grG.Neutral();
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
        if (!group.Contains(g.Neutral()))
            return false;

        if (group.Any(e => !group.Contains(g.Invert(e))))
            return false;

        var cayleyTable = CayleyTable(g, group.ToArray());
        return IsGroup(cayleyTable);
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
        return new ConcreteGroup<T>(g, generators);
    }

    public static ConcreteGroup<T> Generate<T>(string name, IGroup<T> g, params T[] generators)
        where T : struct, IElt<T>
    {
        return new ConcreteGroup<T>(name, g, generators);
    }

    public static ConcreteGroup<T> DirectProduct<T>(ConcreteGroup<T> g1, ConcreteGroup<T> g2) where T : struct, IElt<T>
    {
        if (!g1.BaseGroup.Equals(g2.BaseGroup))
            throw new GroupException(GroupExceptionType.GroupDef);

        if (g1.SuperGroup is null || g2.SuperGroup is null || !g1.SuperGroup.Equals(g2.SuperGroup))
            throw new GroupException(GroupExceptionType.GroupDef);

        var generators = g1.PseudoGenerators.Concat(g2.PseudoGenerators).Distinct().ToArray();
        return new ConcreteGroup<T>($"{g1.Name} x {g2.Name}", g1.SuperGroup, generators);
    }

    public static ConcreteGroup<Coset<T>> Over<T>(this ConcreteGroup<T> g, ConcreteGroup<T> h)
        where T : struct, IElt<T>
    {
        // if (h.SuperGroup is null || !h.SuperGroup.Equals(g))
        //     throw new GroupException(GroupExceptionType.NotSubGroup);
        var quo = new Quotient<T>(g, h);
        return new ConcreteGroup<Coset<T>>(quo);
    }

    public static ConcreteGroup<Coset<T>> Over<T>(this ConcreteGroup<T> g, ConcreteGroup<T> h, string name)
        where T : struct, IElt<T>
    {
        if (h.SuperGroup is null || !h.SuperGroup.Equals(g))
            throw new GroupException(GroupExceptionType.NotSubGroup);
        var quo = new Quotient<T>(g, h);
        return new ConcreteGroup<Coset<T>>(name, quo);
    }

    public static ConcreteGroup<FpPolynom> Galois(char x, int q)
    {
        var cnPoly = PolynomExt.Get(q);
        var p = cnPoly.pm.p;
        var poly = new FpPolynom(p, x, cnPoly.coefs.ToArray());
        var gf = new GFp($"GF({poly})", (x, p), cnPoly.coefs);
        return new ConcreteGroup<FpPolynom>(gf);
    }

    public static ConcreteGroup<FpPolynom> Galois(int q)
    {
        return Galois('x', q);
    }
    
    public static List<WordGroup> Frobenius(int o)
    {
        var ms = IntExt.Dividors(o).Where(d => d > 1 && d % 2 == 1).ToArray();
    
        List<WordGroup> all = new();
        foreach (var m in ms)
        {
            var n = o / m;
            var rs = m.Range().Where(r => IntExt.Gcd(m, n * (r - 1)) == 1 && IntExt.PowMod(r, n, m) == 1).ToArray();
            foreach (var r in rs)
            {
                var wg = new WordGroup($"Frob({m},{n},{r})" ,$"a{m}, b{n}, b-1ab = a{r}");
                if (all.Any(g => g.IsIsomorphicTo(wg)))
                    continue;

                all.Add(wg);
            }
        }

        return all;
    }

    public static ConcreteGroup<Ep2<ZnInt, ZnInt>> MetaCyclicSdp(int m, int n, int r)
    {
        var cm = new Cn(m);
        var cn = new Cn(n);
        var autCm = AutBase(cm);
        var g1 = autCm[(cm[1], cm[r])];
        var aut = Generate(autCm, g1);
        var pMap = PartialMap((cn[0], autCm.Neutral()), (cn[1], g1));
        var theta = Hom(cn, HomomorphismMap(cn, aut, pMap));
        return SemiDirectProd($"Frob({m},{n},{r})" ,cm, theta, cn);
    }

    public static int[] MetaCyclicGetR(int m, int n)
    {
        return m.Range().Where(r => IntExt.Gcd(m, n * (r - 1)) == 1 && IntExt.PowMod(r, n, m) == 1).ToArray();
    }
    public static List<ConcreteGroup<Ep2<ZnInt,ZnInt>>> FrobeniusSdp(int order)
    {
        var ms = IntExt.Dividors(order).Where(d => d > 1 && d % 2 == 1).ToArray();
    
        List<ConcreteGroup<Ep2<ZnInt,ZnInt>>> all = new();
        foreach (var m in ms)
        {
            var n = order / m;
            var rs = MetaCyclicGetR(m, n);
            foreach (var r in rs)
            {
                var sdp = MetaCyclicSdp(m, n, r);
                if (all.Any(g => g.IsIsomorphicTo(sdp)))
                    continue;

                all.Add(sdp);
            }
        }

        return all;
    }

    public static WordGroup DiCyclic(int n) => new WordGroup($"Dic{n}", $"a{n} = b2, b2 = abab");

    public static ConcreteGroup<Mat> Quaternion(int k)
    {
        var m = Math.Log2(k);
        if (Math.Abs(m - 3.0) < 1e-5)
        {
            var gl = new GL(2, 3);
            return Generate("Q8", gl, gl[2, 2, 2, 1], gl[0, 1, 2, 0]);
        }

        if (Math.Abs(m - 4.0) < 1e-5)
        {
            var gl = new GL(2, 7);
            return Generate("Q16", gl, gl[0, 1, 6, 3], gl[6, 1, 5, 1]);
        }

        if (Math.Abs(m - 5.0) < 1e-5)
        {
            var gl = new GL(2, 17);
            return Generate("Q32", gl, gl[14, 0, 0, 11], gl[0, 16, 1, 0]);
        }

        if (Math.Abs(m - 6.0) < 1e-5)
        {
            var gl = new GL(2, 31);
            return Generate("Q64", gl, gl[16, 12, 12, 11], gl[0, 1, 30, 0]);
        }

        if (Math.Abs(m - 7.0) < 1e-5)
        {
            var gl = new GL(2, 193);
            return Generate("Q128", gl, gl[78, 38, 155, 78], gl[71, 13, 13, 122]);
        }

        throw new GroupException(GroupExceptionType.GroupDef);
    }

    public static ConcreteGroup<Ep2<ZnInt, Mat>> DiCyclicSdp(int n)
    {
        var k = IntExt.PrimesDecomposition(n).Count(i => i == 2);

        if (k != 0)
        {
            var m = n / (1 << k);
            var cm = new Cn(m);
            var q = Quaternion(1 << (k + 2));
            if (m == 1)
                return Create($"Dic{n}", Product.Group(new Cn(1), q));

            return SemiDirectProd($"Dic{n}", cm, q);
        }
        
        var gl = new GL(2, 3);
        var g1 = Generate("C4", gl, gl[0, 1, 2, 0]);
        
        var cn = new Cn(n);
        var autCn = AutBase(cn);
        var a = autCn[(cn[1], cn[n - 1])];
        var aut = Generate(autCn, a);
        var pMap = PartialMap((gl[0, 1, 2, 0], a));
        var theta = Hom(g1, HomomorphismMap(g1, aut, pMap));
        return SemiDirectProd($"Dic{n}", cn, theta, g1);
    }
}