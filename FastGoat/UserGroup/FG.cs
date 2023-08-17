using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;

namespace FastGoat.UserGroup;

public static partial class FG
{
    private static Dictionary<int, KPoly<Rational>> CyclotomicPolynomials { get; }

    static FG()
    {
        var x = QPoly('X');
        CyclotomicPolynomials = new() { [1] = x - 1 };
    }

    public static KPoly<Rational> CyclotomicPolynomial(int k)
    {
        if (k < 1)
            throw new Exception();
        
        if (CyclotomicPolynomials.TryGetValue(k, out var polynomial))
            return polynomial;

        var x = QPoly('X');
        var divs = IntExt.Dividors(k);
        var prod = x.One;
        foreach (var k0 in divs)
            prod *= CyclotomicPolynomial(k0); // Recursion

        var poly = x.Pow(k) - 1;
        var phi = CyclotomicPolynomials[k] = poly / prod; // Memoization
        return phi;
    }

    public static ConcreteGroup<Perm> Symmetric(int n) => new Symm(n);

    public static ConcreteGroup<Perm> Alternate(int n)
    {
        if (n < 3)
            throw new GroupException(GroupExceptionType.GroupDef);

        var sn = new Symm(n);
        var gi = (n - 2).Range(3).Select(i => sn[(1, 2, i)]).ToArray();
        return Group.Generate($"Alt{n}", sn, gi);
    }

    public static ConcreteGroup<Perm> Dihedral(int n)
    {
        var m = (n % 2) == 0 ? 1 : 2;
        var sn = new Sn(n);
        var an = Enumerable.Range(1, n).ToArray();
        var a2 = Enumerable.Range(m, n / 2).Select(i => (Tuple2Array)(i, n + m - i)).ToArray();
        var cn = sn.Cycle(an);
        var c2 = sn.ComposesCycles(a2);
        var d2n = Group.Generate($"D{2 * n}", sn, c2, cn);
        return d2n;
    }

    public static ConcreteGroup<Perm> PermGroup(string name, int n, params ValueType[] generators)
    {
        var sn = new Sn(n);
        var gi = generators.Select(g => sn.ComposesCycles(Tuple2Array.ComplexTuples(g))).ToArray();
        return Group.Generate(name, sn, gi);
    }

    public static ConcreteGroup<Perm> PermGroup(int n, params ValueType[] generators) => PermGroup("G", n, generators);

    public static ConcreteGroup<Ep<ZnInt>> Abelian(string name, params int[] seq)
    {
        return Product.GpGenerate(name, seq.Select(i => new Cn(i)).Cast<IGroup<ZnInt>>().ToArray());
    }

    public static ConcreteGroup<Ep<ZnInt>> Abelian(params int[] seq)
    {
        return Product.GpGenerate(seq.Select(i => new Cn(i)).Cast<IGroup<ZnInt>>().ToArray());
    }

    public static WordGroup AbelianWg(string name, params int[] seq)
    {
        if (seq.Length > 5 || seq.Min() <= 1)
            throw new GroupException(GroupExceptionType.GroupDef);
        
        var n = seq.Length.Range();
        var gens = n.Select(i => ((char)('a' + i)).ToString()).ToArray();
        var allCombs = n.SelectMany(i => n.Where(j => j > i).Select(j => (gens[i], gens[j]))).ToArray();
        var relators = gens.Select((g, i) => $"{g}{seq[i]}").Concat(allCombs.Select(e => $"{e.Item1}{e.Item2}={e.Item2}{e.Item1}"))
            .Glue(", ");
        
        return new WordGroup(name, relators);
    }

    public static WordGroup AbelianWg(params int[] seq) => AbelianWg(seq.Glue(" x ", "C{0}"), seq);

    public static SemiDirectProduct<ZnInt, ZnInt> DihedralSdp(int n)
    {
        var cn = new Cn(n);
        var autCn = Group.AutBase(cn);
        var y = autCn[(cn[1], cn[n - 1])];
        var aut = Group.Generate(autCn, y);

        var c2 = new Cn(2);
        var pMap = Group.PartialMap((c2[1], y));
        var theta = Group.Hom(c2, Group.HomomorphismMap(c2, aut, pMap));
        return Group.SemiDirectProd($"D{2 * n}", cn, theta, c2);
    }

    public static WordGroup DihedralWg(int n) => new($"D{2 * n}", $"a{n}, b2, abab");

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
                var wg = new WordGroup($"Frob({m},{n},{r})", $"a{m}, b{n}, b-1ab = a{r}");
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
        var autCm = Group.AutBase(cm);
        var g1 = autCm[(cm[1], cm[r])];
        var aut = Group.Generate(autCm, g1);
        var pMap = Group.PartialMap((cn[0], autCm.Neutral()), (cn[1], g1));
        var theta = Group.Hom(cn, Group.HomomorphismMap(cn, aut, pMap));
        var name = IntExt.Gcd(m, n * (r - 1)) == 1 ? $"Frob({m},{n},{r})" : $"MtCyc({m},{n},{r})";
        return Group.SemiDirectProd(name, cm, theta, cn);
    }

    public static int[] FrobeniusGetR(int m, int n)
    {
        return m.Range().Where(r => IntExt.Gcd(m, n * (r - 1)) == 1 && IntExt.PowMod(r, n, m) == 1).ToArray();
    }

    public static List<ConcreteGroup<Ep2<ZnInt, ZnInt>>> FrobeniusSdp(int order)
    {
        var ms = IntExt.Dividors(order).Where(d => d > 1 && d % 2 == 1).ToArray();

        List<ConcreteGroup<Ep2<ZnInt, ZnInt>>> all = new();
        foreach (var m in ms)
        {
            var n = order / m;
            var rs = FrobeniusGetR(m, n);
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
            return Group.Generate("Q8", gl, gl[2, 2, 2, 1], gl[0, 1, 2, 0]);
        }

        if (Math.Abs(m - 4.0) < 1e-5)
        {
            var gl = new GL(2, 7);
            return Group.Generate("Q16", gl, gl[0, 1, 6, 3], gl[6, 1, 5, 1]);
        }

        if (Math.Abs(m - 5.0) < 1e-5)
        {
            var gl = new GL(2, 17);
            return Group.Generate("Q32", gl, gl[14, 0, 0, 11], gl[0, 16, 1, 0]);
        }

        if (Math.Abs(m - 6.0) < 1e-5)
        {
            var gl = new GL(2, 31);
            return Group.Generate("Q64", gl, gl[16, 12, 12, 11], gl[0, 1, 30, 0]);
        }

        if (Math.Abs(m - 7.0) < 1e-5)
        {
            var gl = new GL(2, 193);
            return Group.Generate("Q128", gl, gl[78, 38, 155, 78], gl[71, 13, 13, 122]);
        }

        throw new GroupException(GroupExceptionType.GroupDef);
    }

    public enum DiCyclicGroupType
    {
        Quaternions,
        Even,
        Odd
    }

    public static DiCyclicGroupType GetDiCyclicType(int n)
    {
        var k = IntExt.PrimesDecomposition(n).Count(i => i == 2);
        var m = n / (1 << k);

        if (k == 0)
            return DiCyclicGroupType.Odd;

        if (m != 1)
            return DiCyclicGroupType.Even;

        return DiCyclicGroupType.Quaternions;
    }

    public static dynamic DiCyclicSdp(int n)
    {
        var k = IntExt.PrimesDecomposition(n).Count(i => i == 2);

        if (k != 0)
        {
            var m = n / (1 << k);
            var q = Quaternion(1 << (k + 2));
            if (m == 1)
                return q;

            var cm = new Cn(m);
            var autCm = Group.AutBase(cm);
            var a0 = autCm[(cm[1], cm[m - 1])];
            var aut0 = Group.Generate(autCm, a0);

            var qi = q.GetGenerators().OrderBy(e => q.ElementsOrders[e]).ToArray();
            var pMap0 = Group.PartialMap((qi[1], aut0.Neutral()), (qi[0], a0));
            var theta0 = Group.Hom(q, Group.HomomorphismMap(q, aut0, pMap0));
            return Group.SemiDirectProd($"Dic{n}", cm, theta0, q);
        }

        var cn = new Cn(n);
        var autCn = Group.AutBase(cn);
        var a1 = autCn[(cn[1], cn[n - 1])];
        var aut1 = Group.Generate(autCn, a1);

        var c4 = new Cn(4);
        var pMap1 = Group.PartialMap((c4[1], a1));
        var theta1 = Group.Hom(c4, Group.HomomorphismMap(c4, aut1, pMap1));
        return Group.SemiDirectProd($"Dic{n}", cn, theta1, c4);
    }

    public static WordGroup SemiDihedral(int n)
    {
        var n1 = 1 << (n - 1);
        var n2 = (1 << (n - 2)) - 1;
        return new WordGroup($"QD{n1 * 2}", $"a{n1}, b2, bab = a{n2}");
    }

    public static SemiDirectProduct<ZnInt, ZnInt> SemiDihedralSdp(int n)
    {
        var n1 = 1 << (n - 1);
        var n2 = (1 << (n - 2)) - 1;

        var cn = new Cn(n1);
        var c2 = new Cn(2);
        var autCn = Group.AutBase(cn);
        var y = autCn[(cn[1], cn[n2])];
        var aut = Group.Generate(autCn, y);
        var pMap = Group.PartialMap((c2[1], y));
        var theta = Group.Hom(c2, Group.HomomorphismMap(c2, aut, pMap));
        return Group.SemiDirectProd($"QD{n1 * 2}", cn, theta, c2);
    }

    public static WordGroup ModularMax(int n)
    {
        var n1 = 1 << (n - 1);
        var n2 = (1 << (n - 2)) + 1;
        return new WordGroup($"QD{n1 * 2}", $"a{n1}, b2, bab = a{n2}");
    }

    public static SemiDirectProduct<ZnInt, ZnInt> ModularMaxSdp(int n)
    {
        var n1 = 1 << (n - 1);
        var n2 = (1 << (n - 2)) + 1;

        var cn = new Cn(n1);
        var c2 = new Cn(2);
        var autCn = Group.AutBase(cn);
        var y = autCn[(cn[1], cn[n2])];
        var aut = Group.Generate(autCn, y);
        var pMap = Group.PartialMap((c2[1], y));
        var theta = Group.Hom(c2, Group.HomomorphismMap(c2, aut, pMap));
        return Group.SemiDirectProd($"QD{n1 * 2}", cn, theta, c2);
    }

    public static WordGroup WordGroup(string relators) => new WordGroup(relators);
    public static WordGroup WordGroup(string name, string relators) => new WordGroup(name, relators);

    public static ConcreteGroup<EPoly<ZnInt>> Galois(int q, char x = 'x')
    {
        var fq = new Fq(q, x);
        return Group.Generate(fq, fq[x]);
    }

    public static KAutGroup<K> KAutGroup<K>(KPoly<K> P) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new KAutGroup<K>(P);
    }

    public static EPoly<ZnInt> FqX(int q, char x = 'x') => new Fq(q, x)[x];

    public static NthRootQ NthRootQ(int n) => new(n);
}