using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;

namespace FastGoat.UserGroup;

public static partial class FG
{
    private static Dictionary<int, KPoly<Rational>> CyclotomicPolynomials { get; }
    private static Dictionary<int, CnfBasis> CnfBasisMap { get; }

    static FG()
    {
        var x = QPoly('X');
        CyclotomicPolynomials = new() { [1] = x - 1 };
        CnfBasisMap = new();
        
        allIds = GroupExt.DBGap.Select(txt => new IdGroup(txt)).GroupBy(e => e.Order).ToDictionary(e => e.Key, e => e.Order().ToArray());
        nbSubGroupsDetails = allIds.ToDictionary(
            e => e.Key,
            e => e.Value.GroupBy(a => a.Infos).ToDictionary(a => a.Key, a => a.Count()));
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

    public static KPoly<K> CyclotomicPolynomial<K>(int k, K scalar) where K : struct, IFieldElt<K>, IRingElt<K>, IElt<K>
    {
        var P = CyclotomicPolynomial(k);
        var one = scalar.One;
        return new(P.x, scalar, P.Coefs.Select(c => (int)c.Num * one).ToArray());
    }

    public static KPoly<Rational>[] CyclotomicPolynomialsSequence() => CyclotomicPolynomials.Values.Order().ToArray();

    public static AbelianDirectSum<T> AbelianDirectSum<T>(this ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        return new(gr);
    }
    
    public static CnfBasis CnfBasis(int ord)
    {
        if (CnfBasisMap.TryGetValue(ord, out var infos))
            return infos;

        return CnfBasisMap[ord] = new CnfBasis(ord);
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

    public static ConcreteGroup<Mat> DihedralGL2p(int n)
    {
        var p = IntExt.Primes10000.First(p => (p - 1) % n == 0);
        var a0 = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(p, n);
        var ord_n = new ZnInt(p, a0);
        var gl = new GL(2, p);
        var m0 = gl[ord_n.K, 0, 0, ord_n.Inv().K];
        var m1 = gl[0, 1, 1, 0];
        return Group.Generate($"D{2 * n}", gl, m0, m1);
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

    public static ConcreteGroup<ZnInt> UnInt(int n) =>
        Group.MulGroup($"U{n}i", IntExt.Coprimes(n).Select(j => new ZnInt(n, j)).ToArray());

    public static ConcreteGroup<Ep<ZnInt>> Abelian(params int[] seq)
    {
        return Product.GpGenerate(seq.Select(i => new Cn(i)).Cast<IGroup<ZnInt>>().ToArray());
    }

    public static ConcreteGroup<Perm> AbelianPerm(params int[] seq)
    {
        var n = seq.Sum();
        var sn = new Sn(n);
        var gens = new List<Perm>();
        for (int i = 0; i < seq.Length; i++)
        {
            var n0 = seq.Take(i).Sum() + 1;
            gens.Add(sn.Cycle(seq[i].Range(n0)));
        }

        var name = seq.Glue(" x ", "C{0}");
        return Group.Generate(name, sn, gens.ToArray());
    }

    public static ConcreteGroup<Ep<ZnInt>> ElementaryAbelian(int q)
    {
        var dec = IntExt.PrimesDec(q);
        if (dec.Count != 1)
            throw new GroupException(GroupExceptionType.GroupDef);

        var (p, n) = dec.First();
        return Abelian(Enumerable.Repeat(p, n).ToArray());
    }

    public static List<ConcreteGroup<Ep<ZnInt>>> AllAbelianGroupsOfOrder(int k)
    {
        if (k == 1)
            return new() { Abelian(1) };
        
        var dec = IntExt.PrimesDec(k);
        return dec.Select(e => IntExt.Partitions32[e.Value].Select(l => l.Select(i => e.Key.Pow(i)).ToArray())).MultiLoop()
            .Select(l => Abelian(l.SelectMany(i => i).OrderDescending().ToArray())).ToList();
    }

    public static WordGroup AbelianWg(string name, params int[] seq)
    {
        if (seq.Length > 6 || seq.Min() < 1)
            throw new GroupException(GroupExceptionType.GroupDef);

        var n = seq.Length.Range();
        var gens = n.Select(i => ((char)('a' + i)).ToString()).ToArray();
        var allCombs = n.SelectMany(i => n.Where(j => j > i).Select(j => (gens[i], gens[j]))).ToArray();
        var relators = gens.Select((g, i) => $"{g}{seq[i]}").Concat(allCombs.Select(e => $"{e.Item1}{e.Item2}={e.Item2}{e.Item1}"))
            .Glue(", ");

        return new WordGroup(name, relators);
    }

    public static WordGroup AbelianWg(params int[] seq) => AbelianWg(seq.Glue(" x ", "C{0}"), seq);

    public static WordGroup ElementaryAbelianWg(int q)
    {
        var dec = IntExt.PrimesDec(q);
        if (dec.Count != 1)
            throw new GroupException(GroupExceptionType.GroupDef);

        var (p, n) = dec.First();
        return AbelianWg(Enumerable.Repeat(p, n).ToArray());
    }

    public static List<WordGroup> AllAbelianGroupsOfOrderWg(int k)
    {
        if (k == 1)
            return new() { AbelianWg(1) };
        
        var dec = IntExt.PrimesDec(k);
        return dec.Select(e => IntExt.Partitions32[e.Value].Select(l => l.Select(i => e.Key.Pow(i)).ToArray())).MultiLoop()
            .Select(l => AbelianWg(l.SelectMany(i => i).OrderDescending().ToArray())).ToList();
    }

    public static WordGroup AlternateWG(int n)
    {
        if (n < 2 || n > 6)
            throw new GroupException(GroupExceptionType.GroupDef);

        var rg = (n - 2).Range();
        string s0(int e) => $"{(char)('a' + e)}";
        var r1 = rg.Select(e => $"{s0(e)}3").Glue(", ");
        var r2 = rg.Grid2D().Where(e => e.t1 < e.t2).Select(e => $"{s0(e.t1)}{s0(e.t2)}{s0(e.t1)}{s0(e.t2)}").Glue(", ");
        return WordGroup($"A{n}", $"{r1}, {r2}");
    }

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

    public static WordGroup MetaCyclicSdpWg(int m, int n, int r)
    {
        if (IntExt.PowMod(r, n, m) != 1 || IntExt.Gcd(r, m) != 1)
            throw new GroupException(GroupExceptionType.GroupDef);
        
        var name = IntExt.Gcd(m, n * (r - 1)) == 1 ? $"Frob({m},{n},{r})" : $"MtCyc({m},{n},{r})";
        return WordGroup(name, $"a{m}, b{n}, b-1ab = a{r}");
    }

    public static ConcreteGroup<Ep2<ZnInt, ZnInt>> MetaCyclicSdp(int m, int n, int r)
    {
        if (IntExt.PowMod(r, n, m) != 1 || IntExt.Gcd(r, m) != 1)
            throw new GroupException(GroupExceptionType.GroupDef);

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

    public static int[] MetaCyclicSdpGetR(int m, int n)
    {
        var rs = IntExt.SolveAll_k_pow_m_equal_one_mod_n(m, n);
        var nCoprimes = IntExt.Coprimes(n).ToArray();
        return rs.Select(r => nCoprimes.Select(k => IntExt.PowMod(r, k, m)).ToHashSet())
            .Distinct(new SetEquality<int>())
            .Select(e => e.Where(c => c != 1).Min())
            .ToArray();
    }

    public static int[] FrobeniusGetR(int m, int n)
    {
        return MetaCyclicSdpGetR(m, n).Where(r => IntExt.Gcd(m, n * (r - 1)) == 1).ToArray();
    }

    public static List<WordGroup> MetaCyclicSdpWg(int o)
    {
        return IntExt.Dividors(o).Where(d => d > 1 && d % 2 == 1)
            .SelectMany(m => MetaCyclicSdpGetR(m, o / m).Select(r => (m, n: o / m, r)))
            .Select(e => MetaCyclicSdpWg(e.m, e.n, e.r))
            .ToList();
    }

    public static List<WordGroup> Frobenius(int order)
    {
        return IntExt.Dividors(order).Where(d => d > 1 && d % 2 == 1)
            .SelectMany(m => FrobeniusGetR(m, order / m).Select(r => (m, n: order / m, r)))
            .Select(e => MetaCyclicSdpWg(e.m, e.n, e.r))
            .ToList();
    }

    public static List<ConcreteGroup<Ep2<ZnInt, ZnInt>>> FrobeniusSdp(int order)
    {
        return IntExt.Dividors(order).Where(d => d > 1 && d % 2 == 1)
            .SelectMany(m => FrobeniusGetR(m, order / m).Select(r => (m, n: order / m, r)))
            .Select(e => MetaCyclicSdp(e.m, e.n, e.r))
            .ToList();
    }

    public static List<ConcreteGroup<Ep2<ZnInt, ZnInt>>> MetaCyclicSdp(int order)
    {
        return IntExt.Dividors(order).Where(d => d > 1 && d % 2 == 1)
            .SelectMany(m => MetaCyclicSdpGetR(m, order / m).Select(r => (m, n: order / m, r)))
            .Select(e => MetaCyclicSdp(e.m, e.n, e.r))
            .ToList();
    }

    public static WordGroup DiCyclic(int n) => new(int.IsPow2(4 * n) ? $"Q{4 * n}" : $"Dic{n}", $"a{n} = b2, b2 = abab");

    public static ConcreteGroup<Mat> Quaternion(int k)
    {
        if (int.IsPow2(k))
            return DicyclicGL2p(k / 4);

        throw new GroupException(GroupExceptionType.GroupDef);
    }

    public static WordGroup QuaternionWg(int k)
    {
        if (int.IsPow2(k))
            return DiCyclic(k / 4);

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

    public static ConcreteGroup<Mat> DicyclicGL2p(int m)
    {
        var (p, a0) = IntExt.Primes10000.Select(p => (p, IntExt.Solve_k_pow_m_equal_one_mod_n_strict(p, 2 * m)))
            .First(e => e.Item2 != -1);
        var a = new ZnInt(p, a0);
        var ai = a.Inv();
        var gl = new GL(2, p);
        var Am = gl[a.K, 0, 0, ai.K];
        var B = gl[0, 1, -1, 0];
        var Dic_m = Group.Generate($"Dic{m}", gl, Am, B);
        if (int.IsPow2(Dic_m.Count()))
            Dic_m.Name = $"Q{m * 4}";
        
        return Dic_m;
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
    
    public static ConcreteGroup<Mat> SemiDihedralGL2p(int n)
    {
        var n1 = 1 << (n - 1);
        var n2 = (1 << (n - 2)) - 1;
    
        var p = IntExt.Primes10000.First(p => (p - 1) % n1 == 0);
        var gl = new GL(2, p);
        var ordns = IntExt.SolveAll_k_pow_m_equal_one_mod_n_strict(p, n1).ToArray();
        var a0 = ordns.First(e => IntExt.PowMod(IntExt.PowMod(e, n2, p), n2, p) == e);
        var a1 = IntExt.PowMod(a0, n2, p);

        var m0 = gl[a0, 0, 0, a1];
        var m1 = gl[0, 1, 1, 0];
    
        return Group.Generate($"QD{2 * n1}", gl, m0, m1);
    }

    public static WordGroup ModularMax(int n)
    {
        var n1 = 1 << (n - 1);
        var n2 = (1 << (n - 2)) + 1;
        return new WordGroup($"MM{n1 * 2}", $"a{n1}, b2, bab = a{n2}");
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
        return Group.SemiDirectProd($"MM{n1 * 2}", cn, theta, c2);
    }

    
    public static ConcreteGroup<Mat> ModularMaxGL2p(int n)
    {
        var n1 = 1 << (n - 1);
        var n2 = (1 << (n - 2)) + 1;
    
        var p = IntExt.Primes10000.First(p => (p - 1) % n1 == 0);
        var gl = new GL(2, p);
        var ordns = IntExt.SolveAll_k_pow_m_equal_one_mod_n_strict(p, n1).ToArray();
        var a0 = ordns.First(e => IntExt.PowMod(IntExt.PowMod(e, n2, p), n2, p) == e);
        var a1 = IntExt.PowMod(a0, n2, p);

        var m0 = gl[a0, 0, 0, a1];
        var m1 = gl[0, 1, 1, 0];
    
        return Group.Generate($"MM{2 * n1}", gl, m0, m1);
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
    public static NthRootFq NthRootFq(int n, int q) => new(n, q);
    public static NthRootFp NthRootFp(int n, int p) => new(n, p);

    public static NthRootFq CycloFq(int q)
    {
        var dec = IntExt.PrimesDec(q);
        if (dec.Count != 1)
            throw new();

        var p = dec.First().Key;
        return new(q - 1, p);
    }

    public static (GFp gf, ConcreteGroup<EPoly<ZnInt>>) GFp(KPoly<ZnInt> f)
    {
        var gf = new GFp(f);
        var gens = f.P.Range().Where(e => !f.Substitute(e).IsZero()).Select(e => gf.X - e).ToArray();
        return (gf, Group.Generate(gf, gens));
    }
}