using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

/*
 * Ideals, Varieties, and Algorithms. chap 7 Invariant Theory of Finite Groups page 345, 385
 */
public static class InvariantTheory
{
    static InvariantTheory()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        ZnInt.Display = ZnDisplay.Signed;
    }

    static Polynomial<K, Xi> ApplyF<K>(Polynomial<K, Xi> f, KMatrix<K> A, params Polynomial<K, Xi>[] xi)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (xi.Length != A.M || xi.Distinct().Count() != xi.Length)
            throw new();

        var nb = xi.Length;
        var monoms = xi.Select(x => x.ExtractIndeterminate).ToArray();

        var kone = f.One;
        var A0 = A.Select(a => a * kone).ToKMatrix(A.M);
        var X = monoms.Select(m => m.ToPolynomial(f.One)).ToKMatrix(A.M);
        var A0X = A0 * X;
        var dicoSubs = nb.Range().ToDictionary(i => monoms[i], i => A0X[i, 0]);
        return f.Substitute(dicoSubs.Select(e => (e.Key, e.Value)).ToList());
    }

    static Polynomial<K, Xi> Reynolds<K>(KMatrix<K>[] G, Polynomial<K, Xi> f, Polynomial<K, Xi>[] xi)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (f.P != 0)
            throw new($"Field characteristic must be 0 but egal {f.P}");

        var gi = (f.KOne * G.Length).Inv();
        return gi * G.Aggregate(f.Zero, (acc, g) => acc + ApplyF(f, g, xi));
    }

    static Polynomial<K, Xi>[] Reynolds<K>(KMatrix<K>[] G,KPoly<Rational> MolienSerie, Polynomial<K, Xi>[] xi)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (M, N) = G[0].Dim;
        if (M != N || G.Select(g => g.Dim).Distinct().Count() != 1)
            throw new();

        if (xi.Length != M || xi.Distinct().Count() != xi.Length)
            throw new();

        var degrees = MolienSerie.Coefs.Select((e, k) => (e, k)).Where(e => !e.e.IsZero()).Select(e => e.k).ToHashSet();
        var f0 = xi[0].One;
        if (f0.P != 0)
            throw new($"Field characteristic must be 0 but egal {f0.P}");
        
        var coefs = G.Length.Range(1).Select(i => xi.Aggregate((a0, a1) => a0 + a1).Pow(i))
            .Aggregate((a0, a1) => a0 + a1);

        var mn = coefs.Coefs.Keys.Where(m => degrees.Contains(m.Degree)).Select(m => m.ToPolynomial(f0)).Order()
            .ToArray();
        var facts = mn.Select(f => (f, Reynolds(G, f, xi))).ToArray();
        facts.Println(
            $"Reynolds Table {xi.Select((xk, k) => $"{xk}^i{k}").Glue(" * ")} when {xi.Select((xk, k) => $"i{k}").Glue(" + ")}<={G.Length}");
        Console.WriteLine($"Nb eqs:{mn.Length}");
        Console.WriteLine();

        var res = facts.Where(e => !e.Item2.IsZero()).Select(e => e.Item2 / e.Item2.LeadingDetails.lc).Distinct()
            .Order().ToArray();
        return res;
    }

    static (Polynomial<K, Xi>[] invGens, Polynomial<K, Xi>[] inv, Polynomial<K, Xi>[] rfs)
        InvariantGLnK<K>(ConcreteGroup<KMatrix<K>> G, KPoly<Rational> MolienSerie, MonomOrder order)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var A = G.Neutral();
        var (M, N) = A.Dim;
        if (M != N)
            throw new();

        var xi0 = Ring.Polynomial(A.KZero, order, M.Range().Select(i => $"x{i}").ToArray());
        var Rfs0 = Reynolds(G.ToArray(), MolienSerie, xi0);
        var one = xi0[0].One;
        one.Indeterminates.ExtendAppend(Rfs0.Length.Range().Select(i => new Xi($"u{i}")).ToArray());
        var xi = one.Indeterminates.Select(u => u.ToPolynomial(one)).ToArray();
        var Rfs1 = Rfs0.Select((f, i) => f - xi[i + M]).ToArray();
        Rfs1.Println("System");

        var dico = Rfs1.GroupBy(e => e.Degree).ToDictionary(e => e.Key, e => e.ToArray());
        var degree = dico.Keys.Order().ToArray();
        var rfs1 = new List<Polynomial<K, Xi>>();
        foreach (var k in degree)
        {
            var part = rfs1.Union(dico[k]).ToArray();
            part.Println("Sys partial");
            var red0 = SimplifyLoop(Ring.ReducedGrobnerBasis(part));
            red0.Println("Red partial");
            Console.WriteLine();
            rfs1.Clear();
            rfs1.AddRange(red0);
        }

        var red1 = SimplifyLoop(rfs1.ToArray());
        rfs1.Clear();
        var ui = red1.SelectMany(s => s.ExtractAllIndeterminates).Distinct().Where(e => e.xi.Contains('u')).Order()
            .ToArray();
        var zip = ui.Zip(xi.Skip(M)).ToArray();
        rfs1.AddRange(red1.Select(p => zip.Aggregate(p, (acc, e) => acc.Substitute(e.Second, e.First))));
        rfs1.Println("Simplifyed generators");
        Console.WriteLine();
        
        var idl = rfs1.Where(f => f.ExtractAllIndeterminates.All(u => u.xi.Contains('u'))).ToArray();
        var inv = Rfs1.Where(f => f.ExtractAllIndeterminates.Any(u => ui.Contains(u)))
            .Select(p => zip.Aggregate(p, (acc, e) => acc.Substitute(e.Second, e.First)))
            .Concat(idl).ToArray();
        inv.Println("Invariant generators");
        Console.WriteLine();
        
        return (inv.ToArray(), rfs1.ToArray(), Rfs1);
    }
    
    static (Polynomial<K, Xi>[] invGens, Polynomial<K, Xi>[] inv, Polynomial<K, Xi>[] rfs)
        InvariantGLnK<K>(ConcreteGroup<KMatrix<K>> G, MonomOrder order = MonomOrder.GrLex)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var A = G.Neutral();
        var (M, N) = A.Dim;
        var xi0 = Ring.Polynomial(Rational.KOne(), MonomOrder.GrLex, M.Range().Select(i => $"x{i}").ToArray());
        var all = G.Count().Range(1).ToDictionary(k => k, k => xi0.Aggregate((a, b) => a + b).Pow(k).Coefs.Count);
        var x = FG.QPoly();
        var serieMax = all.Select(e => e.Value * x.Pow(e.Key)).Aggregate((a, b) => a + b) + 1;
        if (A.P != 0)
            return InvariantGLnK(G, serieMax, order);
        
        var (sum, serie) = MolienSum(G);
        Console.WriteLine(new { serieMax });
        return InvariantGLnK(G, serie, order);
    }

    static Polynomial<K, Xi>[] Simplify<K>(Polynomial<K, Xi>[] sys) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var p0 = sys[0];
        var ind = sys[0].Indeterminates;
        var ui = sys.SelectMany(s => s.ExtractAllIndeterminates).Distinct().Where(e => e.xi.Contains('u')).Order()
            .ToArray();

        var xi = ind.Except(ui).ToArray();
        var subs = sys.Grid2D(ui).Where(e =>
                e.t1.Coefs.Keys.Count(e1 => e1[e.t2] == 1 && e1.Degree == 1) == 1 &&
                xi.All(x => !e.t1.ExtractAllIndeterminates.Contains(x)))
            .Select(e => (e.t2, -e.t1 / e.t1.Coefs[new Monom<Xi>(ind, e.t2, 1)] + e.t2.ToPolynomial(p0.One)))
            .GroupBy(e => e.t2).Select(e => (e.Key, e.First().Item2)).ToArray();

        if (subs.Length == 0)
            return sys;

        var sub = subs.MaxBy(e => e.Key);
        Console.WriteLine($"Sub : {sub}");
        var sys2 = sys.Select(s => s.Substitute(sub.Item2, sub.Key)).Where(e => !e.IsZero()).Distinct().ToArray();
        var sys3 = Ring.ReducedGrobnerBasis(sys2);
        return sys3;
    }

    static Polynomial<K, Xi>[] SimplifyLoop<K>(Polynomial<K, Xi>[] sys)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var sys0 = sys.ToArray();
        var sz = 0;
        while (sys0.Length != sz)
        {
            sz = sys0.Length;
            sys0 = Simplify(sys0);
        }

        return sys0;
    }

    static (GLn<EPoly<Rational>> gl, KMatrix<EPoly<Rational>> A) RealRotation_2Pi_over_k(int k)
    {
        var c = Cnf.Nth(k);
        var n0 = IntExt.Lcm(c.Re.N, c.Im.N);
        var (cos, sin) = (c.Re.ToCnfN(n0).E, c.Im.ToCnfN(n0).E);
        var gl = FG.GLnK($"Q(ξ{n0})", 2, cos);
        return (gl, gl[cos, -sin, sin, cos]);
    }

    static (GLn<Cnf> gl, KMatrix<Cnf> A) RealRotation2_2Pi_over_k(int k)
    {
        var c = Cnf.Nth(k);
        var (cos, sin) = (c.Re, c.Im);
        var gl = FG.GLnK($"Cnf", 2, cos);
        return (gl, gl[cos, -sin, sin, cos]);
    }

    static (Polynomial<EPoly<Rational>, Xi>[] invGens, Polynomial<EPoly<Rational>, Xi>[] inv,
        Polynomial<EPoly<Rational>, Xi>[] rfs) InvariantCn(int n)
    {
        var (gl, A) = RealRotation_2Pi_over_k(n);
        var G = Group.Generate($"C{n}", gl, A);
        DisplayGroup.HeadElements(G);
        return InvariantGLnK(G);
    }

    static (Polynomial<Cnf, Xi>[] invGens, Polynomial<Cnf, Xi>[] inv, Polynomial<Cnf, Xi>[] rfs) Invariant2Cn(int n)
    {
        var (gl, A) = RealRotation2_2Pi_over_k(n);
        var G = Group.Generate($"C{n}", gl, A);
        DisplayGroup.HeadElements(G);
        return InvariantGLnK(G);
    }

    static (EPolynomial<K> sum, KPoly<Rational> serie) MolienSum<K>(ConcreteGroup<KMatrix<K>> G)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var id = G.Neutral();
        if (id.P != 0)
            throw new($"Field characteristic must be 0 but egal {id.P}");
        
        var t = Ring.EPolynomial(id.KOne, MonomOrder.Lex, (1, "t"))[0];
        var sum = t.Zero;
        var dets = new List<KPoly<K>>();
        foreach (var m0 in G)
        {
            var m1 = m0.N.Range().Grid2D().Select(e => id[e.t1, e.t2] - t * m0[e.t1, e.t2]).ToKMatrix(m0.N);
            Console.WriteLine(m0);
            var det = m1.Det;
            dets.Add(det.Num.ToKPoly(det.Num.ExtractIndeterminate) * G.Count());
            var deti = det.Inv();
            Console.WriteLine(deti);
            Console.WriteLine();
            sum += deti;
        }

        sum /= G.Count();
        Console.WriteLine(new { sum });
        
        var derDets = dets.Select(e => (e.One, e)).ToArray();
        var serie = dets[0].Zero;

        for (int i = 0; i <= G.Count(); i++)
        {
            var fi = i == 0 ? sum.KOne : i.Range(1).Aggregate(sum.KOne, (a, b) => a * b);

            for (int k = 0; k < derDets.Length; ++k)
            {
                var (nm, dnm) = derDets[k];
                var s0 = nm[0] / dnm[0];
                serie += s0 * fi.Inv() * nm.X.Pow(i);
                var derNm = nm.Derivative * dnm - nm * dnm.Derivative;
                var derDnm = dnm.Pow(2);
                var gcd = Ring.Gcd(derNm.Monic, derDnm.Monic).Monic;
                derDets[k] = (derNm / gcd, derDnm / gcd);
            }

            Console.WriteLine(new { i, serie });
        }
        
        var T = FG.QPoly('t');
        var serie0 = serie.Coefs.Select((c, k) => int.Parse($"{c}") * T.Pow(k)).Aggregate((a, b) => a + b);
        return (sum, serie0);
    }

    public static void Example_Klein_GL2Q()
    {
        var gl = FG.GLnK("Q", 2, Rational.KOne());
        var A = gl[-1, 0, 0, 1];
        var B = gl[1, 0, 0, -1];
        var G = Group.Generate("V", gl, A, B);
        DisplayGroup.HeadElements(G);
        InvariantGLnK(G);
    }

    public static void Example_C4_GL2C()
    {
        var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
        var A = gl[Cnf.I, 0, 0, -Cnf.I];
        var G = Group.Generate("C4", gl, A);
        DisplayGroup.HeadElements(G);
        InvariantGLnK(G);
    }

    public static void Example_C3_GL2R()
    {
        InvariantCn(3);
        
        {
            var gl = FG.GLnK("R", 2, Rational.KOne());
            var A = gl[0, -1, 1, -1];
            var G = Group.Generate("C3", gl, A);
            DisplayGroup.HeadElements(G);
            Console.WriteLine();

            InvariantGLnK(G);
        }
        
        {
            var gl = FG.GLnK("R", 3, Rational.KOne());
            var A = gl[0, 1, 0, 0, 0, 1, 1, 0, 0];
            var G = Group.Generate("C3", gl, A);
            DisplayGroup.HeadElements(G);
            Console.WriteLine();

            InvariantGLnK(G); // ~4min
        }
    }

    public static void Example_Invariant_GL2K_CyclicGroups()
    {
        InvariantCn(2);
        InvariantCn(3);
        InvariantCn(4);
        InvariantCn(5);
        InvariantCn(6);
        InvariantCn(7);
        InvariantCn(8);
        InvariantCn(9);
    }

    public static void Example_Invariant2_GL2K_CyclicGroups()
    {
        Invariant2Cn(2);
        Invariant2Cn(3);
        Invariant2Cn(4);
        Invariant2Cn(5);
        Invariant2Cn(6);
        Invariant2Cn(7);
        Invariant2Cn(8);
    }

    public static void Example_C4xC2_GL2C()
    {
        var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
        var A = gl[Cnf.I, 0, 0, -Cnf.I];
        var B = gl[-Cnf.I, 0, 0, -Cnf.I];
        var G = Group.Generate("C4 x C2", gl, A, B);
        DisplayGroup.HeadElements(G);
        InvariantGLnK(G);
    }

    public static void Example_Q8_GL2C()
    {
        var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
        var A = gl[Cnf.I, 0, 0, -Cnf.I];
        var B = gl[0, -1, 1, 0];
        var G = Group.Generate("Q8", gl, A, B);
        DisplayGroup.HeadElements(G);
        DisplayGroup.AreIsomorphics(G, FG.Quaternion(8));
        Console.WriteLine();

        InvariantGLnK(G);
    }

    public static void Example_D8_GL2C()
    {
        var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
        var n = 4;
        var e = Cnf.Nth(n);
        var A = gl[e, 0, 0, e.Inv()];
        var B = gl[0, 1, 1, 0];
        var G = Group.Generate($"D{2 * n}", gl, A, B);
        DisplayGroup.HeadElements(G);
        DisplayGroup.AreIsomorphics(G, FG.Dihedral(4));
        Console.WriteLine();

        InvariantGLnK(G);
    }

    public static void Examples_MolienTheorem()
    {
        {
            var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
            var A = gl[Cnf.I, 0, 0, -Cnf.I];
            var B = gl[0, 1, 1, 0];
            var G = Group.Generate("D8", gl, A, B);
            DisplayGroup.HeadElements(G);
            DisplayGroup.AreIsomorphics(G, FG.Dihedral(4));
            Console.WriteLine();

            var (sum, serie) = MolienSum(G);
            Console.WriteLine($"MolienSum({G}) = {sum}");
            Console.WriteLine($"MolienSerie({G}) = {serie}");
            Console.WriteLine();
            // MolienSum(D8) = -1/(-t0^6 + t0^4 + t0^2 - 1)
            // MolienSerie(D8) = 3*t^8 + 2*t^6 + 2*t^4 + t^2 + 1
        }

        {
            var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
            var A = gl[Cnf.I, 0, 0, -Cnf.I];
            var B = gl[0, -1, 1, 0];
            var G = Group.Generate("Q8", gl, A, B);
            DisplayGroup.HeadElements(G);
            DisplayGroup.AreIsomorphics(G, FG.Quaternion(8));
            Console.WriteLine();

            var (sum, serie) = MolienSum(G);
            Console.WriteLine($"MolienSum({G}) = {sum}");
            Console.WriteLine($"MolienSerie({G}) = {serie}");
            Console.WriteLine();
            // MolienSum(Q8) = (4/3*t0^4 - 4/3*t0^2 + 4/3)/(4/3*t0^6 - 4/3*t0^4 - 4/3*t0^2 + 4/3)
            // MolienSerie(Q8) = 3*t^8 + t^6 + 2*t^4 + 1
        }
        
        {
            var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
            var c = Cnf.Nth(6);
            var A = gl[c, 0, 0, c.Inv()];
            var B = gl[0, 1, 1, 0];
            var G = Group.Generate("D12", gl, A, B);
            DisplayGroup.HeadElements(G);
            DisplayGroup.AreIsomorphics(G, FG.Dihedral(6));
            Console.WriteLine();

            var (sum, serie) = MolienSum(G);
            Console.WriteLine($"MolienSum({G}) = {sum}");
            Console.WriteLine($"MolienSerie({G}) = {serie}");
            Console.WriteLine();
            // MolienSum(D12) = -1/(-t0^8 + t0^6 + t0^2 - 1)
            // MolienSerie(D12) = 3*t^12 + 2*t^10 + 2*t^8 + 2*t^6 + t^4 + t^2 + 1
        }
    }

    public static void Example_C2xC2xC2_GL3R()
    {
        var gl = FG.GLnK("R", 3, Rational.KOne());
        var A = gl[-1, 0, 0, 0, 1, 0, 0, 0, 1];
        var B = gl[1, 0, 0, 0, -1, 0, 0, 0, 1];
        var C = gl[1, 0, 0, 0, 1, 0, 0, 0, -1];
        var G = Group.Generate("C2 x C2 x C2", gl, A, B, C);
        DisplayGroup.HeadElements(G);
        Console.WriteLine();

        InvariantGLnK(G);
    }

    public static void Example_Dic3_GL3C()
    {
        var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
        var c = Cnf.Nth(3);
        var A = gl[c, 0, 0, c.Inv()];
        var B = gl[0, -1, 1, 0];
        var G = Group.Generate("Dic3", gl, A, B);
        DisplayGroup.HeadElements(G);
        DisplayGroup.AreIsomorphics(G, FG.DiCyclic(3));
        Console.WriteLine();

        InvariantGLnK(G);
    }
    
    public static void Example_C3xC3_GL3C()
    {
        var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
        var c = Cnf.Nth(3);
        var A = gl[c, 0, 0, 1];
        var B = gl[1, 0, 0, c];
        var G = Group.Generate("C3 x C3", gl, A, B);
        DisplayGroup.HeadElements(G);
        Console.WriteLine();

        InvariantGLnK(G);
    }
}