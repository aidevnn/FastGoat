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

    static Polynomial<K, Xi>[] Reynolds<K>(KMatrix<K>[] G, params Polynomial<K, Xi>[] xi)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (M, N) = G[0].Dim;
        if (M != N || G.Select(g => g.Dim).Distinct().Count() != 1)
            throw new();

        if (xi.Length != M || xi.Distinct().Count() != xi.Length)
            throw new();

        var f0 = xi[0].One;
        var coefs = G.Length.Range(1).Select(i => xi.Aggregate((a0, a1) => a0 + a1).Pow(i))
            .Aggregate((a0, a1) => a0 + a1);

        var mn = coefs.Coefs.Keys.Select(m => m.ToPolynomial(f0)).Order().ToArray();
        var facts = mn.Select(f => (f, G.Aggregate(f0.Zero, (acc, g) => acc + ApplyF(f, g, xi)))).ToArray();
        facts.Println(
            $"Reynolds Table {xi.Select((xk, k) => $"{xk}^i{k}").Glue(" * ")} when {xi.Select((xk, k) => $"i{k}").Glue(" + ")}<={G.Length}");
        Console.WriteLine($"Nb eqs:{mn.Length}");
        Console.WriteLine();

        var res = facts.Where(e => !e.Item2.IsZero()).Select(e => e.Item2 / e.Item2.LeadingDetails.lc).Distinct()
            .Order().ToArray();
        return res;
    }

    static (Polynomial<K, Xi>[] invGens, Polynomial<K, Xi>[] inv, Polynomial<K, Xi>[] rfs)
        InvariantGLnK<K>(ConcreteGroup<KMatrix<K>> G, MonomOrder order = MonomOrder.GrLex)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var A = G.Neutral();
        var (M, N) = A.Dim;
        if (M != N)
            throw new();

        var xi0 = Ring.Polynomial(A.KZero, order, M.Range().Select(i => $"x{i}").ToArray());
        var Rfs0 = Reynolds(G.ToArray(), xi0);
        xi0[0].Indeterminates.ExtendAppend(Rfs0.Length.Range().Select(i => new Xi($"u{i}")).ToArray());
        var xi = xi0[0].Indeterminates.Select(u => u.ToPolynomial(xi0[0].One)).ToArray();
        var Rfs1 = Rfs0.Select((f, i) => f - xi[i + M]).ToArray();
        Rfs1.Println("System");
        var red = Ring.ReducedGrobnerBasis(Rfs1);
        red.Println("Reduced System");
        Console.WriteLine();
        var red2 = SimplifyLoop(red);
        red2.Println("Simplifyed generators");
        Console.WriteLine();
        var ui = red2.SelectMany(s => s.ExtractAllIndeterminates).Distinct().Where(e => e.xi.Contains('u')).Order()
            .ToArray();
        var idl = red2.Where(f => f.ExtractAllIndeterminates.All(u => ui.Contains(u))).ToArray();
        var invGens = new Polynomial<K, Xi>[0];
        if (idl.Length != 0)
        {
            var uis = idl.SelectMany(f => f.ExtractAllIndeterminates).Distinct().Order().ToArray();
            var inv = Rfs1.Where(f => f.ExtractAllIndeterminates.Distinct().Any(u => uis.Contains(u))).Concat(idl)
                .ToArray();
            inv.Println("Invariant generators");
            invGens = inv.ToArray();
            Console.WriteLine();
        }

        return (invGens, red2, Rfs1);
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
        return (gl, gl[cos, sin, -sin, cos]);
    }

    static (GLn<Cnf> gl, KMatrix<Cnf> A) RealRotation2_2Pi_over_k(int k)
    {
        var c = Cnf.Nth(k);
        var (cos, sin) = (c.Re, c.Im);
        var gl = FG.GLnK($"Cnf", 2, cos);
        return (gl, gl[cos, sin, -sin, cos]);
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

    static (EPolynomial<K>, Polynomial<K, Xi> serie) MolienSum<K>(ConcreteGroup<KMatrix<K>> G)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var id = G.Neutral();
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

        for (int i = 0; i < G.Count(); i++)
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
        
        return (sum, serie.Substitute(t.Num));
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

    public static void Example_Invariant_GL2K_CyclicGroups()
    {
        InvariantCn(2);
        InvariantCn(3);
        InvariantCn(4);
        InvariantCn(5);
        InvariantCn(6);
        InvariantCn(7);
        InvariantCn(8);
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

    public static void Example_C4xC2_GL2Z()
    {
        var gl = FG.GLnK("F5", 2, new ZnInt(5, 0));
        var A = gl[2, 0, 0, 3];
        var B = gl[3, 0, 0, 3];
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

    public static void Example_Q8_GL2Z()
    {
        var gl = FG.GLnK("F5", 2, new ZnInt(5, 0));
        var A = gl[2, 0, 0, -2];
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

    public static void Example_D12_GL2Z()
    {
        var gl = FG.GLnK("F7", 3, new ZnInt(7, 0));
        var A = gl[2, 0, 0, 0, 4, 0, 0, 0, 6];
        var B = gl[0, 1, 0, 1, 0, 0, 0, 0, 1];
        var G = Group.Generate("D12", gl, A, B);
        DisplayGroup.HeadElements(G);
        DisplayGroup.AreIsomorphics(G, FG.Dihedral(6));
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
            // MolienSum(D8) = -1/(-t^6 + t^4 + t^2 - 1)
            // MolienSerie(D8) = 2*t^6 + 2*t^4 + t^2 + 1
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
            // MolienSum(Q8) = (4/3*t^4 - 4/3*t^2 + 4/3)/(4/3*t^6 - 4/3*t^4 - 4/3*t^2 + 4/3)
            // MolienSerie(Q8) = t^6 + 2*t^4 + 1
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
            // MolienSum(D12) = -1/(-t^8 + t^6 + t^2 - 1)
            // MolienSerie(D12) = 2*t^10 + 2*t^8 + 2*t^6 + t^4 + t^2 + 1
        }
    }
}