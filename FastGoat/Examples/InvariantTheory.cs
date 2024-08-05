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
        GlobalStopWatch.Restart();
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

    static Polynomial<K, Xi>[] Reynolds<K>(KMatrix<K>[] G, Polynomial<K, Xi>[] xi)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (M, N) = G[0].Dim;
        if (M != N || G.Select(g => g.Dim).Distinct().Count() != 1)
            throw new();

        if (xi.Length != M || xi.Distinct().Count() != xi.Length)
            throw new();

        var f0 = xi[0].One;
        if (f0.P != 0)
            throw new($"Field characteristic must be 0 but egal {f0.P}");

        var coefs = G.Length.Range(1).Select(i => xi.Aggregate((a0, a1) => a0 + a1).Pow(i))
            .Aggregate((a0, a1) => a0 + a1);

        var mn = coefs.Coefs.Keys.Select(m => m.ToPolynomial(f0)).Order().ToArray();
        var facts = mn.Select(f => (f, Rf: Reynolds(G, f, xi))).ToArray();
        facts.Println(
            $"Reynolds Table {xi.Select((xk, k) => $"{xk}^i{k}").Glue(" * ")} when {xi.Select((xk, k) => $"i{k}").Glue(" + ")}<={G.Length}");
        Console.WriteLine($"Nb eqs:{mn.Length}");
        Console.WriteLine();

        return facts.Where(e => !e.Rf.IsZero()).Select(e => e.Rf).DistinctBy(p => p.Monic()).Order().ToArray();
    }

    static bool IdealContains<K>(List<Polynomial<K, Xi>> I, Polynomial<K, Xi> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var allXis = I.Append(f).SelectMany(p => p.ExtractAllIndeterminates).Distinct().ToArray();
        var xi = f.Indeterminates.First(xi => !allXis.Contains(xi));
        var s = I.Append(f * f.X(xi) - 1).ToArray();
        var red = Ring.ReducedGrobnerBasis(s);
        return red.Any(p => p.Degree == 0);
    }

    static void Coefs<K>(int o, List<Polynomial<K, Xi>> coefs, Polynomial<K, Xi> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var k = o / f.Degree;
        var lf = k.Range(1).Select(i => f.Pow(i)).ToArray();
        var lc = coefs.SelectMany(g1 => lf.Select(g2 => g1 * g2).Where(g3 => g3.Degree <= o)).ToArray();
        var tmp = coefs.Concat(lf).Concat(lc).Order().ToList();
        coefs.Clear();
        coefs.AddRange(tmp);
    }

    static (Polynomial<K, Xi>[] inv, (Polynomial<K, Xi> p, Polynomial<K, Xi> mod)[] mods)
        InvariantGLnK<K>(ConcreteGroup<KMatrix<K>> G, KPoly<Rational> MolienSerie, MonomOrder order)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var A = G.Neutral();
        var (M, N) = A.Dim;
        if (M != N)
            throw new();

        var molien = (MolienSerie.Degree + 1).Range().ToDictionary(i => i, i => (int)MolienSerie.Coefs[i].Num);
        var xi0 = Ring.Polynomial(A.KZero, order, M.Range().Select(i => $"x{i}").ToArray());
        var Rfs0 = Reynolds(G.ToArray(), xi0);
        
        var one = xi0[0].One;
        one.Indeterminates.ExtendAppend(Rfs0.Length.Range().Select(i => new Xi($"u{i}")).ToArray());
        var xi = one.Indeterminates.Select(u => u.ToPolynomial(one)).ToArray();
        var ui = xi.Skip(M).ToArray();
        
        Rfs0.Println("System");

        var invariants = new List<Polynomial<K, Xi>>();
        var modulos = new List<(Polynomial<K, Xi> p, Polynomial<K, Xi> mod)>();
        var set = Rfs0.ToHashSet();
        var coefs = new List<Polynomial<K, Xi>>();
        foreach (var p in Rfs0.Where(p => molien[p.Degree] != 0))
        {
            var contains = IdealContains(invariants, p);
            if (!contains)
            {
                invariants.Add(p);
                Console.WriteLine("Invariant added");
                Console.WriteLine(p);
                set.Remove(p);
                Coefs(G.Count(), coefs, p);
                Console.WriteLine();
            }
        }

        var invs = Ring.ReducedGrobnerBasis(invariants.Zip(ui).Select(f => f.First.Monic() - f.Second).ToArray());
        foreach (var p in set.Where(p => molien[p.Degree] != 0))
        {
            if (coefs.Where(f => f.Degree <= p.Degree).Any(f => p.Div(f).rem.IsZero()))
                continue;

            var uf = p.Indeterminates.Where(u => u.xi.Contains('u'))
                .Except(invs.Append(p).SelectMany(e => e.ExtractAllIndeterminates))
                .First();
            var pu = p.Monic() - p.X(uf);
            var sys = invs.Append(pu).ToArray();
            sys.Println("System Modulo");
            var red = Ring.ReducedGrobnerBasis(sys);
            var mods = red.Select(f => (f, xi: f.ExtractAllIndeterminates))
                .Where(e => e.xi.Contains(uf) && e.xi.All(x => x.xi.Contains('u')))
                .Select(e => e.f)
                .ToArray();
            Console.WriteLine();
            if (mods.Length == 0)
            {
                Console.WriteLine("Continue");
                Console.WriteLine();
            }
            
            var mod = mods.Order().FirstOrDefault(f => f.Coefs.Keys.Any(mn => mn.Contains(uf)), p.One);
            Console.WriteLine($"Modulo:{mod}");
            if (mod.Coefs.Keys.Count(mn => mn.Equals(uf)) == 1 && mod.Coefs.Keys.Count(mn => mn.Contains(uf)) == 1)
            {
                Console.WriteLine("Continue");
                Console.WriteLine();
                continue;
            }

            Console.WriteLine("Add modulo");
            Console.WriteLine();
            modulos.Add((pu, mod));
            break;
        }

        var invsf = invariants.Zip(ui).Select(f => f.First.Monic() - f.Second).ToArray();
        invsf.Println("Invariant generators");
        Console.WriteLine();
        if (modulos.Count > 0)
        {
            Console.WriteLine("Generators and modulos");
            foreach (var (p, mod) in modulos)
            {
                Console.WriteLine($"    {p}");
                Console.WriteLine($"    {mod}");
                Console.WriteLine();
            }
        }

        return (invsf, modulos.ToArray());
    }

    // TODO fix taylor serie
    static (Polynomial<K, Xi>[] inv, (Polynomial<K, Xi> p, Polynomial<K, Xi> mod)[] mods)
        InvariantGLnK<K>(ConcreteGroup<KMatrix<K>> G, MonomOrder order = MonomOrder.GrLex, bool molien = false)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var A = G.Neutral();
        if (A.P != 0)
            throw new($"Field characteristic must be 0 but egal {A.P}");

        GlobalStopWatch.AddLap();
        var x = FG.QPoly();
        var serie = molien ? MolienSum(G).serie : (x + 1).Pow(G.Count());
        var r = InvariantGLnK(G, serie, order);
        GlobalStopWatch.Show();
        Console.WriteLine();

        return r;
    }
    
    static (Polynomial<EPoly<Rational>, Xi>[] inv, (Polynomial<EPoly<Rational>, Xi> p,
        Polynomial<EPoly<Rational>, Xi>mod)[] mods) Invariant_Cn_GL2R_EPoly(int n)
    {
        var c = Cnf.Nth(n);
        var n0 = IntExt.Lcm(c.Re.N, c.Im.N);
        var (cos, sin) = (c.Re.ToCnfN(n0).E, c.Im.ToCnfN(n0).E);
        var gl = FG.GLnK($"Q(ξ{n0})", 2, cos);
        var A = gl[cos, -sin, sin, cos];
        var G = Group.Generate($"C{n}", gl, A);
        DisplayGroup.HeadElements(G);
        return InvariantGLnK(G);
    }

    static (Polynomial<Cnf, Xi>[] inv, (Polynomial<Cnf, Xi> p,
        Polynomial<Cnf, Xi> mod)[] mods) Invariant_Cn_GL2R_Cnf(int n)
    {
        var c = Cnf.Nth(n);
        var (cos, sin) = (c.Re, c.Im);
        var gl = FG.GLnK($"Cnf", 2, cos);
        var A = gl[cos, -sin, sin, cos];
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

        var derDets = dets.Select(e => (c: e.KOne, nm: e.One, dnm: e)).ToArray();
        var serie = dets[0].Zero;

        for (int i = 0; i <= G.Count(); i++)
        {
            var fi = i == 0 ? sum.KOne : i * sum.KOne;
            for (int k = 0; k < derDets.Length; ++k)
            {
                var (c, nm, dnm) = derDets[k];
                var s0 = c * nm[0] / dnm[0];
                serie += s0 * fi.Inv() * nm.X.Pow(i);
                var derNm = nm.Derivative * dnm - nm * dnm.Derivative;
                var derDnm = dnm.Pow(2);
                var gcd = Ring.Gcd(derNm.Monic, derDnm.Monic).Monic;
                (derNm, derDnm) = (derNm / gcd, derDnm / gcd);
                derDets[k] = (c * fi.Inv() * derNm.LT / derDnm.LT, derNm / derNm.LT, derDnm / derDnm.LT);
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

    public static void Example_C4_GLnC()
    {
        {
            var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
            var A = gl[Cnf.I, 0, 0, -Cnf.I];
            var G = Group.Generate("C4", gl, A);
            DisplayGroup.HeadElements(G);
            InvariantGLnK(G);
        }

        {
            var gl = FG.GLnK("R", 2, Rational.KOne());
            var A = gl[0, -1, 1, 0];
            var G = Group.Generate("C4", gl, A);
            DisplayGroup.HeadElements(G);
            InvariantGLnK(G);
        }

        {
            var gl = FG.GLnK("R", 3, Rational.KOne());
            var A = gl[0, 1, 0, -1, 0, 0, 0, 0, -1];
            var G = Group.Generate("C4", gl, A);
            DisplayGroup.HeadElements(G);
            InvariantGLnK(G);
            
            // Invariant generators
            //     x2^2 - u0
            //     x0^2 + x1^2 - u1
            //     x0^2*x1^2 - u2
            // 
            // Generators and modulos
            //     x0*x1*x2 - u3
            //     u0*u2 - u3^2
            // 
            // #  Time:203ms
        }
    }

    public static void Example_C3_GLnR()
    {
        Invariant_Cn_GL2R_EPoly(3);

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

            InvariantGLnK(G); // Time:1.735s
        }
    }

    public static void Example_C6_GL2R()
    {
        var gl = FG.GLnK("R", 2, Rational.KOne());
        var A = gl[0, 1, -1, 1];
        var G = Group.Generate("C6", gl, A);
        DisplayGroup.HeadElements(G);
        Console.WriteLine();

        InvariantGLnK(G);
    }

    public static void Example_Invariant_GL2R_EPoly_CyclicGroups()
    {
        GlobalStopWatch.AddLap();
        Invariant_Cn_GL2R_EPoly(2);
        Invariant_Cn_GL2R_EPoly(3);
        Invariant_Cn_GL2R_EPoly(4);
        Invariant_Cn_GL2R_EPoly(5);
        Invariant_Cn_GL2R_EPoly(6);
        Invariant_Cn_GL2R_EPoly(7);
        Invariant_Cn_GL2R_EPoly(8);
        Invariant_Cn_GL2R_EPoly(9);
        GlobalStopWatch.Show("End");
    }

    public static void Example_Invariant_GL2R_Cnf_CyclicGroups()
    {
        GlobalStopWatch.AddLap();
        Invariant_Cn_GL2R_Cnf(2);
        Invariant_Cn_GL2R_Cnf(3);
        Invariant_Cn_GL2R_Cnf(4);
        Invariant_Cn_GL2R_Cnf(5);
        Invariant_Cn_GL2R_Cnf(6);
        Invariant_Cn_GL2R_Cnf(7);
        Invariant_Cn_GL2R_Cnf(8);
        GlobalStopWatch.Show("End");
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
        var c = Cnf.Nth(4);
        var A = gl[c, 0, 0, c.Inv()];
        var B = gl[0, 1, 1, 0];
        var G = Group.Generate("D4", gl, A, B);
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

    #region Groups of order 12
    public static void Example_C12_GL3C()
    {
        var gl = FG.GLnK("Cnf", 2, Cnf.CnfOne);
        var c = Cnf.Nth(12);
        var A = gl[c.Re, -c.Im, c.Im, c.Re];
        var G = Group.Generate("C12", gl, A);
        DisplayGroup.HeadElements(G);
        Console.WriteLine();

        InvariantGLnK(G);
    }

    public static void Example_C2xC6_GL3C()
    {
        var gl = FG.GLnK("Cnf", 3, Cnf.CnfOne);
        var c = Cnf.Nth(6);
        var A = gl[c.Re, -c.Im, 0, c.Im, c.Re, 0, 0, 0, 1];
        var B = gl[1, 0, 0, 0, 1, 0, 0, 0, -1];
        var G = Group.Generate("C2 x C6", gl, A, B);
        DisplayGroup.HeadElements(G);
        Console.WriteLine();

        InvariantGLnK(G);
    }

    public static void Example_A4_GL2R()
    {
        var gl = FG.GLnK("R", 3, Rational.KOne());
        var A = gl[0, 1, 0, 0, 0, 1, 1, 0, 0];
        var B = gl[1, 0, 0, 0, -1, 0, 0, 0, -1];
        var G = Group.Generate("A4", gl, A, B);
        DisplayGroup.HeadElements(G);
        DisplayGroup.AreIsomorphics(G, FG.Alternate(4));
        Console.WriteLine();

        InvariantGLnK(G); // Time:10.993s
    }

    public static void Example_Dic3_GL2C()
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

    public static void Example_D12_GL2C()
    {
        var gl = FG.GLnK("Cnf", 3, Cnf.CnfOne);
        var c = Cnf.Nth(3);
        var A = gl[c, 0, 0, 0, c.Pow(2), 0, 0, 0, -1];
        var B = gl[0, 1, 0, 1, 0, 0, 0, 0, 1];
        var G = Group.Generate("D12", gl, A, B);
        DisplayGroup.HeadElements(G);
        DisplayGroup.AreIsomorphics(G, FG.Dihedral(6));
        Console.WriteLine();

        InvariantGLnK(G); // Time:32.421s
    }

    public static void Example_D12_GL2R()
    {
        var gl = FG.GLnK("Cnf", 3, Cnf.CnfOne);
        var c = Cnf.Nth(6);
        var A = gl[c.Re, -c.Im, 0, c.Im, c.Re, 0, 0, 0, 1];
        var B = gl[1, 0, 0, 0, -1, 0, 0, 0, -1];
        var G = Group.Generate("D12", gl, A, B);
        DisplayGroup.HeadElements(G);
        DisplayGroup.AreIsomorphics(G, FG.Dihedral(6));
        Console.WriteLine();

        InvariantGLnK(G); // Time:54.042s
    }

    public static void Example_Dn_GL2C()
    {
        GlobalStopWatch.AddLap();
        for (int n = 3; n < 9; n++)
        {
            var c = Cnf.Nth(n);
            var gl = FG.GLnK("Cnf", 2, c);
            var A = gl[c, 0, 0, c.Inv()];
            var B = gl[0, 1, 1, 0];
            var G = Group.Generate($"D{2 * n}", gl, A, B);
            DisplayGroup.HeadElements(G);
            DisplayGroup.AreIsomorphics(G, FG.Dihedral(n));
            Console.WriteLine();

            InvariantGLnK(G);
        }

        GlobalStopWatch.Show("End"); // Time:1m25s

        // D2n = <[ξn, 0, 0, ξn^-1], [0, 1, 1, 0]>
        // C[x, y]D2n = C[xy, x^n + y^n] TODO proof
    }
    #endregion

    public static void Example_Dn_GL2R()
    {
        GlobalStopWatch.AddLap();
        for (int n = 3; n < 9; n++)
        {
            var c = Cnf.Nth(n);
            var gl = FG.GLnK("Cnf", 2, c);
            var A = gl[c.Re, -c.Im, c.Im, c.Re];
            var B = gl[0, 1, 1, 0];
            var G = Group.Generate($"D{2 * n}", gl, A, B);
            DisplayGroup.HeadElements(G);
            DisplayGroup.AreIsomorphics(G, FG.Dihedral(n));
            Console.WriteLine();

            InvariantGLnK(G);
        }

        GlobalStopWatch.Show("End"); // Time:3m17s

        // D2n = <[cos(2π/n), -sin(2π/n), sin(2π/n), cos(2π/n)], [0, 1, 1, 0]>
        // C[x, y]D2n TODO Invariant ring
    }
    
    public static void Example_M4sdp4_3_GL3C()
    {
        var I = Cnf.I;
        var gl = FG.GLnK("Cnf", 3, I);
        var A = gl[1, 0, 0, 0, I, 0, 0, 0, -I];
        var B = gl[I, 0, 0, 0, 0, 1, 0, 1, 0];
        var G = Group.Generate("M(4x:4)3", gl, A, B);
        DisplayGroup.HeadElements(G);
        DisplayGroup.AreIsomorphics(G, FG.MetaCyclicSdpWg(4, 4, 3));
        Console.WriteLine();

        InvariantGLnK(G);
        
        // Invariant generators
        //     x1*x2 - u0
        //     x1^4 + x2^4 - u1
        //     x0^4 - u2
        // 
        // Generators and modulos
        //     x0^2*x1^4 - x0^2*x2^4 - u3
        //     u0^4*u2 - 1/4*u1^2*u2 + 1/4*u3^2
        // 
        // #  Time:16.015s
    }
    
    public static void Example_Q16_GL2C()
    {
        var c = Cnf.Nth(8);
        var gl = FG.GLnK("Cnf", 2, c);
        var A = gl[0, 1, -1, 0];
        var B = gl[c.Pow(3), 0, 0, -c];
        var G = Group.Generate("Q16", gl, A, B);
        DisplayGroup.HeadElements(G);
        DisplayGroup.AreIsomorphics(G, FG.Quaternion(16));
        Console.WriteLine();

        InvariantGLnK(G);
        
        // Invariant generators
        //     x0^2*x1^2 - u0
        //     x0^8 + x1^8 - u1
        // 
        // Generators and modulos
        //     x0^9*x1 - x0*x1^9 - u2
        //     u0^5 - 1/4*u0*u1^2 + 1/4*u2^2
        // 
        // #  Time:3.524s
    }
    
    public static void Example_C2xC2sdpC4_GL3C()
    {
        var c = Cnf.Nth(4);
        var gl = FG.GLnK("Cnf", 3, c);
        var A = gl[1, 0, 0, 0, 1, 0, 0, 0, -1];
        var B = gl[1, 0, 0, 0, -1, 0, 0, 0, -1];
        var C = gl[c, 0, 0, 0, 0, 1, 0, 1, 0];
        var G = Group.Generate("(C2 x C2) x: C4", gl, A, B, C);
        DisplayGroup.HeadElements(G);
        DisplayGroup.AreIsomorphics(G, FG.WordGroup("(C2 x C2) x: C4", "a4, b2, c2, bcbc, caca-1, abca-1b"));
        Console.WriteLine();

        InvariantGLnK(G);
        // Invariant generators
        //     x1^2 + x2^2 - u0
        //     x1^2*x2^2 - u1
        //     x0^4 - u2
        // 
        // Generators and modulos
        //     x0^2*x1^2 - x0^2*x2^2 - u3
        //     u0^2*u2 - 4*u1*u2 - u3^2
        // #  Time:13.617s
    }
}