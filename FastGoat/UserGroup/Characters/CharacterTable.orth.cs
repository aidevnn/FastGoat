using System.ComponentModel.DataAnnotations;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.Characters;

public partial class CharacterTable<T> where T : struct, IElt<T>
{
    public void PrefillOrthogonality()
    {
        var clIdx = Classes.Count.Range();
        var Ocl = Classes.ToDictionary(e => e, e => Classes.GetClassStabx(e).Count());
        foreach (var g in Classes)
        {
            var sum = Cnf.CnfZero;
            foreach (var i in clIdx)
            {
                var c = AllCharacters[i][g];
                if (!c.HasValue)
                    continue;

                sum += c.Value!.Pow(2);
            }

            // if (sum.Simplify().Equals(Ocl[g] * Cnf.CnfOne))
            if ((sum - Ocl[g]).IsZero())
            {
                foreach (var i in clIdx)
                {
                    var chiMap = AllCharacters[i].Map.ToDictionary(kv => kv.Key, kv => kv.Value);
                    var c = chiMap[g];
                    if (!c.HasValue)
                    {
                        chiMap[g] = Cnf.CnfZero;
                        AllCharacters[i] = new(Classes, chiMap);
                    }
                }
            }
        }
    }


    public void SolveOrthogonality()
    {
        PrefillOrthogonality();
        var nbOcl = AllCharacters.Count(chi => !chi.HasAllValues);
        if (nbOcl == 0)
            return;
        if (nbOcl > 3)
        {
            Console.WriteLine($"Missing {nbOcl} characters or more");
            // return;
        }

        var cells = new List<(int, int)>();
        for (int i = 0; i < NbClasses; i++)
        {
            for (int j = 0; j < NbClasses; j++)
            {
                if (!AllCharacters[i][Classes.GetRepresentative(j)].HasValue)
                    cells.Add((i, j));
            }
        }

        var xis = Ring.Polynomial(Cnf.CnfZero, MonomOrder.Lex, (cells.Count + 1, "x"));
        var mapCells = xis.SkipLast(1).Select((e, i) => (e, i)).ToDictionary(e => e.e, e => cells[e.i]);
        var mapSymb = mapCells.ToDictionary(e => e.Value, e => e.Key);
        var mapInd = mapCells.ToDictionary(e => e.Key.ExtractIndeterminate, e => e.Value);
        var xz = xis.Last();
        var mat = new Polynomial<Cnf, Xi>[NbClasses, NbClasses];

        foreach (var gi in Classes)
        {
            var i = Classes.GetIndex(gi);
            for (int j = 0; j < NbClasses; j++)
            {
                var c = AllCharacters[j][gi];
                if (!c.HasValue)
                    mat[j, i] = mapSymb[(j, i)];
                else
                    mat[j, i] = c.Value * xz.One;
            }
        }

        var rg = NbClasses.Range();
        var Ocl = Classes.ToDictionary(e => e, e => Classes.GetClassStabx(e).Count());
        Console.WriteLine();
        Console.WriteLine(Ring.Matrix2String(mat));
        Console.WriteLine();
        
        var allCombs = rg.SelectMany(i => rg.Where(j => j > i).Select(j => (i, j))).ToArray();
        var ggi = Gr.Select(e => (g: Classes.GetIndex(e), gi: Classes.GetIndex(Gr.Invert(e)))).ToArray();
        // var clggi = rg.ToDictionary(i => i, i => IndexesTint[CClasses.GetRepresentative(Gr.Invert(Repr[i]))]);
        var clggi = rg.ToDictionary(i => i, i => Classes.GetIndex(Gr.Invert(Classes.GetRepresentative(i))));
        
        var eqs = new List<Polynomial<Cnf, Xi>>();

        // Regular character
        eqs.AddRange(rg.Skip(1).Select(j => rg.Aggregate(xz.Zero, (sum, i) => sum + mat[i, 0] * mat[i, j])));
        
        // All i, Sum[g](Xi(g)Xi(g^−1))= |G|
        eqs.AddRange(rg.Select(i =>
            ggi.Aggregate(xz.Zero, (sum, kp) => sum + mat[i, kp.g] * mat[i, kp.gi]) -
            (xz.One * Gr.Count())));
        
        // All i <> j, Sum[g](Xi(g)Xj(g^−1))=  0
        eqs.AddRange(allCombs.Select(e =>
            ggi.Aggregate(xz.Zero,
                (sum, kp) => sum + mat[e.i, kp.g] * mat[e.j, kp.gi])));
        
        // All g, h in Cl(g), Sum[r](Xr(g)Xr(h^−1))= |Cl(g)|
        eqs.AddRange(clggi.Select(kp =>
            rg.Aggregate(xz.Zero,
                (sum, r) => sum + mat[r, kp.Key] * mat[r, kp.Value]) -
            (xz.One * Ocl[Classes.GetRepresentative(kp.Key)])));
        // All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1))=  0
        eqs.AddRange(
            allCombs.Select(e => rg.Aggregate(xz.Zero,
                (sum, r) => sum + mat[r, e.i] * mat[r, clggi[e.j]])));
        
        eqs.RemoveAll(e => e.IsZero());
        eqs.Distinct().OrderBy(eq => eq.Degree).ThenBy(eq => eq.Coefs.Count).Println("System");
        
        // var mapCnf = mapSymb.ToDictionary(e => e.Value.ExtractIndeterminate, e => Cnf.Nth(NbStabx[e.Key.Item2]));
        var mapCnf = mapSymb.ToDictionary(e => e.Value.ExtractIndeterminate,
            e => Cnf.Nth(Classes.GetClassStabx(e.Key.Item2).Count()));
        
        var solDegreeOne = SolveDegreOne(eqs.Where(eq => eq.Degree == 1).ToArray(), mapCnf, xz);
        var subsEq = eqs.Select(eq => solDegreeOne.Aggregate(eq, (eq0, s) => eq0.Substitute(s.Value, s.Key)))
            .Where(eq => !eq.IsZero()).Distinct().ToArray();
        subsEq.Println("Substitute Eqs");
        var (subs, subsEq2) = SubstitutionDegreeOne(subsEq);
        subs.Println();
        subsEq2.Println("Substitute Eqs2");
        var redEqs = Ring.ReducedGrobnerBasis(subsEq2);
        
        redEqs.Println("Reduced System");
        Console.WriteLine();
        var allSolutions = SolveSystem(solDegreeOne, redEqs, mapCnf, xz)
            .Select(sol => ReverseSubstitution(sol, subs))
            .ToArray();
        foreach (var solution in allSolutions)
        {
            solution.Println("Solution");
            Console.WriteLine();
            // break;
        }

        var firstSol = allSolutions[0];
        foreach (var (xi, cnf) in firstSol)
        {
            var (i, j) = mapInd[xi];
            var chiMap = AllCharacters[i].Map.ToDictionary(kv => kv.Key, kv => kv.Value);
            chiMap[Classes.GetRepresentative(j)] = cnf;
            AllCharacters[i] = new(Classes, chiMap);
            // Cells[i + 3, j + 2] = CnfCells[i, j] = new CnfCell(cnf);
        }
    }

    private Dictionary<Xi, Cnf> SolveDegreOne(Polynomial<Cnf, Xi>[] eqs, Dictionary<Xi, Cnf> mapCnf,
        Polynomial<Cnf, Xi> xz)
    {
        var oneIndeterminate = eqs.Where(eq => eq.Degree == 1 && eq.NbIndeterminates == 1)
            .OrderBy(eq => eq.Coefs.Keys.Max(m => m.Degree))
            .ToArray();

        var mapSol = new Dictionary<Xi, Cnf>();
        foreach (var p in oneIndeterminate)
        {
            var ind = p.ExtractIndeterminate;
            var cnfOrd = mapCnf[ind];
            mapSol[ind] = Solve(p, ind, cnfOrd)[0];
        }

        return mapSol;
    }

    private (List<(Xi, Polynomial<Cnf, Xi>)>, Polynomial<Cnf, Xi>[]) SubstitutionDegreeOne(Polynomial<Cnf, Xi>[] eqs)
    {
        var eqs0 = eqs.ToArray();
        var lt = new List<(Xi, Polynomial<Cnf, Xi>)>();
        while (eqs0.Any(eq => eq.Degree == 1))
        {
            var eq = eqs0.Where(eq => eq.Degree == 1).OrderBy(eq => eq.Coefs.Count).First(eq => eq.Degree == 1);
            var xi = eq.ExtractAllIndeterminates.Order().First();
            var ci = eq[new(eq.Indeterminates, xi)];
            var yi = -eq / ci + eq.X(xi);
            lt.Add((xi, yi));
            Console.WriteLine(new { eq, xi, yi });
            eqs0 = eqs0.Select(eq0 => eq0.Substitute(yi, xi)).Where(eq0 => !eq0.IsZero()).Distinct().ToArray();
        }

        return (lt, eqs0);
    }

    private Dictionary<Xi, Cnf> ReverseSubstitution(Dictionary<Xi, Cnf> solutions, List<(Xi, Polynomial<Cnf, Xi>)> subs)
    {
        var mapSol = solutions.ToDictionary(kv => kv.Key, kv => kv.Value);
        var subs0 = subs.ToList();
        while (true)
        {
            var s0 = subs0.Where(e => e.Item2.ExtractAllIndeterminates.ToHashSet().IsSubsetOf(mapSol.Keys)).ToArray();
            if (s0.Length == 0)
                break;
            
            var dic0 = s0.Select(e =>
                    (e.Item1, e.Item2.Substitute(mapSol.Select(f => (f.Key, f.Value * e.Item2.One)).ToList())))
                .ToDictionary(e => e.Item1, e => e.Item2.ConstTerm);

            foreach (var (k, v) in dic0)
                mapSol[k] = v;

            var xis = s0.Select(e => e.Item1).ToHashSet();
            subs0.RemoveAll(e => xis.Contains(e.Item1));
        }

        return mapSol;
    }

    private IEnumerable<Dictionary<Xi, Cnf>> SolveSystem(Dictionary<Xi, Cnf> solutions, Polynomial<Cnf, Xi>[] eqs,
        Dictionary<Xi, Cnf> mapCnf, Polynomial<Cnf, Xi> xz)
    {
        var oneIndeterminate = eqs.Where(eq => eq.NbIndeterminates == 1).OrderBy(eq => eq.Coefs.Keys.Max(m => m.Degree))
            .ToArray();
        // eqs.Println("Reduced System");
        // oneIndeterminate.Println("One Indeterminate");
        // Console.WriteLine();

        var mapSol = new Dictionary<Xi, Cnf[]>();
        foreach (var p in oneIndeterminate)
        {
            var ind = p.ExtractIndeterminate;
            var cnfOrd = mapCnf[ind];
            mapSol[ind] = Solve(p, ind, cnfOrd);
        }

        // mapSol.Select(e => $"{e.Key} in {{{e.Value.Glue("; ")}}}").Println("Solutions One Indeterminates");

        var allPos = mapSol.SelectMany(e => e.Value.Select(c => (e.Key, c * xz.One))).GroupBy(e => e.Item1)
            .Select(e => e.ToArray()).ToArray();

        foreach (var pos in allPos.MultiLoop())
        {
            var pos0 = pos.ToArray();
            pos0.Println("Possibility");
            var subsEq = eqs.Select(eq => pos0.Aggregate(eq, (eq0, s) => eq0.Substitute(s.Item2, s.Item1)))
                .Where(eq => !eq.IsZero()).ToArray();

            var nsol = new Dictionary<Xi, Cnf>(solutions);
            foreach (var kp in pos0)
                nsol.Add(kp.Key, kp.Item2.Coefs.Values.First());

            if (subsEq.Length == 0)
            {
                // Console.WriteLine("######### End");
                yield return new Dictionary<Xi, Cnf>(nsol);
            }
            else if (subsEq.All(eq => eq.NbIndeterminates != 0))
            {
                // subsEq.Println("Substituate System");
                if (subsEq.ToHashSet().SetEquals(eqs))
                {
                    yield return nsol;
                    yield break;
                }

                foreach (var nsol0 in SolveSystem(nsol, subsEq, mapCnf, xz))
                {
                    yield return nsol0;
                }
            }
            else
            {
                subsEq.Println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Rejected System");
                Console.WriteLine();
            }
        }
    }

    private Cnf[] Solve(Polynomial<Cnf, Xi> p, Xi ind, Cnf c)
    {
        List<Cnf> algSol = new();
        var f0 = p.ToKPoly(ind);
        if (f0.Degree == 1)
        {
            var sol = -f0[0] / f0[1];
            // Console.WriteLine($"Eq {f} = 0 with solutions {{{sol}}}");
            algSol.Add(sol * c.One);
            return algSol.ToArray();
        }

        var x = FG.QPoly();
        if (f0.Coefs.Any(c0 => c0.E.Poly.Degree > 0))
        {
            throw new Exception($"############# Faillure to solve {p} = 0");
        }

        var f = p.ToKPoly(ind).Coefs.Select((c0, i) => c0.E[0] * x.Pow(i)).Aggregate(x.Zero, (sum, xi) => sum + xi);
        // var facts0 = PolynomialFactorizationPart2.FirrZ(f, details: true);
        var facts0 = IntFactorisation.FirrZ(f);
        var degreeOne = facts0.Where(fi => fi.Degree == 1).ToArray();
        var sols = degreeOne.Select(fi => -fi[0] / fi[1]).ToArray();

        var irrs = facts0.Where(fi => fi.Degree > 1).ToArray();
        if (irrs.Length == 0)
        {
            // Console.WriteLine($"Eq {f} = 0 with solutions {{{sols.Glue("; ")}}}");
            algSol.AddRange(sols.Select(s => c.One * s));
            return algSol.ToArray();
        }

        algSol.AddRange(sols.Select(s => s * c));
        foreach (var fi in irrs)
        {
            var (X, a) = FG.EPolyXc(c.E.F, 'a');
            var P = fi.Substitute(X);
            Console.WriteLine($"Factors of {P} in splitting field {a.F} of Q({c})[x]");
            // var facts = AlgebraicFactorization.AlgebraicFactors(P, details: true);
            var facts = IntFactorisation.AlgebraicFactors(P);
            var degreeOneAlg = facts.Where(fj => fj.Degree == 1).ToArray();
            var irrsAlg = facts.Where(fj => fj.Degree > 1).ToArray();
            if (irrsAlg.Length > 0)
            {
                throw new Exception($"############# Faillure to solve {p} = 0");
            }

            algSol.AddRange(degreeOneAlg.Select(fj =>
                (-fj[0] / fj[1]).Poly.Coefs.Select((k, i) => k * c.Pow(i)).Aggregate(c.Zero, (sum, ci) => sum + ci)));
        }

        // Console.WriteLine($"Eq {f} = 0 with solutions {{{algSol.Glue("; ")}}}");
        return algSol.ToArray();
    }
}