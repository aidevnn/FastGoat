using System.ComponentModel.DataAnnotations;
using System.Numerics;
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

    public void SolveOrthogonality(params (int dim, int[] linIdx)[] infos)
    {
        PrefillOrthogonality();
        AllCharacters = AllCharacters.Order().ToArray();

        var todoChis = AllCharacters.Select((chi, k) => (chi, k)).Where(e => !e.chi.HasAllValues).ToArray();
        if (todoChis.Length == 0)
            return;

        var cells = new List<(int, int)>();
        for (int i = 0; i < NbClasses; i++)
        {
            var chi = AllCharacters[i];
            for (int j = 0; j < NbClasses; j++)
            {
                var e = Classes.GetRepresentative(j);
                if (!chi[e].HasValue)
                {
                    cells.Add((i, j));
                }
            }
        }

        var xis = Ring.Polynomial(Cnf.CnfZero, MonomOrder.Lex, (cells.Count + 1, "x"));
        var mapCells = cells.Take(cells.Count).Select((ij, k) => (ij, k)).ToDictionary(e => xis[e.k], e => e.ij);

        var mapSymb = mapCells.GroupBy(e => e.Value)
            .ToDictionary(e => e.Key, e => e.Select(e0 => e0.Key).Aggregate((a0, a1) => a0 + a1));
        var mapInd = mapCells.ToDictionary(e => e.Key.ExtractIndeterminate, e => e.Value);
        var xz = xis[cells.Count];
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

        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine();
            Console.WriteLine(Ring.Matrix2String(mat));
            Console.WriteLine();
        }

        var list = new List<Polynomial<Cnf, Xi>[]>();
        for (int i = 0; i < NbClasses; i++)
            list.Add(NbClasses.Range().Select(j => mat[i, j]).ToArray());

        var eqs = new List<Polynomial<Cnf, Xi>>();
        var firstIdx = todoChis.Max(e => e.k);

        var rmSeq = todoChis.GroupBy(e => e.chi[Gr.Neutral()]!.Value)
            .ToDictionary(l => (int)l.Key.E[0].Num, l => l.Select(e => e.k).ToArray());

        var allIdx = infos.SelectMany(l => l.linIdx).Distinct().ToArray();
        var linChis = AllCharacters.Where(chi => chi.IsLinear)
            .Select((chi, k) => (chi, k))
            .Where(e => allIdx.Contains(e.k))
            .ToDictionary(e => e.k, e => e.chi);

        var infosDim = infos.GroupBy(e => e.dim).ToDictionary(e => e.Key, e => e.Select(l => l.linIdx).ToArray());
        foreach (var (dim, linIdxs) in infosDim)
        {
            var rm = rmSeq[dim];
            int nb = 0;
            foreach (var linIdx in linIdxs)
            {
                var k0 = rm[nb];
                nb += linIdx.Length;
                foreach (var (idx, i) in linIdx.Select((idx, i) => (l0: idx, i)))
                {
                    var k1 = k0 + i;
                    var linChi = linChis[idx];
                    var pol = NbClasses.Range().Select(k => mat[k0, k] * linChi[Classes.GetRepresentative(k)]!.Value)
                        .ToArray();
                    foreach (var i0 in NbClasses.Range())
                        eqs.Add(list[k1][i0] - pol[i0]);

                    list[k1] = pol;
                }
            }
        }

        // Regular Character
        for (int i = 1; i < NbClasses; ++i)
        {
            var xi0 = mat[firstIdx, i];
            var xi1 = -NbClasses.Range()
                .Where(k => k != firstIdx)
                .Aggregate(xz.Zero, (sum, k) => sum + list[k][0] * list[k][i]);

            eqs.Add(xi0 * list[firstIdx][0] - xi1);
        }

        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine();
            Console.WriteLine(Ring.Matrix2String(Ring.Matrix(list.Count, list.SelectMany(l => l).ToArray())));
            Console.WriteLine();
        }

        var rg = NbClasses.Range();
        var Ocl = Classes.ToDictionary(e => e, e => Classes.GetClassStabx(e).Count());

        var allCombs = rg.SelectMany(i => rg.Where(j => j > i).Select(j => (i, j))).ToArray();
        var ggi = Gr.Select(e => (g: Classes.GetIndex(e), gi: Classes.GetIndex(Gr.Invert(e)))).ToArray();
        var clggi = rg.ToDictionary(i => i, i => Classes.GetIndex(Gr.Invert(Classes.GetRepresentative(i))));

        var mapXiDegDim = mapCells.ToDictionary(e => e.Key.ExtractIndeterminate,
            e => (Gr.ElementsOrders[Classes.GetRepresentative(e.Value.Item2)], AllCharacters[e.Value.Item1].Dim));

        // All i, Sum[g](Xi(g)Xi(g^−1))= |G|
        eqs.AddRange(rg.Select(i =>
            ggi.Aggregate(xz.Zero, (sum, kp) => sum + list[i][kp.g] * list[i][kp.gi]) - xz.One * Gr.Count()));

        // All i <> j, Sum[g](Xi(g)Xj(g^−1))=  0
        eqs.AddRange(allCombs.Select(e =>
            ggi.Aggregate(xz.Zero,
                (sum, kp) => sum + list[e.i][kp.g] * list[e.j][kp.gi])));

        // All g, h in Cl(g), Sum[r](Xr(g)Xr(h^−1))= |Cl(g)|
        eqs.AddRange(clggi.Select(kp =>
            rg.Aggregate(xz.Zero, (sum, r) => sum + list[r][kp.Key] * list[r][kp.Value]) -
            xz.One * Ocl[Classes.GetRepresentative(kp.Key)]));

        // All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1))=  0
        eqs.AddRange(
            allCombs.Select(e => rg.Aggregate(xz.Zero,
                (sum, r) => sum + list[r][e.i] * list[r][clggi[e.j]])));

        eqs.RemoveAll(e => e.IsZero());
        if (Logger.Level != LogLevel.Off)
            eqs.Distinct().OrderBy(eq => eq.Degree).ThenBy(eq => eq.Coefs.Count).Println("System");

        var solDegreeOne = SolveDegreOne(eqs.Where(eq => eq.Degree == 1).ToArray(), mapXiDegDim, xz);
        var subsEq = eqs.Select(eq => solDegreeOne.Aggregate(eq, (eq0, s) => eq0.Substitute(s.Value, s.Key)))
            .Where(eq => !eq.IsZero()).Distinct().ToArray();

        if (Logger.Level != LogLevel.Off)
        {
            solDegreeOne.Println(e=>$"{e.Key} <- {e.Value}", "Substitution1");
            subsEq.Println("Substitute Eqs");
        }

        var (subs, subsEq2) = SubstitutionDegreeOne(subsEq, mapXiDegDim);

        if (Logger.Level != LogLevel.Off)
        {
            subs.Println(e=>$"{e.Item1} <- {e.Item2}", "Substitution2");
            subsEq2.Println("Substitute Eqs2");
        }

        var redEqs = Ring.ReducedGrobnerBasis(subsEq2);

        if (Logger.Level != LogLevel.Off)
        {
            redEqs.Println("Reduced System");
            Console.WriteLine();
        }

        var allSolutions = SolveSystem(solDegreeOne, redEqs, mapXiDegDim, xz)
            .Select(sol => ReverseSubstitution(sol, subs));

        foreach (var sols in allSolutions.Where(sol => sol.Count + 1 == xis.Length))
        {
            var chis = AllCharacters.Select(chi => new Character<T>(Classes, chi.Map)).ToArray();
            foreach (var groupXis in sols.Where(e => mapInd.ContainsKey(e.Key)).GroupBy(e => mapInd[e.Key].Item1))
            {
                var chiMap = chis[groupXis.Key].Map.ToDictionary(kv => kv.Key, kv => kv.Value);
                foreach (var (xi, cnf) in groupXis)
                {
                    var j = mapInd[xi].Item2;
                    chiMap[Classes.GetRepresentative(j)] = cnf;
                }

                chis[groupXis.Key] = new Character<T>(Classes, chiMap);
                if (chis.All(chi => chi.HasAllValues))
                    break;
            }

            if (CheckExtSymmDecomposition(chis))
            {
                if (Logger.Level != LogLevel.Off)
                {
                    sols.OrderBy(e => e.Key).Println("Solutions");
                    Console.WriteLine();
                }

                foreach (var (_, idx) in todoChis)
                {
                    var chi = chis[idx];
                    var state = AddCharacter(chi);
                    if (state == AddCharacterState.TableFull)
                        return;
                }
            }
            
            if (Logger.Level != LogLevel.Off)
            {
                sols.OrderBy(e => e.Key).Println("Rejected");
                Console.WriteLine();
            }
        }
    }

    private Dictionary<Xi, Cnf> SolveDegreOne(Polynomial<Cnf, Xi>[] eqs, Dictionary<Xi, (int, int)> mapXiDegDim,
        Polynomial<Cnf, Xi> xz)
    {
        var oneIndeterminate = eqs.Where(eq => eq.Degree == 1 && eq.NbIndeterminates == 1)
            .OrderBy(eq => eq.Coefs.Keys.Max(m => m.Degree))
            .ToArray();

        var mapSol = new Dictionary<Xi, Cnf>();
        foreach (var p in oneIndeterminate)
        {
            var (ind, cs) = Solve(p, mapXiDegDim).First()[0];
            mapSol[ind] = cs;
        }

        return mapSol;
    }

    private (List<(Xi, Polynomial<Cnf, Xi>)>, Polynomial<Cnf, Xi>[]) SubstitutionDegreeOne(Polynomial<Cnf, Xi>[] eqs,
        Dictionary<Xi, (int deg, int dim)> mapXiDegDim)
    {
        var eqs0 = eqs.ToArray();
        var lt = new List<(Xi, Polynomial<Cnf, Xi>)>();
        while (eqs0.Any(eq => eq.Degree == 1))
        {
            var eq = eqs0.Where(eq => eq.Degree == 1)
                .OrderByDescending(eq => eq.ExtractAllIndeterminates.Max(xi => mapXiDegDim[xi].deg))
                .ThenBy(eq => eq.Coefs.Count)
                .First(eq => eq.Degree == 1);
            var xi = eq.ExtractAllIndeterminates.MaxBy(xi => mapXiDegDim[xi].deg);
            var ci = eq[new(eq.Indeterminates, xi)];
            var yi = -eq / ci + eq.X(xi);
            lt.Add((xi, yi));
            // Console.WriteLine(new { eq, xi, yi });
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
        Dictionary<Xi, (int deg, int dim)> mapXiDegDim, Polynomial<Cnf, Xi> xz)
    {
        var eqFirst = eqs.OrderBy(eq => eq.NbIndeterminates)
            .ThenBy(eq => eq.Degree)
            .ThenBy(eq => eq.ExtractAllIndeterminates.Select(xi => mapXiDegDim[xi])
                .Select(e => BigInteger.Pow(e.deg + 1, e.dim)).Aggregate((a0, a1) => a0 * a1))
            .ThenBy(eq => eq.Coefs.Count())
            .First();
        foreach (var pos in Solve(eqFirst, mapXiDegDim))
        {
            if (pos.Length == 0)
                yield break;

            var pos0 = pos.ToArray();
            if (Logger.Level != LogLevel.Off)
                pos0.Println("Possibility");

            var subsEq = eqs.Select(eq => pos0.Aggregate(eq, (eq0, s) => eq0.Substitute(s.Item2, s.Item1)))
                .Where(eq => !eq.IsZero()).ToArray();

            var nsol = new Dictionary<Xi, Cnf>(solutions);
            foreach (var kp in pos0)
                nsol.Add(kp.Item1, kp.Item2);

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

                foreach (var nsol0 in SolveSystem(nsol, subsEq, mapXiDegDim, xz))
                {
                    yield return nsol0;
                }
            }
            else
            {
                if (Logger.Level != LogLevel.Off)
                {
                    subsEq.Println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Rejected System");
                    Console.WriteLine();
                }
            }
        }
    }

    private IEnumerable<(Xi xi, Cnf)[]> Solve(Polynomial<Cnf, Xi> p, Dictionary<Xi, (int deg, int dim)> mapXi)
    {
        var xis = p.ExtractAllIndeterminates;
        if (xis.Length == 1 && p.Degree == 1)
        {
            List<Cnf> algSol = new();
            var ind = xis[0];
            var (deg, dim) = mapXi[ind];
            var c = Cnf.Nth(deg);
            var f0 = p.ToKPoly(ind);
            var sol = -f0[0] / f0[1];
            algSol.Add(sol * c.One);
            return algSol.Select(c0 => new[] { (ind, c0) });
        }

        if (xis.Length == 1 && p.Degree == 2)
        {
            var ind = xis[0];
            var deg0 = mapXi[ind].deg;
            var p0 = p.ToKPoly(ind);
            var coefs0 = p0.Coefs.ToArray();
            var deg1 = IntExt.Lcm(p0.Coefs.Select(c0 => c0.N).ToArray());

            var deg = IntExt.Lcm(deg0, deg1);
            var e = FG.CyclotomicEPoly(deg);
            var coefs = p0.ToCnfPoly(deg).Coefs.Select(c0 => c0.E).ToArray();

            var X = FG.KPoly('X', e);
            var P = coefs.Select((c0, i) => c0 * X.Pow(i)).Aggregate(X.Zero, (sum, xi) => sum + xi);

            var lvl = Logger.SetOff();
            var roots = IntFactorisation.AlgebraicRoots(P).Select(r => new Cnf(deg, r).Simplify())
                .OrderByDescending(r => r.E.Poly.Coefs.All(c0 => c0.Denom == 1)) // alg. int. first
                .ToArray();
            Logger.Level = lvl;

            if (Logger.Level != LogLevel.Off)
            {
                Console.WriteLine(
                    $"Factors of {P} in splitting field {e.F} of Q({e.F.x})[x] with {e.F.x} = {Cnf.Nth(deg)}");
                roots.Select(r => (ind, r, p.Substitute(r, ind), p.Substitute(r, ind).IsZero()))
                    .Println($"Nb possibilities:{roots.Length}");
            }

            if (roots.Length == 0)
                throw new();

            return roots.Select(r => new[] { (ind, r) });
        }

        if (xis.Length > 0)
        {
            var cidim = xis.Select(xi => (xi, mapXi[xi])).Select(e => (e.xi, e.Item2.deg, e.Item2.dim))
                .ToArray();

            var seqs = cidim.Select(e => (e.xi, e.deg, e.dim,
                    seq: e.deg.Range().Select(k => Cnf.Nth(e.deg).Pow(k).Simplify())
                        .Append(Cnf.CnfZero)
                        // .OrderBy(c => c) // uncomment to change ordering, for SL(2,3) log details and C2 x SL(2,3) counterexample
                        .OrderBy(c => c.IsZero()) // comment to change ordering
                        .ThenBy(c => c.N) // comment to change ordering
                        .ToHashSet()))
                .ToArray();

            var total = seqs.Select(e => BigInteger.Pow(e.deg + 1, e.dim)).Aggregate((a0, a1) => a0 * a1);

            if (Logger.Level != LogLevel.Off)
                seqs.Select(e => (e.xi, e.dim, e.deg, e.seq.Count)).Println($"Nb possibilities:{total}");

            var xseqs = seqs.Select(e => e.seq.MultiLoop(e.dim).Select(l => l.Aggregate((a0, a1) => a0 + a1))
                    .Distinct()
                    .Select(c => (e.xi, e.deg, e.dim, tr: c.Simplify()))
                    .ToArray())
                .ToArray();

            return xseqs.MultiLoop().Where(l => l.Aggregate(p, (p0, e) => p0.Substitute(e.tr, e.xi)).IsZero())
                .Select(l => l.Select(l0 => (l0.xi, l0.tr.Simplify())).ToArray());
        }

        throw new();
    }

    public static (Character<T> ext, Character<T> symm) ExtSymm(Character<T> chi)
    {
        var gr = chi.Gr;
        var classes = chi.Classes;
        var extMap = chi.Zero.Map.ToDictionary(kv => kv.Key, kv => kv.Value);
        var symmMap = chi.Zero.Map.ToDictionary(kv => kv.Key, kv => kv.Value);
        foreach (var g in classes)
        {
            var g2 = gr.Op(g, g);
            var cg = chi[g]!.Value;
            var cg2 = chi[g2]!.Value;
            extMap[g] = (cg * cg - cg2) / 2;
            symmMap[g] = (cg * cg + cg2) / 2;
        }

        return (new(classes, extMap), new(classes, symmMap));
    }

    public static void QuickDisplayTable(Character<T>[] list)
    {
        var classes = list[0].Classes;
        var ordCls = classes.Count.Range().Select(i => classes.GetRepresentative(i)).ToArray();
        var head = ordCls.Select(e => classes.GetClassName(e)).Prepend("Cl").ToArray();
        var rows = list.SelectMany((chi, i) => ordCls.Select(e => $"{chi[e]!.Value}").Prepend($"X.{i + 1}")).ToArray();
        Ring.DisplayMatrix(Ring.Matrix(list.Length + 1, head.Concat(rows).ToArray()), "  ");
        Console.WriteLine();
    }

    public static bool CheckCenterAndKernel(Character<T>[] list)
    {
        var chis = list.Order().ToArray();
        foreach (var (chi, i) in chis.Select((e, i) => (e, i)))
        {
            var g = chi.Gr;
            var z = Group.Generate($"Z[{i + 1}]", g, chi.Centre().ToArray());
            var ker = Group.Generate($"Ker[{i + 1}]", g, chi.Kernel().ToArray());
            Console.WriteLine($"X.{i + 1} = {chi}");
            Console.WriteLine($"Is Faithfull:{chi.IsFaithfull}");
            DisplayGroup.HeadOrders(z);
            DisplayGroup.HeadOrders(ker);
            if (!Group.IsNormalSubgroup(g, z) || !Group.IsNormalSubgroup(g, ker) || !Group.IsNormalSubgroup(z, ker))
                return false;

            var h = z.Over(ker);
            if (h.GroupType == GroupType.NonAbelianGroup)
                return false;
            
            var hType = Group.AbelianGroupType(h);
            if (hType.Length != 1)
                return false;

            var k = g.Over(ker);
            var kType = Group.AbelianGroupType(Group.Zentrum(k));
            if (kType.Length != 1 || kType[0] != hType[0])
                return false;
        }
        
        return true;
    }

    public static bool CheckExtSymmDecomposition(Character<T>[] list)
    {
        var chis = list.Order().ToArray();
        foreach (var (chi, i) in chis.Select((e, i) => (e, i)))
        {
            var (ext, symm) = ExtSymm(chi);
            var symmCheck = chis.Select((chi0, k) => (c: FG.InnerProduct(chi0, symm), Xi: k + 1))
                .Where(e => !e.c.IsZero() && !e.c.IsPositiveInteger)
                .ToArray();
            if (symmCheck.Length != 0)
            {
                if (Logger.Level != LogLevel.Off)
                {
                    QuickDisplayTable(chis);
                    Console.WriteLine($"X.{i + 1} = {chis[i]}");
                    Console.WriteLine($"X.{i + 1}_Symm = {symm}");
                    Console.WriteLine($"        = {symmCheck.Select(e => $"{e.c} * X.{e.Xi}").Glue(" + ")}");
                    Console.WriteLine();
                }
                
                return false;
            }

            var extCheck = chis.Select((chi0, k) => (c: FG.InnerProduct(chi0, ext), Xi: k + 1))
                .Where(e => !e.c.IsZero() && !e.c.IsPositiveInteger)
                .ToArray();
            if (extCheck.Length != 0)
            {
                if (Logger.Level != LogLevel.Off)
                {
                    QuickDisplayTable(chis);
                    Console.WriteLine($"X.{i + 1} = {chis[i]}");
                    Console.WriteLine($"X.{i + 1}_Ext = {ext}");
                    Console.WriteLine($"        = {extCheck.Select(e => $"{e.c} * X.{e.Xi}").Glue(" + ")}");
                    Console.WriteLine();
                }
                
                return false;
            }
        }

        return true;
    }
}