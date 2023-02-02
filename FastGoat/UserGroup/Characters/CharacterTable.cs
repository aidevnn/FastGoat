using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.Characters;

public class CharacterTable<T> where T : struct, IElt<T>
{
    public ConcreteGroup<T> Gr { get; }

    public CharacterTable(ConcreteGroup<T> g)
    {
        Gr = g;
        Classes = Group.AllConjugacyClassesNames(Gr);
        CClasses = new ConjugacyClasses<T>(Gr);
        Cells = new ACell[Classes.Length + 3, Classes.Length + 2];
        CnfCells = new CnfCell[Classes.Length, Classes.Length];
        IndexesTint = Classes.Select((e, i) => (e, i)).ToDictionary(e => e.e.repr, e => e.i);
        Repr = Classes.Select((e, i) => (e, i)).ToDictionary(e => e.i, e => e.e.repr);
        ReprOrbx = Classes.SelectMany(e => e.orbx.Select(ei => (e, ei))).ToDictionary(e => e.ei, e => e.e.repr);
        NbOrbx = Classes.Select((e, i) => (e, i)).ToDictionary(e => e.i, e => e.e.orbx.Count);
        NbStabx = Classes.Select((e, i) => (e, i)).ToDictionary(e => e.i, e => e.e.stabx.Count);
        InitializeCells();
        if (Gr.GroupType == GroupType.AbelianGroup)
            AbelianCase();
        else
        {
            DerivedGroupLift();
            DerivedGroupCharacters();
            SolveSquareSum();
            SolveOrthogonality();
        }
    }

    public ConjugacyClasses<T> CClasses { get; }
    public int R => Classes.Length;
    private ACell[,] Cells { get; }
    public CnfCell[,] CnfCells { get; }
    public Dictionary<int, int> NbStabx { get; }
    public Dictionary<int, int> NbOrbx { get; }
    public Dictionary<T, int> IndexesTint { get; }
    public Dictionary<int, T> Repr { get; }
    private Dictionary<T, T> ReprOrbx { get; }
    public (string name, T repr, HashSet<T> stabx, HashSet<T> orbx)[] Classes { get; }
    private ConcreteGroup<T> DerivedGroup { get; set; }
    private ConcreteGroup<Coset<T>> QuotientGroup { get; set; }

    private void InitializeCells()
    {
        for (int i = 0; i < Classes.Length + 3; i++)
        {
            for (int j = 0; j < Classes.Length + 2; j++)
            {
                if (i == 0)
                    Cells[i, j] = j == 0 ? new Label("Class") : j == 1 ? new Label("   ") : new Label(Classes[j - 2].name);
                else if (i == 1)
                    Cells[i, j] = j == 0 ? new Label("Size") : j == 1 ? new Label("   ") : new Label($"{Classes[j - 2].orbx.Count}");
                else if (i == 2)
                    Cells[i, j] = new Label(" ");

                if (j == 0 && i > 2)
                    Cells[i, j] = new Label($"X.{i - 2}");

                if (j > 1 && i > 2)
                {
                    if (i == 3)
                        Cells[3, j] = CnfCells[0, j - 2] = new CnfCell(Cnf.CnfOne);
                    else
                        Cells[i, j] = CnfCells[i - 3, j - 2] = new CnfCell();
                }
            }
        }
    }

    private void AbelianCase()
    {
        var facts = AbelianInvariantsFactors.Reduce(Gr).ToArray();
        var gr0 = FG.Abelian(facts);
        var map = Group.AllMorphisms(Gr, gr0, Group.MorphismType.Isomorphism).First();
        var o = Gr.Count();
        var w = new Cnf(o);
        var wi = facts.Select(i => w.Pow(o / i)).ToArray();
        var nb = facts.Length.Range();

        Cnf Chi(int r, Ep<ZnInt> g)
        {
            var r0 = r;
            var prod = w.One;
            for (int i = 0; i < facts.Length; ++i)
            {
                var (q, ri) = Int32.DivRem(r0, facts[i]);
                prod *= wi[i].Pow(g.Ei[i].K * ri);
                r0 = q;
            }

            return prod;
        }

        for (int i = 0; i < o; i++)
        {
            for (int j = 0; j < o; j++)
            {
                var g = map[Classes[j].repr];
                Cells[i + 3, j + 2] = CnfCells[i, j] = new(Chi(i, g));
            }
        }
    }

    private void DerivedGroupLift()
    {
        var Og = Gr.Count();
        DerivedGroup = Group.Derived(Gr);
        var Odg = DerivedGroup.Count();
        if (Og == Odg)
            return;
        
        QuotientGroup = Gr.Over(DerivedGroup);
        var quo = (Quotient<T>)QuotientGroup.BaseGroup;
        var n = QuotientGroup.Count();

        var ctGoDg = FG.CharactersTable(QuotientGroup);

        var Ocl = Classes.ToDictionary(e => e.repr, e => e.stabx.Count);
        var derivedClasses = QuotientGroup.ToDictionary(e => ReprOrbx[e.X], e => ctGoDg.IndexesTint[e]);
        var otherClasses = IndexesTint.Where(e => !derivedClasses.ContainsKey(e.Key)).ToDictionary(e => e.Key, e => e.Value);
        foreach (var (gk, dk) in derivedClasses)
        {
            var k = IndexesTint[gk];
            foreach (var (gj, j) in IndexesTint)
            {
                var dgj = quo.GetRepresentative(gj);
                var dj = ctGoDg.IndexesTint[dgj];
                Cells[k + 3, j + 2] = CnfCells[k, j] = ctGoDg.CnfCells[dk, dj];
            }

            if (Ocl[gk] == Og / Odg)
            {
                foreach (var (_, i) in otherClasses)
                    Cells[i + 3, k + 2] = CnfCells[i, k] = new CnfCell(Cnf.CnfZero);
            }
        }
    }

    public void SolveSquareSum()
    {
        var otherClasses = IndexesTint.Where(e => IndexesTint.Any(f => CnfCells[e.Value, f.Value].IsEmpty))
            .ToDictionary(e => e.Key, e => e.Value);

        var Og = Gr.Count();
        var nbOcl = otherClasses.Count;
        DisplayCells();
        if (nbOcl > 0 && nbOcl < 4)
        {
            var sum = IndexesTint.Where(e => !otherClasses.ContainsKey(e.Key))
                .Sum(e => (int)(CnfCells[e.Value, 0].E.Pow(2).E[0].Num));

            var x = Og - sum;
            var sol = IntExt.SolveSquareInt[nbOcl][x][0];
            foreach (var (e, i) in otherClasses.Select((e, i) => (e, i)))
                Cells[e.Value + 3, 2] = CnfCells[e.Value, 0] = new CnfCell(sol[i] * Cnf.CnfOne);
        }
    }

    public void DerivedGroupCharacters()
    {
        if(DerivedGroup.Count() == Gr.Count())
            return;
        
        var clIdx = R.Range();
        var liftIdx = clIdx.Where(i => clIdx.All(j => !CnfCells[i, j].IsEmpty)).ToArray();

        var ctDg = FG.CharactersTable(DerivedGroup);
        var induced = new Dictionary<int, Dictionary<T, Cnf>>();
        foreach (var i in ctDg.IndexesTint.Values)
        {
            var indGH = IndexesTint.ToDictionary(e => e.Key, _ => Cnf.CnfZero);
            foreach (var g in CClasses.GetRepresentatives())
            {
                foreach (var y in Gr)
                {
                    var g0 = Gr.Op(Gr.Invert(y), Gr.Op(g, y));
                    if (!DerivedGroup.Contains(g0))
                        continue;

                    var j = ctDg.IndexesTint[ctDg.CClasses.GetRepresentative(g0)];
                    indGH[g] += ctDg.CnfCells[i, j].E;
                }

                indGH[g] /= DerivedGroup.Count();
            }

            induced[i] = indGH;
        }

        var restricted = IndexesTint.Values.ToDictionary(
            i => i,
            i => ctDg.IndexesTint.ToDictionary(e => e.Key, e => CnfCells[i, IndexesTint[CClasses.GetRepresentative(e.Key)]].E)
        );

        Cnf HProd(Dictionary<T, Cnf> chi1, Dictionary<T, Cnf> chi2)
        {
            var sum = Cnf.CnfZero;
            foreach (var g in Gr)
            {
                var rg = CClasses.GetRepresentative(g);
                sum += chi1[rg] * chi2[rg].Conj;
            }

            return sum / Gr.Count();
        }

        var chars = liftIdx.ToList();
        var deg = IntExt.Dividors(Gr.Count() / DerivedGroup.Count()).ToArray();
        foreach (var (a, ind) in induced)
        {
            var other = new List<int>(R.Range().Except(chars));
            if (other.Count == 0)
                break;

            if (chars.Select(r => R.Range().ToDictionary(j => Repr[j], j => CnfCells[r, j].E)).Any(r => !HProd(r, ind).IsZero()))
                continue;

            var e0 = ind[Gr.Neutral()].Simplify();
            if (e0.N != 1)
                throw new($"{e0.N} {e0} {e0.E.F}");

            var r0 = other.Min();
            var nb0 = Gr.Count() - other.Count * 4;
            var irrs = deg.Select(d => ind.ToDictionary(e => e.Key, e => e.Value / d))
                .Where(ind0 =>
                    HProd(ind0, ind0).Equals(Cnf.CnfOne) &&
                    ( ind0[Gr.Neutral()].Pow(2)).Module <= nb0)
                .ToArray();
            if (irrs.Length == 0)
                continue;

            chars.Add(r0);
            foreach (var (h, cnf) in irrs.Last())
            {
                var j = IndexesTint[h];
                Cells[r0 + 3, j + 2] = CnfCells[r0, j] = new(cnf);
            }
        }
    }

    public void SolveOrthogonality()
    {
        var clIdx = R.Range();
        var nbOcl = clIdx.Count(i => clIdx.Any(j => CnfCells[i, j].IsEmpty));
        if (nbOcl == 0)
            return;
        if (nbOcl > 2)
        {
            Console.WriteLine("Missing 3 characters or more");
            return;
        }

        var cells = new List<(int, int)>();
        for (int i = 0; i < Classes.Length; i++)
        {
            for (int j = 0; j < Classes.Length; j++)
            {
                if (CnfCells[i, j].IsEmpty)
                    cells.Add((i, j));
            }
        }

        var xis = Ring.EPolynomial(Cnf.CnfZero, MonomOrder.Lex, (cells.Count + 1, "x"));
        var mapCells = xis.SkipLast(1).Select((e, i) => (e, i)).ToDictionary(e => e.e, e => cells[e.i]);
        var mapSymb = mapCells.ToDictionary(e => e.Value, e => e.Key);
        var mapInd = mapCells.ToDictionary(e => e.Key.Num.ExtractIndeterminate, e => e.Value);
        var xz = xis.Last();
        var table = new Dictionary<T, KMatrix<EPolynomial<Cnf>>>();

        foreach (var (gi, i) in IndexesTint)
        {
            var mat = new KMatrix<EPolynomial<Cnf>>(xz, Classes.Length, 1);
            for (int j = 0; j < Classes.Length; j++)
            {
                var c = CnfCells[j, i];
                if (c.IsEmpty)
                    mat.Coefs[j, 0] = mapSymb[(j, i)];
                else
                    mat.Coefs[j, 0] = c.E * xz.One;
            }

            table[gi] = mat;
        }

        var rg = table.Count.Range();
        var keys = table.Keys.ToArray();
        var allCombs = rg.SelectMany(i => rg.Where(j => j > i).Select(j => (keys[i], Gr.Invert(keys[j])))).ToArray();
        var orth = allCombs.Select(e => (table[e.Item1].T * table[ReprOrbx[e.Item2]])[0, 0]).ToArray();
        var ord = keys.Select(gi => (table[gi].T * table[ReprOrbx[Gr.Invert(gi)]])[0, 0] - NbStabx[IndexesTint[gi]] * xz.One)
            .ToArray();
        var eqs = orth.Concat(ord).Select(p => p.Num).Where(p => !p.IsZero()).ToArray();
        var redEqs = Ring.ReducedGrobnerBasis(eqs);

        var sys = KMatrix<EPolynomial<Cnf>>.MergeSameRows(table.Values.ToArray());
        DisplayGroup.Head(Gr);
        // Ring.DisplayMatrix(sys.Coefs, "  ");

        eqs.Println("System");
        redEqs.Println("Reduced System");
        Console.WriteLine();

        var mapCnf = mapSymb.ToDictionary(e => e.Value.Num.ExtractIndeterminate, e => Cnf.Nth(NbStabx[e.Key.Item2]));
        var allSolutions = SolveSystem(new Dictionary<Xi, Cnf>(), redEqs, mapCnf, xz.Num).ToArray();
        foreach (var solution in allSolutions)
        {
            solution.Println("Solution");
            Console.WriteLine();
            break;
        }

        var firstSol = allSolutions[0];
        foreach (var (xi, cnf) in firstSol)
        {
            var (i, j) = mapInd[xi];
            Cells[i + 3, j + 2] = CnfCells[i, j] = new CnfCell(cnf);
        }
    }

    private IEnumerable<Dictionary<Xi, Cnf>> SolveSystem(Dictionary<Xi, Cnf> solutions, Polynomial<Cnf, Xi>[] eqs,
        Dictionary<Xi, Cnf> mapCnf, Polynomial<Cnf, Xi> xz)
    {
        var oneIndeterminate = eqs.Where(eq => eq.NbIndeterminates == 1).OrderBy(eq => eq.Coefs.Keys.Max(m => m.Degree)).ToArray();
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
            // pos0.Println("Possibility");
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
                // Console.WriteLine();
                foreach (var nsol0 in SolveSystem(nsol, subsEq, mapCnf, xz))
                    yield return nsol0;
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
        var facts0 = PolynomialFactorizationPart2.FirrZ(f);
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
            var facts = AlgebraicFactorization.AlgebraicFactors(P);
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

    public void CheckProperties()
    {
        var n = Classes.Length;
        var mat = new KMatrix<Cnf>(Cnf.CnfZero, n, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                var cell = CnfCells[i, j];
                if (cell.IsEmpty)
                    return;

                mat.Coefs[i, j] = cell.E;
            }
        }

        var rg = Classes.Length.Range();
        var Ocl = Classes.ToDictionary(e => e.repr, e => e.stabx.Count);
        var e0 = Cnf.CnfZero;
        var allCombs = rg.SelectMany(i => rg.Where(j => j > i).Select(j => (i, j))).ToArray();
        var ggi = Gr.Select(e => (g: CClasses.GetIndex(e), gi: CClasses.GetIndex(Gr.Invert(e)))).ToArray();
        var clggi = rg.ToDictionary(i => i, i => IndexesTint[CClasses.GetRepresentative(Gr.Invert(Repr[i]))]);

        Console.WriteLine("All i,                 Sum[g](Xi(g)Xi(g^−1))= |G|      : {0}",
            rg.All(i => ggi.Aggregate(e0.Zero, (sum, kp) => sum + mat[i, kp.g] * mat[i, kp.gi]).Equals(e0.One * Gr.Count())));
        Console.WriteLine("All i <> j,            Sum[g](Xi(g)Xj(g^−1))=  0       : {0}",
            allCombs.All(e => ggi.Aggregate(e0.Zero, (sum, kp) => sum + mat[e.i, kp.g] * mat[e.j, kp.gi]).IsZero()));

        Console.WriteLine("All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1))= |Cl(g)|  : {0}",
            clggi.All(kp =>
                rg.Aggregate(e0.Zero, (sum, r) => sum + mat[r, kp.Key] * mat[r, kp.Value]).Equals(e0.One * Ocl[Repr[kp.Key]])));
        Console.WriteLine("All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1))=  0       : {0}",
            allCombs.All(e => rg.Aggregate(e0.Zero, (sum, r) => sum + mat[r, e.i] * mat[r, clggi[e.j]]).IsZero()));
    }

    public void DisplayCells()
    {
        Console.WriteLine(Gr.ShortName);
        Ring.DisplayMatrix(Cells, " ");
        CheckProperties();
        Console.WriteLine();
    }
}