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
        Cells = new ACell[Classes.Length + 3, Classes.Length + 2];
        CnfCells = new CnfCell[Classes.Length, Classes.Length];
        IndexesTint = Classes.Select((e, i) => (e, i)).ToDictionary(e => e.e.repr, e => e.i);
        ReprOrbx = Classes.SelectMany(e => e.orbx.Select(ei => (e, ei))).ToDictionary(e => e.ei, e => e.e.repr);
        InitializeCells();
        if (Gr.GroupType == GroupType.AbelianGroup)
            AbelianCase();
        else
        {
            DerivedGroupLift();
        }
    }

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
        var Dg = Group.Derived(Gr);
        var Odg = Dg.Count();
        var GoDg = Gr.Over(Dg);
        var quo = (Quotient<T>)GoDg.BaseGroup;
        var n = GoDg.Count();
        if (n == 1)
            return;

        var ctGoDg = FG.CharactersTable(GoDg);

        var Ocl = Classes.ToDictionary(e => e.repr, e => e.stabx.Count);
        var derivedClasses = GoDg.ToDictionary(e => ReprOrbx[e.X], e => ctGoDg.IndexesTint[e]);
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

        var nbOcl = otherClasses.Count;
        if (nbOcl > 0 && nbOcl < 4)
        {
            var sum = IndexesTint.Where(e => !otherClasses.ContainsKey(e.Key))
                .Sum(e => (int)(CnfCells[e.Value, 0].E.Pow(2).E[0].Num));

            var x = Og - sum;
            var sol = IntExt.SolveSquareInt[nbOcl][x][0];
            foreach (var (e, i) in otherClasses.Select((e, i) => (e, i)))
                Cells[e.Value + 3, 2] = CnfCells[e.Value, 0] = new CnfCell(sol[i] * Cnf.CnfOne);
        }

        if (nbOcl > 0 && nbOcl < 4)
        {
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
            var ord = keys.Select(gi => (table[gi].T * table[ReprOrbx[Gr.Invert(gi)]])[0, 0] - Ocl[gi] * xz.One).ToArray();
            var eqs = orth.Concat(ord).Select(p => p.Num).Where(p => !p.IsZero()).ToArray();

            var sys = KMatrix<EPolynomial<Cnf>>.MergeSameRows(table.Values.ToArray());
            Console.WriteLine("#######################################################");
            DisplayGroup.Head(Gr);
            Ring.DisplayMatrix(sys.Coefs, "  ");

            Console.WriteLine(Gr.ShortName);
            eqs.Println("System");
            Ring.ReducedGrobnerBasis(eqs).Println("Solve");
        }
    }

    public int R => Classes.Length;
    private ACell[,] Cells { get; }
    private CnfCell[,] CnfCells { get; }
    public Dictionary<T, int> IndexesTint { get; }
    private Dictionary<T, T> ReprOrbx { get; }
    public (string name, T repr, HashSet<T> stabx, HashSet<T> orbx)[] Classes { get; }

    public void CheckProperties()
    {
        var n = Gr.Count();
        var mat = new KMatrix<Cnf>(Cnf.CnfZero, n, n);
        for (int i = 0; i < Classes.Length; i++)
        {
            for (int j = 0; j < Classes.Length; j++)
            {
                var cell = CnfCells[i, j];
                if (cell.IsEmpty)
                    return;

                mat.Coefs[i, j] = cell.E;
            }
        }

        var table = Classes.Select((e, i) => (i, e.repr)).ToDictionary(e => e.repr, e => mat.GetCol(e.i));
        var Ocl = Classes.ToDictionary(e => e.repr, e => e.stabx.Count);
        var e0 = table.Values.First().ToArray().First();
        var rg = Classes.Length.Range();
        var keys = table.Keys.ToArray();
        var allCombs = rg.SelectMany(i => rg.Where(j => j > i).Select(j => (keys[i], Gr.Invert(keys[j])))).ToArray();
        Console.WriteLine("Sum[Xi(g)Xj(g^−1 )]= 0  : {0}", allCombs.All(e => (table[e.Item1].T * table[ReprOrbx[e.Item2]]).IsZero()));
        Console.WriteLine("Sum[Xi(g)Xi(g^−1 )]=|G| : {0}",
            table.Keys.All(gi => (table[gi].T * table[ReprOrbx[Gr.Invert(gi)]])[0, 0].Equals(e0.One * Ocl[gi])));
    }

    public void DisplayCells()
    {
        Console.WriteLine(Gr.ShortName);
        Ring.DisplayMatrix(Cells, " ");
        CheckProperties();
        Console.WriteLine();
    }
}