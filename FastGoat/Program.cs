using System.Collections;
using System.ComponentModel;
using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

(List<(Xi, Polynomial<EPoly<Rational>, Xi>)> allSols, Polynomial<EPoly<Rational>, Xi>[] sys0)
    Substitutions(Polynomial<EPoly<Rational>, Xi>[] sys)
{
    var sys0 = sys.ToArray();
    var k = 1;
    var allSols = new List<(Xi, Polynomial<EPoly<Rational>, Xi>)>();

    while (true)
    {
        var sols = new List<(Xi, Polynomial<EPoly<Rational>, Xi>)>();
        foreach (var eq in sys0.Where(eq => eq.Coefs.Any(m => m.Key.Degree == 1))
                     .OrderBy(eq => eq.NbIndeterminates)
                     .ThenByDescending(eq => eq.ExtractAllIndeterminates.Max(x0 => x0.xi)).ToArray())
        {
            var xi0 = eq.ExtractAllIndeterminates.Where(m => eq.Coefs.Keys.Count(m1 => m1[m] == 1 && m1.Degree == 1) == 1).ToArray();
            if (xi0.Length == 0)
                continue;

            var xi = xi0[0];
            var xeq = eq.X(xi).ExtractMonom;
            var a = eq[xeq];
            var b = eq - a * eq.X(xi);
            var c = -b / a;
            if (c.ExtractAllIndeterminates.Contains(xi))
                continue;
            
            sols.Add((xi, c));
            allSols.Add((xi, c));
        }

        sols.Println($"Subs{k++}");
        var sys1 = sys0.Select(eq => eq.Substitute(sols)).Where(eq => !eq.IsZero()).ToArray();
        sys1.Println("Sys1");
        sys0 = Ring.ReducedGrobnerBasis(sys1);
        if (sols.Count == 0)
            break;
    }

    allSols.OrderBy(e => e.Item1).Select(e => $"{e.Item1} = {e.Item2}").Println();
    return (allSols, sys0);
}

void Orthogonality<T>(CharacterTable<T> ct) where T : struct, IElt<T>
{
    var cells = new List<(int, int)>();
    var cls = ct.Classes.OrderBy(e => ct.Classes.GetClassName(e)).ToArray();
    var lcm = Lcm(cls.Select(g => ct.Gr.ElementsOrders[g]).ToArray());
    var e0 = FG.CyclotomicEPoly(lcm).X;
    cls.Select(e => ct.Classes.GetClassName(e)).Println($"lcm:{lcm} e0:{e0} F:{e0.F}");
    for (int i = 0; i < ct.Classes.Count; i++)
    {
        var chi = ct.AllCharacters[i];
        for (int j = 0; j < ct.Classes.Count; j++)
        {
            var e = cls[j];
            if (!chi[e].HasValue)
                cells.Add((i, j));
        }
    }

    var xis = Ring.EPolynomial(e0, MonomOrder.Lex, (cells.Count + 1, "x"));
    var mapCells = xis.SkipLast(1).Select((e, i) => (e, i)).ToDictionary(e => e.e, e => cells[e.i]);
    var mapSymb = mapCells.ToDictionary(e => e.Value, e => e.Key);
    var mapInd = mapCells.ToDictionary(e => e.Key.Num.ExtractIndeterminate, e => e.Value);
    var xz = xis.Last();
    var table = new Dictionary<T, KMatrix<EPolynomial<EPoly<Rational>>>>();

    foreach (var (gi, i) in cls.Select((gi, i) => (gi, i)))
    {
        var mat = new KMatrix<EPolynomial<EPoly<Rational>>>(xz, ct.Classes.Count, 1);
        for (int j = 0; j < ct.Classes.Count; j++)
        {
            var c = ct.AllCharacters[j][gi];
            if (!c.HasValue)
                mat.Coefs[j, 0] = mapSymb[(j, i)];
            else
            {
                var ek = e0.Pow(lcm / c.Value.N);
                var c0 = c.Value.E.Substitute(ek);
                mat.Coefs[j, 0] = c0 * xz.One;
            }
        }

        table[gi] = mat;
    }

    var rg = table.Count.Range();
    var keys = table.Keys.ToArray();
    var NbStabx = ct.Classes.GetStabxCount().ToArray();
    var allCombs = rg.SelectMany(i => rg.Where(j => j > i).Select(j => (keys[i], ct.Gr.Invert(keys[j])))).ToArray();
    var orth = allCombs.Select(e => (table[e.Item1].T * table[ct.Classes.GetRepresentative(e.Item2)])[0, 0]).ToArray();
    var ord = keys.Select(gi =>
            (table[gi].T * table[ct.Classes.GetRepresentative(ct.Gr.Invert(gi))])[0, 0] - NbStabx[ct.Classes.GetIndex(gi)] * xz.One)
        .ToArray();
    var eqs = orth.Concat(ord).Select(p => p.Num).Where(p => !p.IsZero()).ToArray();
    Console.WriteLine(ct.Gr.ShortName);
    eqs.Println("System");

    var sys = KMatrix<EPolynomial<EPoly<Rational>>>.MergeSameRows(table.Values.ToArray());
    DisplayGroup.Head(ct.Gr);
    Ring.DisplayMatrix(sys.Coefs, "  ");

    GlobalStopWatch.Restart();
    var redEqs = Ring.ReducedGrobnerBasis(eqs);
    redEqs.Println("Reduced System");
    var (sols, sys2) = Substitutions(redEqs);
    sys2.Println("Sys2");
    GlobalStopWatch.Show();
    Console.WriteLine();
    foreach (var (xi, c) in sols.Where(e => e.Item2.Degree == 0))
    {
        var cnf = new Cnf(lcm, c.ConstTerm);
        var (i, j) = mapInd[xi];
        var chi = ct.AllCharacters[i];
        chi.Map[cls[j]] = cnf;
        ct.AllCharacters[i] = new(ct.Classes, chi.Map);
    }
}

{
    var g = FG.SL2p(3);
    var ct = FG.CharacterTableEmpty(g);
    ct.DerivedSubGroupLift();
    ct.DisplayCells();
    Orthogonality(ct);
    ct.DisplayCells();
}