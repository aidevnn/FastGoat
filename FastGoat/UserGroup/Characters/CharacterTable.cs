using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.Characters;

public class CharacterTable<T> where T : struct, IElt<T>
{
    public CharacterTable(ConcreteGroup<T> gr, bool empty = false)
    {
        Gr = gr;
        var one = new Character<T>(gr);
        Classes = one.Classes;
        NbClasses = Classes.Count();
        AllCharacters = NbClasses.Range().Select(i => one.Zero).ToArray();
        ChiE = AllCharacters[0] = one;
        if (!empty)
        {
            if (gr.GroupType == GroupType.AbelianGroup)
            {
                AbelianTable();
            }
            else
            {
                DerivedSubGroupLift();
                InductionFromSubGroups();
            }
        }
    }

    public Character<T> ChiE { get; }

    public ConcreteGroup<T> Gr { get; }
    public ConjugacyClasses<T> Classes { get; }
    public Character<T>[] AllCharacters { get; set; }
    public int NbClasses { get; }

    public void ClearTable()
    {
        AllCharacters = NbClasses.Range().Select(i => ChiE.Zero).ToArray();
        AllCharacters[0] = ChiE;
    }

    public void AbelianTable()
    {
        if (Gr.GroupType != GroupType.AbelianGroup)
            throw new GroupException(GroupExceptionType.GroupDef);

        AllCharacters = FG.LinearCharacters(Gr).ToArray();
    }

    public void DerivedSubGroupLift()
    {
        var dg = Group.Derived(Gr);
        if (dg.Count() < Gr.Count())
        {
            var quo = Gr.Over(dg);
            var ctQuo = new CharacterTable<Coset<T>>(quo);
            foreach (var chi in ctQuo.AllCharacters.Where(chi => chi.HasAllValues))
            {
                var lift = FG.Lift(chi, Classes);
                var liftState = AddCharacter(lift);
                if (liftState == AddCharacterState.TableFull)
                    return;
            }

            InductionFromSubGroup(dg);
        }
    }

    public void InductionFromSubGroups()
    {
        if (AllCharacters.All(chi => chi.HasAllValues))
            return;

        var og = Gr.Count();
        var all = Group.AllSubGroups(Gr);
        var allSubs = all.Keys.Where(g => g.Count() != 1 && g.Count() != og).ToArray();
        foreach (var sg in allSubs.OrderBy(sg => sg.Count()))
        {
            InductionFromSubGroup(sg);
            if (AllCharacters.All(chi => chi.HasAllValues))
                return;
        }

        var allSubs2 = all.Values.SelectMany(e => e).Where(g => g.Count() != 1 && g.Count() != og).ToArray();
        foreach (var sg in allSubs2.OrderBy(sg => sg.Count()))
        {
            InductionFromSubGroup(sg);
            if (AllCharacters.All(chi => chi.HasAllValues))
                return;
        }
    }

    public void SolveSumSquare()
    {
        var ne = Gr.Neutral();
        var doneChis = AllCharacters.Where(chi => chi[ne].HasValue).ToArray();
        var oG = Gr.Count();
        var sum = (int)doneChis.Sum(chi => Double.Pow(chi[ne]!.Value.Module, 2));
        var rem = Gr.Count() - sum;
        var nbChis = NbClasses - doneChis.Length;
        if (nbChis == 0 || !IntExt.SolveSquareInt.ContainsKey(nbChis) || !IntExt.SolveSquareInt[nbChis].ContainsKey(rem))
            return;

        var sq = IntExt.SolveSquareInt[nbChis][rem];
        var cd = sq.Where(l => l.All(i => oG % i == 0)).ToArray();
        if (cd.Length == 1)
        {
            var lt = cd[0];
            var todoChis = AllCharacters.Where(e => !e.HasAllValues).ToArray();
            for (int i = 0; i < nbChis; i++)
            {
                todoChis[i].Map[ne] = lt[i] * Cnf.CnfOne;
            }
        }
    }

    public void InductionFromSubGroup(ConcreteGroup<T> subGr)
    {
        if (!subGr.SubSetOf(Gr))
            return;

        SolveSumSquare();
        AllCharacters = AllCharacters.Order().ToArray();
        if (AllCharacters.All(chi => chi.HasAllValues))
            return;

        var ctSubGr = new CharacterTable<T>(subGr);
        foreach (var chi in ctSubGr.AllCharacters.Where(chi => chi.HasAllValues))
        {
            var ind = FG.Induction(chi, Classes);
            var indState = AddCharacter(ind);
            if (indState == AddCharacterState.TableFull)
                return;
        }
    }

    public void InductionFromSubGroup(CharacterTable<T> ctSubGr)
    {
        if (!ctSubGr.Gr.SubSetOf(Gr))
            return;

        SolveSumSquare();
        AllCharacters = AllCharacters.Order().ToArray();
        if (AllCharacters.All(chi => chi.HasAllValues))
            return;

        foreach (var chi in ctSubGr.AllCharacters.Where(chi => chi.HasAllValues))
        {
            var ind = FG.Induction(chi, Classes);
            var indState = AddCharacter(ind);
            if (indState == AddCharacterState.TableFull)
                return;
        }
    }

    public void InductionFromSubGroup<T2>(ConcreteGroup<T2> sg) where T2 : struct, IElt<T2>
    {
        var sg0 = Group.IsomorphicSubgroup(Gr, sg);
        foreach (var sg1 in Group.SubGroupsConjugates(Gr, sg0))
            InductionFromSubGroup(sg1);
    }

    public void RestrictionFromSuperGroup(ConcreteGroup<T> superGr)
    {
        if (!Gr.SubSetOf(superGr))
            return;

        var ne = Gr.Neutral();
        AllCharacters = AllCharacters.Order().ToArray();
        if (AllCharacters.All(chi => chi.HasAllValues))
            return;

        var ctSuperGr = new CharacterTable<T>(superGr);
        foreach (var chi in ctSuperGr.AllCharacters.Where(chi => chi.HasAllValues))
        {
            var res = FG.Restriction(chi, Classes);
            var resState = AddCharacter(res);
            if (resState == AddCharacterState.TableFull)
                return;
        }
    }

    public void RestrictionFromSuperGroup()
    {
        var superGr = Gr.SuperGroup;
        if (superGr is not null && superGr.Count() > Gr.Count())
            RestrictionFromSuperGroup(superGr);
    }

    public void RestrictionFromSuperGroup(CharacterTable<T> ctSuperGr)
    {
        if (!Gr.SubSetOf(ctSuperGr.Gr))
            return;

        var ne = Gr.Neutral();
        AllCharacters = AllCharacters.Order().ToArray();
        if (AllCharacters.All(chi => chi.HasAllValues))
            return;

        foreach (var chi in ctSuperGr.AllCharacters.Where(chi => chi.HasAllValues))
        {
            var res = FG.Restriction(chi, Classes);
            var resState = AddCharacter(res);
            if (resState == AddCharacterState.TableFull)
                return;
        }
    }

    public AddCharacterState AddCharacter(Character<T> chi1)
    {
        var (chiState, chi2) = AddCharacterInternal(chi1);
        if (chiState == AddCharacterState.Done)
        {
            var e = Gr.Neutral();
            if (chi2[e]!.Value.Module > 1)
            {
                var linears = AllCharacters.Where(chi => chi.HasAllValues && !chi.Equals(ChiE) && chi[e]!.Value.Equals(Cnf.CnfOne))
                    .ToArray();
                foreach (var chiLinear in linears)
                {
                    var chi3 = chi2 * chiLinear;
                    var (chiLinState, _) = AddCharacterInternal(chi3);
                    // Console.WriteLine(new { chiLinState, chiLinear, chi2, chi3 });
                    if (chiLinState == AddCharacterState.TableFull)
                        return AddCharacterState.TableFull;
                }
            }
        }

        return chiState;
    }

    private (AddCharacterState, Character<T>) AddCharacterInternal(Character<T> chi1)
    {
        var e = Gr.Neutral();
        var doneChis = AllCharacters.Where(c => c.HasAllValues).ToList();
        var todoChis = AllCharacters.Select((c, i) => (c, i)).Where(c => !c.c.HasAllValues).ToList();

        if (todoChis.Count == 0)
            return (AddCharacterState.TableFull, chi1);

        var coefs = doneChis.Select(c => (c, FG.InnerProduct(chi1, c))).ToArray();
        var chi2 = coefs.Aggregate(chi1, (id, e0) => id - e0.Item2 * e0.c).Simplify();
        if (chi2.IsZero())
            return (AddCharacterState.Rejected, chi2);

        var r2 = chi2[e]!.Value.Module;
        if (!Double.IsInteger(r2) || (int)r2 < 1)
            return (AddCharacterState.Rejected, chi2);

        var nbDone = doneChis.Count;
        var shortName = Gr.ShortName;
        // Console.WriteLine(new { shortName, NbClasses, nbDone, chi2 });

        if (doneChis.All(chi => FG.InnerProduct(chi2, chi).IsZero()))
        {
            // if (!chi1.Equals(chi2))
            //     Console.WriteLine($"Orthogonal   Xi = {chi1} => {chi2}");

            var prod = FG.InnerProduct(chi2, chi2).Simplify();
            if (prod.Equals(Cnf.CnfOne))
            {
                var dim = chi2[e]!.Value.Module;
                var idx = todoChis.FindIndex(c => c.c[e].HasValue && c.c[e]!.Value.Module.Equals(dim));
                var name = Gr.ShortName;
                // Console.WriteLine("New Character {0}", new { idx, dim, name, chi2 });
                var (_, i) = idx == -1 ? todoChis.First() : todoChis[idx];
                AllCharacters[i] = chi2.Simplify();
                return (todoChis.Count > 1 ? AddCharacterState.Done : AddCharacterState.TableFull, chi2);
            }

            // Console.WriteLine($"{Gr.ShortName} Orthogonal Not Irreductible Xi = {chi1} => {chi2}");
            return (AddCharacterState.NotIrr, chi2);
        }

        return (AddCharacterState.NotOrth, chi2);
    }

    public void CheckProperties()
    {
        if (AllCharacters.Any(chi => !chi.HasAllValues))
            return;

        var rg = NbClasses.Range();
        var Ocl = Classes.ToDictionary(e => e, e => Classes.GetClassStabx(e).Count());
        var e0 = Cnf.CnfZero;
        var allCombs = rg.SelectMany(i => rg.Where(j => j > i)
            .Select(j => (gi: Classes.GetRepresentative(i), gj: Classes.GetRepresentative(j), i, j))).ToArray();
        var ggi = Gr.Select(e => (g: e, gi: Gr.Invert(e))).ToArray();
        var clggi = Classes.ToDictionary(e => e, e => Classes.GetRepresentative(Gr.Invert(e)));

        var chis = AllCharacters;
        Console.WriteLine("All i,                 Sum[g](Xi(g)Xi(g^−1)) = |G|      : {0}",
            rg.All(i => ggi.Aggregate(e0.Zero, (sum, kp) => sum + chis[i][kp.g]!.Value * chis[i][kp.gi]!.Value)
                .Equals(e0.One * Gr.Count())));
        Console.WriteLine("All i <> j,            Sum[g](Xi(g)Xj(g^−1)) =  0       : {0}",
            allCombs.All(e => ggi.Aggregate(e0.Zero, (sum, kp) => sum + chis[e.i][kp.g]!.Value * chis[e.j][kp.gi]!.Value).IsZero()));

        Console.WriteLine("All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : {0}",
            clggi.All(kp =>
                rg.Aggregate(e0.Zero, (sum, r) => sum + chis[r][kp.Key]!.Value * chis[r][kp.Value]!.Value)
                    .Equals(e0.One * Ocl[kp.Key])));
        Console.WriteLine("All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : {0}",
            allCombs.All(e =>
                rg.Aggregate(e0.Zero, (sum, r) => sum + chis[r][e.gi]!.Value * chis[r][clggi[e.gj]]!.Value).IsZero()));
    }

    public static (string c, string minus_c, string conj_c, string minus_conj) PrettyPrintCnf(Cnf c)
    {
        var cStr = $"{c}";
        if (!cStr.Contains(' ') || cStr.Contains('I'))
            return (cStr, $"{-c}", $"{c.Conj}", $"{-c.Conj}");

        if (c.Im.IsZero())
        {
            var p0 = IntFactorisation.PrettyPrintCnf(c);
            if (!$"{p0}".Contains('ξ'))
            {
                var p1 = -p0;
                var p2 = p0;
                var p3 = -p2;
                return ($"{p0}", $"{p1}", $"{p2}", $"{p3}");
            }
        }
        else
        {
            var p0 = IntFactorisation.PrettyPrintCnf(c);
            if (!$"{p0}".Contains('ξ'))
            {
                var p1 = -p0;
                var p2 = IntFactorisation.PrettyPrintCnf(c.Conj);
                var p3 = -p2;
                return ($"{p0}", $"{p1}", $"{p2}", $"{p3}");
            }
        }

        return (cStr, $"{-c}", $"{c.Conj}", $"{-c.Conj}");
    }

    private string[,] PrepareTable()
    {
        var idxs = NbClasses.Range().Grid2D().ToArray();
        var allCnfs = idxs.Select(e =>AllCharacters[e.t1][Classes.GetRepresentative(e.t2)])
            .Where(e => e.HasValue).Select(e =>  e!.Value.Simplify()).ToHashSet();
        var cnf2str = new Dictionary<Cnf, string>(NbClasses.Pow(2));

        foreach (var c in allCnfs)
        {
            if (!cnf2str.ContainsKey(c))
            {
                var (c0, c1, c2, c3) = PrettyPrintCnf(c);
                cnf2str[c] = c0;
                if (!c.IsZero())
                {
                    cnf2str[-c] = c1;
                    if (!c.Im.IsZero())
                    {
                        cnf2str[c.Conj] = c2;
                        cnf2str[-c.Conj] = c3;
                    }
                }
            }
        }

        var cellStrs = new string[NbClasses, NbClasses];
        foreach (var (i, j) in idxs)
        {
            var c = AllCharacters[i][Classes.GetRepresentative(j)];
            cellStrs[i, j] = !c.HasValue ? "#" : " " + cnf2str[c.Value];
        }

        return cellStrs;
    }

    public void DisplayCells()
    {
        AllCharacters = AllCharacters.Order().ToArray();
        var form = Ring.MatrixDisplayForm;
        Ring.MatrixDisplayForm = Ring.MatrixDisplay.Table;
        var Cells = new string[NbClasses + 3, NbClasses + 2];
        var StrTable = PrepareTable();
        for (int i = 0; i < NbClasses + 3; i++)
        {
            for (int j = 0; j < NbClasses + 2; j++)
            {
                if (i == 0)
                    Cells[i, j] = j == 0 ? "Class" : j == 1 ? "   " : Classes.GetClassName(j - 2);
                else if (i == 1)
                    Cells[i, j] = j == 0 ? "Size" : j == 1 ? "   " : $"{Classes.GetClassOrbx(j - 2).Count()}";
                else if (i == 2)
                    Cells[i, j] = " ";

                if (j == 0 && i > 2)
                    Cells[i, j] = $"Ꭓ.{i - 2}";

                if (j > 1 && i > 2)
                    Cells[i, j] = StrTable[i - 3, j - 2];
            }
        }

        DisplayGroup.Head(Gr);
        Ring.DisplayMatrix(Cells, " ");
        CheckProperties();
        Console.WriteLine();

        Ring.MatrixDisplayForm = form;
    }
}