using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Subgroups;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.Characters;

public partial class CharacterTable<T> where T : struct, IElt<T>
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
                InductionFromStabilizers();
                InductionFromSubGroups();
            }
        }
    }

    public Character<T> ChiE { get; }

    public ConcreteGroup<T> Gr { get; }
    public ConjugacyClasses<T> Classes { get; }
    public Character<T>[] AllCharacters { get; set; }
    public int NbClasses { get; }
    public bool TableComplete => AllCharacters.All(chi => chi.HasAllValues);
    public Character<T>[] DoneChis => AllCharacters.Where(chi => chi.HasAllValues).ToArray();
    public Character<T>[] ToDoChis => AllCharacters.Where(chi => !chi.HasAllValues).ToArray();

    public void ClearTable()
    {
        AllCharacters = NbClasses.Range().Select(i => ChiE.Zero).ToArray();
        AllCharacters[0] = ChiE;
    }

    public void AbelianTable()
    {
        if (Gr.GroupType != GroupType.AbelianGroup)
            throw new GroupException(GroupExceptionType.GroupDef);

        AllCharacters = FG.LinearCharacters(Gr).Order().ToArray();
    }

    public void DerivedSubGroupLift()
    {
        if (TableComplete)
            return;

        var der = Group.Derived(Gr);
        var od = der.Count();
        if (od == Gr.Count())
            return;

        NormalSubGroupLift(der);
        InductionFromSubGroup(der);
    }

    public void NormalSubGroupLift(ConcreteGroup<T> normal)
    {
        var quo = Gr.Over(normal);
        var ctQuo = new CharacterTable<Coset<T>>(quo);
        foreach (var chi in ctQuo.DoneChis)
        {
            var lift = FG.Lift(chi, Classes);
            var liftState = AddCharacter(lift);
            if (liftState == AddCharacterState.TableFull)
                return;
        }
    }

    public void InductionFromStabilizers()
    {
        if (TableComplete)
            return;

        var stabilizers = Classes.Skip(1)
            .Select(e => Group.Generate($"Stab({Classes.GetClassName(e)})", Gr, Classes.GetClassStabx(e).ToArray()))
            .Where(sg => sg.Count() != Gr.Count());
        foreach (var sg in stabilizers)
        {
            if (TableComplete)
                break;

            if (Gr.All(x => sg.Select(y => Gr.Op(x, y)).ToHashSet().SetEquals(sg.Select(y => Gr.Op(y, x)))))
                NormalSubGroupLift(sg);

            InductionFromSubGroup(sg);
        }
    }

    public void InductionFromSubGroups()
    {
        if (TableComplete)
            return;

        InductionFromSubGroups(Gr.AllSubgroups());
    }

    public void InductionFromSubGroups(AllSubgroups<T> subgroups)
    {
        if (TableComplete)
            return;

        foreach (var cj in subgroups.AllSubgroupConjugates.Where(cj => !cj.IsTrivial && cj.IsProper))
        {
            if (cj.IsNormal)
                NormalSubGroupLift(cj.Representative);

            InductionFromSubGroup(cj.Representative);
            if (TableComplete)
                return;
        }

        foreach (var sg in subgroups.AllSubgroupConjugates.Where(cj => !cj.IsNormal && cj.IsProper)
                     .SelectMany(cj => cj.Conjugates.Skip(1)))
        {
            InductionFromSubGroup(sg);
            if (TableComplete)
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
        if (nbChis == 0 || !IntExt.SolveSquareInt.ContainsKey(nbChis) ||
            !IntExt.SolveSquareInt[nbChis].ContainsKey(rem))
            return;

        var sq = IntExt.SolveSquareInt[nbChis][rem];
        var cd = sq.Where(l => l.All(i => oG % i == 0)).ToArray();
        if (cd.Length == 1)
        {
            var lt = cd[0];
            var todoChis = ToDoChis;
            for (int i = 0; i < nbChis; i++)
            {
                todoChis[i].Map[ne] = lt[i] * Cnf.CnfOne;
            }
        }
    }

    private bool RegularCharacter()
    {
        if (AllCharacters.Count(chi => !chi.HasAllValues) != 1)
            return false;

        var e = Gr.Neutral();
        var doneChis = DoneChis;
        var idx = AllCharacters.ToList().FindIndex(chi => !chi.HasAllValues);

        Cnf? Starting(T e0) => e0.Equals(e) ? Cnf.CnfOne * Gr.Count() : Cnf.CnfZero;
        var map = Classes.ToDictionary(
            cl => cl,
            cl => doneChis.Aggregate(Starting(cl), (sum, chi) => sum - chi[e] * chi[cl])
        );

        var c0 = AllCharacters[idx][e];
        foreach (var (k, v) in map.ToArray())
            map[k] = v / c0;

        AllCharacters[idx] = new(Classes, map);
        return true;
    }

    public void InductionFromSubGroup(ConcreteGroup<T> subGr)
    {
        if (!subGr.SubSetOf(Gr))
            return;

        SolveSumSquare();
        AllCharacters = AllCharacters.Order().ToArray();
        if (TableComplete)
            return;

        var ctSubGr = new CharacterTable<T>(subGr);
        foreach (var chi in ctSubGr.DoneChis)
        {
            var ind = FG.Induction(chi, Classes);
            var indState = AddCharacter(ind);
            if (indState == AddCharacterState.TableFull || RegularCharacter())
                return;
        }
    }

    public void InductionFromSubGroup(CharacterTable<T> ctSubGr)
    {
        if (!ctSubGr.Gr.SubSetOf(Gr))
            return;

        SolveSumSquare();
        AllCharacters = AllCharacters.Order().ToArray();
        if (TableComplete)
            return;

        foreach (var chi in ctSubGr.DoneChis)
        {
            var ind = FG.Induction(chi, Classes);
            var indState = AddCharacter(ind);
            if (indState == AddCharacterState.TableFull || RegularCharacter())
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

        AllCharacters = AllCharacters.Order().ToArray();
        if (TableComplete)
            return;

        var ctSuperGr = new CharacterTable<T>(superGr);
        foreach (var chi in ctSuperGr.DoneChis)
        {
            var res = FG.Restriction(chi, Classes);
            var resState = AddCharacter(res);
            if (resState == AddCharacterState.TableFull || RegularCharacter())
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

        AllCharacters = AllCharacters.Order().ToArray();
        if (TableComplete)
            return;

        foreach (var chi in ctSuperGr.DoneChis)
        {
            var res = FG.Restriction(chi, Classes);
            var resState = AddCharacter(res);
            if (resState == AddCharacterState.TableFull || RegularCharacter())
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
                var linears = AllCharacters.Where(chi => chi.IsLinear && !chi.Equals(ChiE)).ToArray();
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
        var doneChis = DoneChis;
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

        if (doneChis.All(chi => FG.InnerProduct(chi2, chi).IsZero()))
        {
            var prod = FG.InnerProduct(chi2, chi2).Simplify();
            if (prod.Equals(Cnf.CnfOne))
            {
                var dim = chi2[e]!.Value.Module;
                var idx = todoChis.FindIndex(c => c.c[e].HasValue && c.c[e]!.Value.Module.Equals(dim));
                var (_, i) = idx == -1 ? todoChis.First() : todoChis[idx];
                AllCharacters[i] = chi2;
                return (todoChis.Count > 1 ? AddCharacterState.Done : AddCharacterState.TableFull, chi2);
            }

            return (AddCharacterState.NotIrr, chi2);
        }

        return (AddCharacterState.NotOrth, chi2);
    }

    public void CheckProperties()
    {
        if (!TableComplete)
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
            allCombs.All(e =>
                ggi.Aggregate(e0.Zero, (sum, kp) => sum + chis[e.i][kp.g]!.Value * chis[e.j][kp.gi]!.Value).IsZero()));

        Console.WriteLine("All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : {0}",
            clggi.All(kp =>
                rg.Aggregate(e0.Zero, (sum, r) => sum + chis[r][kp.Key]!.Value * chis[r][kp.Value]!.Value)
                    .Equals(e0.One * Ocl[kp.Key])));
        Console.WriteLine("All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : {0}",
            allCombs.All(e =>
                rg.Aggregate(e0.Zero, (sum, r) => sum + chis[r][e.gi]!.Value * chis[r][clggi[e.gj]]!.Value).IsZero()));
    }

    private string[,] PrepareTable()
    {
        var idxs = NbClasses.Range().Grid2D().ToArray();
        var allCnfs = idxs.Where(e => !AllCharacters[e.t1].IsLinear)
            .Select(e => AllCharacters[e.t1][Classes.GetRepresentative(e.t2)])
            .Where(e => e.HasValue).Select(e => e!.Value.Simplify()).ToHashSet();
        var cnf2str = new Dictionary<Cnf, string>(NbClasses.Pow(2));

        foreach (var c in allCnfs)
        {
            if (!cnf2str.ContainsKey(c))
            {
                var (c0, c1, c2, c3) = FG.PrettyPrintCnf(c);
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
            var chi = AllCharacters[i];
            var c = chi[Classes.GetRepresentative(j)];
            if (!chi.IsLinear)
                cellStrs[i, j] = !c.HasValue ? "#" : $" {cnf2str[c.Value]}";
            else
                cellStrs[i, j] = !c.HasValue ? "#" : $" {c.Value}";
        }

        return cellStrs;
    }

    public void DisplayCells(bool tableOnly = false, char Chi = 'Ꭓ')
    {
        AllCharacters = AllCharacters.Order().ToArray();
        var form = Ring.MatrixDisplayForm;
        var digits = $"{NbClasses}".Length;
        var fmtChi = $"{Chi}.{{0,-{digits}}}";
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
                    Cells[i, j] = string.Format(fmtChi, i - 2);

                if (j > 1 && i > 2)
                    Cells[i, j] = StrTable[i - 3, j - 2];
            }
        }

        if (!tableOnly)
            DisplayGroup.Head(Gr);

        Console.WriteLine("Characters Table");
        Console.WriteLine();
        Ring.DisplayMatrix(Cells, " ");

        if (!tableOnly)
            CheckProperties();

        Console.WriteLine();

        Ring.MatrixDisplayForm = form;
    }
}