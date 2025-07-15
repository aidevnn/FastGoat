using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Polynoms;

namespace CraftsAndExamples.Examples;

public static class CharacterTableExamplesPart2
{
    static CharacterTable<T> TensorTable<T>(CharacterTable<T> ctG, ConcreteGroup<T> H, ConcreteGroup<T> K)
        where T : struct, IElt<T>
    {
        if (ctG.TableComplete)
            return ctG;

        var G = ctG.Gr;
        var ctH = FG.CharacterTable(H);
        var ctK = FG.CharacterTable(K);
        if (H.Name == "SL(2,3)")
            ctH.SolveOrthogonality();
        if (K.Name == "SL(2,3)")
            ctK.SolveOrthogonality();
        
        var keys = H.Grid2D(K).Select(e => (h: e.t1, k: e.t2, hk: ctG.Classes.GetRepresentative(G.Op(e.t1, e.t2))))
            .DistinctBy(e => e.hk)
            .ToArray();

        foreach (var (chi1, chi2) in ctH.AllCharacters.Grid2D(ctK.AllCharacters))
        {
            var map = keys.ToDictionary(e => e.hk, e => (Cnf?)(chi1[e.h]!.Value * chi2[e.k]!.Value));
            var state = ctG.AddCharacter(new(ctG.Classes, map));
            if (state == AddCharacterState.TableFull)
                break;
        }

        if (Logger.Level != LogLevel.Off)
        {
            ctH.DisplayCells();
            ctK.DisplayCells();
        }

        return ctG;
    }

    static void TensorTable<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var gSubgrs = g.AllSubgroups();
        gSubgrs.Naming();
        var dirProd = gSubgrs.DecomposeProducts(gSubgrs.ProperNonTrivialNormalSubgroups())
            .Where(e => e.isDirectProduct).Take(1).ToArray();
        if (dirProd.Length == 0)
            return;

        var (H, K) = (dirProd[0].lhs.Representative, dirProd[0].rhs.Representative);

        Console.WriteLine($"#### {g.Name} = {H.Name} x {K.Name} ####");
        var ctG = FG.CharacterTableEmpty(g);
        ctG = TensorTable(ctG, H, K);
        ctG.DisplayCells();

        var ctG2 = FG.CharacterTableEmpty(g);
        ctG2.DerivedSubGroupLift();
        ctG2.InductionFromStabilizers();
        ctG2.InductionFromSubGroups(gSubgrs);
        ctG2.OrderCharacters();
        
        if (!ctG.AllCharacters.Zip(ctG2.AllCharacters).All(e=> e.First.Equals(e.Second)))
        {
            Console.WriteLine("----------------------------------------------------------------------------------------------");
            ctG2.DisplayCells();
            
            Console.WriteLine("#### Error");
            Console.Beep();
        }
        
        Console.WriteLine($"END {g.ShortName}");
        
        Console.WriteLine();
    }

    public static void Example1()
    {
        Logger.Level = LogLevel.Level1;
        TensorTable(FG.Abelian(2, 3));
        TensorTable(Product.Generate(FG.AbelianPerm(3), FG.Dihedral(4)));
    }

    public static void Example2()
    {
        // Logger.Level = LogLevel.Level1;
        foreach (var g in FG.AllGroupsOfOrder(1, 48).Where(g => g.GroupType == GroupType.NonAbelianGroup))
            TensorTable(g);
    }
}