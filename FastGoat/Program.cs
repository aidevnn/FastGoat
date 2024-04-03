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

Perm.Style = DisplayPerm.CyclesComma;

Ep<ZnInt> Image(Perm g, Ep<ZnInt> x) => new(x.Ei.Select(e => new ZnInt(g.Sn.N, g.Table[e.K])).ToArray());
HashSet<Ep<ZnInt>> Images(Perm g, Ep<ZnInt>[] Y) => Y.Select(y => Image(g, y)).ToHashSet();

void PrimitivesTransitivesOfSn(int n)
{
    var sn = new Symm(n);
    Console.WriteLine($"Transitives groups of order {n}");
    var Xn = FG.Abelian(n).ToArray();
    var allSubgSn = sn.AllSubgroups();
    allSubgSn.Naming();
    var Ys = Xn.AllCombinationsFromM(2).Where(y => y.Length != n).ToArray();
    var dic = new[] { "Primitive", "Imprimitive" }.Grid2D(["Transitve", "NonTransitive"])
        .ToDictionary(e => e, _ => new List<SubgroupConjugates<Perm>>());
    foreach (var g1 in allSubgSn.AllRepresentatives)
    {
        var allOrbits = Group.AllOrbits(g1, Xn, Image);
        var blocks = Ys.Where(Y => g1.All(g =>
        {
            var Yg = Images(g, Y);
            return Yg.SetEquals(Y) || !Yg.Overlaps(Y);
        })).ToArray();
        var isPrimitve = blocks.Length == 0 ? "Primitive" : "Imprimitive";
        var isTransitive = allOrbits.Count == 1 ? "Transitve" : "NonTransitive";
        var cj = allSubgSn.AllSubgroupConjugates.First(cj => cj.Contains(g1));
        dic[(isPrimitve, isTransitive)].Add(cj);
    }

    dic.Select(e => ($"{e.Key} => {e.Value.Count}/{e.Value.Sum(cj => cj.Size)}")).Println("Count");
    dic[("Primitive", "Transitve")].Println(e => e.FullName, "Primitive, Transitve");
    Console.WriteLine();
}

void Test()
{
    PrimitivesTransitivesOfSn(3);
    PrimitivesTransitivesOfSn(4);
    PrimitivesTransitivesOfSn(5);
    PrimitivesTransitivesOfSn(6);
    PrimitivesTransitivesOfSn(7);
    
    Console.Beep();
}

void WreathProduct<T>(ConcreteGroup<T> H, ConcreteGroup<Perm> S) where T : struct, IElt<T>
{
    var n = S.Neutral().Sn.N;
    var Hn = Product.GpGenerate(H, n);
    var autHn = Group.AutBase(Hn);
    var op = S.ToDictionary(s => s, s => new Automorphism<Ep<T>>(autHn, Hn.ToDictionary(ep => ep, ep => new Ep<T>(s.Apply(ep.Ei)))));
    var theta = new Homomorphism<Perm, Automorphism<Ep<T>>>(S, op);
    var wr = Group.SemiDirectProd($"{H.NameParenthesis()} wr {S.NameParenthesis()}", Hn, theta, S);
    
    FG.DisplayBox(wr.AllSubgroups(), 0);
    // var subgs = wr.AllSubgroups().ToGroupWrapper();
    // var names = NamesTree.BuildName(subgs);
    // FG.DisplayName(subgs.Parent, subgs, names);
}

{
    WreathProduct(FG.AbelianPerm(2), FG.AbelianPerm(2));
    WreathProduct(FG.AbelianPerm(3), FG.AbelianPerm(2));
    WreathProduct(FG.AbelianPerm(2), FG.AbelianPerm(3));
    WreathProduct(FG.AbelianPerm(4), FG.AbelianPerm(2));
    WreathProduct(FG.AbelianPerm(2, 2), FG.AbelianPerm(2));
    WreathProduct(FG.AbelianPerm(2), FG.Symmetric(3));
    WreathProduct(FG.AbelianPerm(5), FG.AbelianPerm(2));
    WreathProduct(FG.AbelianPerm(2), FG.AbelianPerm(4));
    WreathProduct(FG.AbelianPerm(2), FG.AbelianPerm(2, 2)); // [64,226]
    var V = FG.PermGroup("C2 x C2", 4, ((1, 2), (3, 4)), ((1, 3), (2, 4)));
    WreathProduct(FG.AbelianPerm(2), V); // [64,138]
    WreathProduct(FG.AbelianPerm(6), FG.AbelianPerm(2));
    WreathProduct(FG.Symmetric(3), FG.AbelianPerm(2));
    WreathProduct(FG.AbelianPerm(3), FG.AbelianPerm(3));
    WreathProduct(FG.AbelianPerm(7), FG.AbelianPerm(2));
}