using System.Text;
using Craft;
using Craft.Craft;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.Tools;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Perm.Style = DisplayPerm.CyclesComma;

XSet<Perm> InnerAut(Automorphism<Perm> aut)
{
    var G = aut.Domain;
    var act = Group.ByConjugate(G);
    var gens = G.GetGenerators().ToArray();
    var g0 = gens.OrderBy(e => e.Orbits.Length).First();
    if (!Perm.TypeEquals(g0, aut[g0]))
        return new();
    
    return UGCraft.InnerAut(g0, aut[g0]).Where(b => gens.All(x => act(b, x).Equals(aut[x]))).ToXSet();
}

ConcreteGroup<Perm> SdpWr(ConcreteGroup<Perm> N, Homomorphism<Perm, Automorphism<Perm>> theta,
    ConcreteGroup<Perm> H)
{
    if (theta.Kernel().Count() != 1)
        throw new();
    
    var snN = N.Neutral().Sn;
    var nGens = N.GetGenerators().ToArray();
    var hGens = H.GetGenerators().SelectMany(e => InnerAut(theta[e])).OrderBy(e => e.Order).ToArray();
    return Group.Generate("K", snN, nGens.Concat(hGens).ToArray());
}

ConcreteGroup<Perm> ProductPermGroup(params ConcreteGroup<Perm>[] Gs)
{
    var HnBase = Product.Gp(Gs.Cast<IGroup<Perm>>().ToArray());
    var HnGens = HnBase.GetGenerators().Select(e => FG.ConcatPerm(e.Ei)).ToArray();
    return Group.Generate(HnBase.Name, HnGens[0].Sn, HnGens);
}

void WreathProductPerm(ConcreteGroup<Perm> H, ConcreteGroup<Perm> S)
{
    var id0 = WreathProduct(H, S);
    var n = S.Neutral().Sn.N;
    var HnBase = Product.Gp(H, n);
    var HnGens = HnBase.GetGenerators().Select(e => FG.ConcatPerm(e.Ei)).ToArray();
    var Hn = Group.Generate(HnBase.Name, HnGens[0].Sn, HnGens);
    var autHn = Group.AutBase(Hn);
    var hom = S.ToDictionary(
        s => s,
        s => new Automorphism<Perm>(autHn,
            Group.AutomorphismMap(Hn,
                HnBase.GetGenerators().ToDictionary(ep => FG.ConcatPerm(ep.Ei), ep => FG.ConcatPerm(s.Apply(ep.Ei)))
            )
        )
    );

    var theta = new Homomorphism<Perm, Automorphism<Perm>>(S, hom);
    var wr = SdpWr(Hn, theta, S);
    wr.Name = $"{H.NameParenthesis()} wr {S.NameParenthesis()}";
    DisplayGroup.Generators(wr);

    var subgs = wr.AllSubgroups().ToGroupWrapper();
    var id1 = FG.FindIdGroup(wr, subgs.Infos);
    if (!id0.SequenceEqual(id1))
        throw new();
}

IdGroup[] WreathProduct<T>(ConcreteGroup<T> H, ConcreteGroup<Perm> S) where T : struct, IElt<T>
{
    var n = S.Neutral().Sn.N;
    var Hn = Product.GpGenerate(H, n);
    var autHn = Group.AutBase(Hn);
    var hom = S.ToDictionary(
        s => s,
        s => new Automorphism<Ep<T>>(autHn,
            Group.AutomorphismMap(Hn,
                Hn.GetGenerators().ToDictionary(ep => ep, ep => new Ep<T>(s.Apply(ep.Ei)))
            )
        )
    );
    var theta = new Homomorphism<Perm, Automorphism<Ep<T>>>(S, hom);
    var wr = Group.SemiDirectProd($"{H.NameParenthesis()} wr {S.NameParenthesis()}", Hn, theta, S);

    var subgs = wr.AllSubgroups().ToGroupWrapper();
    var names = NamesTree.BuildName(subgs);
    FG.DisplayName(subgs.Parent, subgs, names);
    return FG.FindIdGroup(wr, subgs.Infos);
}

{
    WreathProductPerm(FG.AbelianPerm(2), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(3), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(2), FG.AbelianPerm(3));
    WreathProductPerm(FG.AbelianPerm(4), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(2, 2), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(2), FG.Symmetric(3));
    WreathProductPerm(FG.AbelianPerm(5), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(2), FG.AbelianPerm(4));
    WreathProductPerm(FG.AbelianPerm(2), FG.AbelianPerm(2, 2)); // [64,226]
    var V = FG.PermGroup("C2 x C2", 4, ((1, 2), (3, 4)), ((1, 3), (2, 4)));
    WreathProductPerm(FG.AbelianPerm(2), V); // [64,138]
    WreathProductPerm(FG.AbelianPerm(6), FG.AbelianPerm(2));
    WreathProductPerm(FG.Symmetric(3), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(3), FG.AbelianPerm(3));
    WreathProductPerm(FG.AbelianPerm(7), FG.AbelianPerm(2));
}

