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

IEnumerable<Perm> InnerAutMatrix(Automorphism<Perm> aut)
{
    var G = aut.Domain;
    var sn = G.Neutral().Sn;
    var n = sn.N;
    var xis = Ring.Polynomial(ZnInt.ZnZero(2), (n * n).SeqLazy(1).Select(i => $"x{i}").ToArray());
    var Mx = xis.ToKMatrix(n);
    var rowsM = Mx.Rows.Select(r => r.ToXSet()).ToArray();
    var colsM = Mx.Cols.Select(c => c.ToXSet()).ToArray();
    var sys = G.GetGenerators().Select(a =>
        {
            var b = aut[a];
            var Ma = MatrixExt.Permutation(a.Table).Select(i => i * Mx.KOne).ToKMatrix(n);
            var Mb = MatrixExt.Permutation(b.Table).Select(i => i * Mx.KOne).ToKMatrix(n);
            return Mx * Ma - Mb * Mx;
        }).SelectMany(m => m.Where(e => !e.IsZero()))
        .ToHashSet();

    var bagInit = sys.Select(e => e.Variables.ToXSet()).ToHashSet();
    var bagEnd = bagInit.Take(0).ToHashSet();
    while (bagInit.Count != 0)
    {
        var set = bagInit.Max();
        var inter = bagInit.Where(s => s.Overlaps(set)).ToArray();
        bagInit.ExceptWith(inter);

        if (inter.Length == 1)
            bagEnd.Add(set);
        else
            bagInit.Add(inter.Union());
    }

    var dependantXis = bagEnd.ToDictionary(e => e.Min(), e => e);
    var zeros = dependantXis.Where(e =>
            rowsM.Any(r => r.Intersect(e.Value).Count() > 1) || colsM.Any(c => c.Intersect(e.Value).Count() > 1))
        .SelectMany(e => e.Value).ToXSet();
    var M1 = Mx.Select(e => zeros.Contains(e) ? e.Zero : e).ToKMatrix(n);
    if (Logger.IsOn)
    {
        Console.WriteLine($"aut:{aut}");
        Console.WriteLine(Mx);
        dependantXis.OrderBy(e => e.Value)
            .Println(e => $"IsZero:{zeros.Contains(e.Key),-6} {e.Key} => {e.Value}", "Solve x*a = aut(a)*x");
        Console.WriteLine(M1);
    }
    
    if (M1.Rows.Any(r => r.All(e => e.IsZero())) || M1.Cols.Any(c => c.All(e => e.IsZero())))
        yield break;
    
    var pos = n.Range().Grid2D().ToDictionary(e => Mx[e.t1, e.t2], e => e);
    var blocks = dependantXis.Where(e => !zeros.Contains(e.Key))
        .ToDictionary(e => e.Value, e =>
        {
            var coords = e.Value.Select(xi => pos[xi]).ToArray();
            var tl = (coords.Min(p => p.t1), coords.Min(p => p.t2));
            var br = (coords.Max(p => p.t1), coords.Max(p => p.t2));
            return (tl, br);
        })
        .GroupBy(e => e.Value).ToDictionary(e => e.Key, e => e.Select(f => f.Key).ToXSet());
    
    if (Logger.IsOn)
        blocks.Println("Blocks");
    var rg = n.Range();
    var act = Group.ByConjugate(sn);
    foreach (var perms in blocks.Values.Select(l => l.ToArray()).MultiLoop().Select(l => l.ToArray()))
    {
        var ones = perms.SelectMany(e => e).ToHashSet();
        var M2 = M1.Select(e => ones.Contains(e) ? e.KOne : e.KZero).ToKMatrix(n);
        if (Logger.IsOn)
        {
            Console.WriteLine($"perms:{perms.ToXSet()}");
            Console.WriteLine(M2);
        }
        var arr = rg.Select(i => rg.FirstOrDefault(j => M2[i, j].IsOne(), 0)).ToArray();
        if (arr.Order().SequenceEqual(rg))
        {
            var x = sn.CreateElementTable(arr);
            if (Logger.IsOn)
            {
                Console.WriteLine($"    x = {x}");
                Console.WriteLine();
            }
            if (G.GetGenerators().All(a => act(x, a).Equals(aut[a])))
                yield return x;
            else
                throw new(); // never reached
        }
        else
        {
            if (Logger.IsOn)
                Console.WriteLine($"bad-arr:[{arr.Glue(", ")}]");
        }
    }
}

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
    var inn = H.GetGenerators().Select(e => InnerAutMatrix(theta[e]).Select(f => (e, f)).ToArray()).ToArray();
    var subGroupAutN = Group.Generate($"SubGr-Aut[{N}]", snN, inn.SelectMany(e => e.Select(a => a.f)).ToArray());
    if (subGroupAutN.Order == N.Order * H.Order && subGroupAutN.SuperSetOf(N))
        return subGroupAutN;

    DisplayGroup.HeadOrders(subGroupAutN);
    Console.WriteLine();

    var hom = Group.Hom(H,
        inn.MultiLoop().Select(e => Group.HomomorphismMap(H, subGroupAutN, e.ToDictionary(f => f.e, f => f.f)))
            .First(e => e.Count(f => f.Value.Order == 1) == 1)
    );

    var hGens = H.GetGenerators().Select(e => hom[e]).ToArray();
    var sdp = Group.Generate("K", snN, nGens.Concat(hGens).ToArray());
    Console.WriteLine($"inter:{sdp.Intersect(subGroupAutN).Count()}");
    Console.WriteLine();
    return sdp;
}

ConcreteGroup<Perm> ProductPermGroup(params ConcreteGroup<Perm>[] Gs)
{
    var HnBase = Product.Gp(Gs.Cast<IGroup<Perm>>().ToArray());
    var HnGens = HnBase.GetGenerators().Select(e => FG.ConcatPerm(e.Ei)).ToArray();
    return Group.Generate(HnBase.Name, HnGens[0].Sn, HnGens);
}

void WreathProductPerm(ConcreteGroup<Perm> H, ConcreteGroup<Perm> S)
{
    // var (wr0, id0) = WreathProduct(H, S);
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
    DisplayGroup.HeadOrdersGenerators(wr);
    var exportGap = wr.GetGenerators()
        .OrderBy(e => e.DisjoinCycles.Length).ThenBy(e => e)
        .Select(p => p.DisjoinCycles.Glue()).Glue(", ")
        .Replace("[", "").Replace("]", "");
    Console.WriteLine("ExportGap.");
    Console.WriteLine($"Group([{exportGap}]);");
    Console.WriteLine();

    // var subgs = wr.AllSubgroups().ToGroupWrapper();
    // var id1 = FG.FindIdGroup(wr, subgs.Infos);
    // if (!id0.SequenceEqual(id1))
    //     throw new();
}

(SemiDirectProduct<Ep<T>, Perm> wr, IdGroup[]) WreathProduct<T>(ConcreteGroup<T> H, ConcreteGroup<Perm> S)
    where T : struct, IElt<T>
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

    var subgs = wr.AllSubgroups();
    FG.DisplayBox(subgs, 0);
    return (wr, FG.FindIdGroup(wr, subgs.Infos));
}

void ExamplesWreathProduct()
{
    GlobalStopWatch.Restart();

    // order 8
    WreathProductPerm(FG.AbelianPerm(2), FG.AbelianPerm(2));

    // order 18
    WreathProductPerm(FG.AbelianPerm(3), FG.AbelianPerm(2));

    // order 24
    WreathProductPerm(FG.AbelianPerm(2), FG.AbelianPerm(3));

    // order 32
    WreathProductPerm(FG.AbelianPerm(4), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(2, 2), FG.AbelianPerm(2));

    // order 48
    WreathProductPerm(FG.AbelianPerm(2), FG.Symmetric(3));

    // order 50
    WreathProductPerm(FG.AbelianPerm(5), FG.AbelianPerm(2));

    // order 64
    WreathProductPerm(FG.AbelianPerm(2), FG.AbelianPerm(4));
    WreathProductPerm(FG.AbelianPerm(2), FG.AbelianPerm(2, 2)); // [64,226]
    var V = FG.PermGroup("C2 x C2", 4, ((1, 2), (3, 4)), ((1, 3), (2, 4)));
    WreathProductPerm(FG.AbelianPerm(2), V); // [64,138]

    // order 72
    WreathProductPerm(FG.AbelianPerm(6), FG.AbelianPerm(2));
    WreathProductPerm(FG.Symmetric(3), FG.AbelianPerm(2));

    // order 81
    WreathProductPerm(FG.AbelianPerm(3), FG.AbelianPerm(3));

    // order 98
    WreathProductPerm(FG.AbelianPerm(7), FG.AbelianPerm(2));

    // order 128
    WreathProductPerm(FG.AbelianPerm(8), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(2, 4), FG.AbelianPerm(2));
    WreathProductPerm(FG.Dihedral(4), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(2), FG.Dihedral(4));
    WreathProductPerm(FG.QuaternionPg(8), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(2, 2, 2), FG.AbelianPerm(2));

    // order 160
    WreathProductPerm(FG.AbelianPerm(2), FG.AbelianPerm(5));

    // order 162
    WreathProductPerm(FG.AbelianPerm(9), FG.AbelianPerm(2));
    WreathProductPerm(FG.AbelianPerm(3), FG.Symmetric(3));
    WreathProductPerm(FG.AbelianPerm(3, 3), FG.AbelianPerm(2));
    // |(C3 x C3) wr C2| = 162
    // Type        NonAbelianGroup
    // BaseGroup   S12
    // 
    // Elements Orders : [1]:1, [2]:9, [3]:80, [6]:72
    // Generators of (C3 x C3) wr C2
    // gen1 of order 2
    // [(1, 7), (2, 8), (3, 9), (4, 10), (5, 11), (6, 12)]
    // gen2 of order 6
    // [(1, 7), (2, 8), (3, 9), (4, 10, 5, 11, 6, 12)]
    // gen3 of order 6
    // [(1, 7, 2, 8, 3, 9), (4, 10), (5, 11), (6, 12)]
    // 
    // ExportGap.
    // Group([(1, 7)(2, 8)(3, 9)(4, 10, 5, 11, 6, 12), (1, 7, 2, 8, 3, 9)(4, 10)(5, 11)(6, 12), (1, 7)(2, 8)(3, 9)(4, 10)(5, 11)(6, 12)]);

    // order 192
    WreathProductPerm(FG.AbelianPerm(4), FG.AbelianPerm(3));
    WreathProductPerm(FG.AbelianPerm(2), FG.Alternate(4));
    WreathProductPerm(FG.AbelianPerm(2, 2), FG.AbelianPerm(3));

    GlobalStopWatch.Show("End"); // # End Time:4.859s
}

{
    Logger.SetOff();
    ExamplesWreathProduct();
    
    Logger.SetLevel1();
    WreathProductPerm(FG.AbelianPerm(2), FG.AbelianPerm(2));
    // aut:[(3, 4)]->[[(1, 2)]]; [(1, 2)]->[[(3, 4)]]
    // [ x1,  x2,  x3,  x4]
    // [ x5,  x6,  x7,  x8]
    // [ x9, x10, x11, x12]
    // [x13, x14, x15, x16]
    // Solve x*a = aut(a)*x
    //     IsZero:False  x14 => { x14, x9 }
    //     IsZero:False  x13 => { x13, x10 }
    //     IsZero:False  x8 => { x8, x3 }
    //     IsZero:False  x7 => { x7, x4 }
    //     IsZero:True   x16 => { x16, x15, x12, x11 }
    //     IsZero:True   x6 => { x6, x5, x2, x1 }
    // [  0,   0, x3, x4]
    // [  0,   0, x7, x8]
    // [ x9, x10,  0,  0]
    // [x13, x14,  0,  0]
    // Blocks
    //     [((0, 2), (1, 3)), { { x8, x3 }, { x7, x4 } }]
    //     [((2, 0), (3, 1)), { { x14, x9 }, { x13, x10 } }]
    // perms:{ { x13, x10 }, { x7, x4 } }
    // [0, 0, 0, 1]
    // [0, 0, 1, 0]
    // [0, 1, 0, 0]
    // [1, 0, 0, 0]
    //     x = [(1, 4), (2, 3)]
    // 
    // perms:{ { x14, x9 }, { x7, x4 } }
    // [0, 0, 0, 1]
    // [0, 0, 1, 0]
    // [1, 0, 0, 0]
    // [0, 1, 0, 0]
    //     x = [(1, 4, 2, 3)]
    // 
    // perms:{ { x13, x10 }, { x8, x3 } }
    // [0, 0, 1, 0]
    // [0, 0, 0, 1]
    // [0, 1, 0, 0]
    // [1, 0, 0, 0]
    //     x = [(1, 3, 2, 4)]
    // 
    // perms:{ { x14, x9 }, { x8, x3 } }
    // [0, 0, 1, 0]
    // [0, 0, 0, 1]
    // [1, 0, 0, 0]
    // [0, 1, 0, 0]
    //     x = [(1, 3), (2, 4)]
    // 
    // |C2 wr C2| = 8
    // Type        NonAbelianGroup
    // BaseGroup   S4
    // 
    // Elements Orders : [1]:1, [2]:5, [4]:2
    // Generators of C2 wr C2
    // gen1 of order 2
    // [(1, 3), (2, 4)]
    // gen2 of order 4
    // [(1, 3, 2, 4)]
    // 
    // ExportGap.
    // Group([(1, 3, 2, 4), (1, 3)(2, 4)]);
    // 
}