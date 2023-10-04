using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Gp<Tg> CreatePr<Tg>(ConcreteGroup<Tg> G, int r) where Tg : struct, IElt<Tg>
{
    return Product.Gp(Enumerable.Repeat(G, r + 1).Cast<IGroup<Tg>>().ToArray());
}

HashSet<Ep<Tg>> GetBasis<Tg>(ConcreteGroup<Tg> G, int r) where Tg : struct, IElt<Tg>
{
    var n0 = G.Neutral();
    var set = G.MultiLoop(r).Select(l => new Ep<Tg>(l.Prepend(n0).ToArray())).Order().ToHashSet();
    set.Println($"P{r} of Z[{G}]-module elts:{set.Count}");
    return set;
}

HashSet<Ep<Tg>> GetAll<Tg>(ConcreteGroup<Tg> G, int r) where Tg : struct, IElt<Tg>
{
    var set = G.MultiLoop(r + 1).Select(l => new Ep<Tg>(l.ToArray())).Order().ToHashSet();
    set.Println($"P{r} of Z[{G}]-module elts:{set.Count}");
    return set;
}

HashSet<Ep<Tg>> ApplyDr<Tg>(ConcreteGroup<Tg> G, HashSet<Ep<Tg>> elts) where Tg : struct, IElt<Tg>
{
    var set = new HashSet<Ep<Tg>>();
    var rs = elts.Select(e => e.Ei.Length).Distinct().ToArray();
    if (rs.Length != 1)
        throw new();

    var r = rs[0] - 1;
    var gr0 = CreatePr(G, r - 1);
    foreach (var ep in elts)
    {
        var n0 = gr0.Neutral();
        for (int i = 0; i <= r; i++)
        {
            var ep0 = ep.SkipAt(i);
            n0 = gr0.Op(n0, gr0.Times(ep0, (-1).Pow(i)));
        }

        Console.WriteLine($"d{r}({ep}) = {n0}");
        set.Add(n0);
    }

    Console.WriteLine();
    return set;
}

{
    var (c2, c4, c2c2, c5, c2c2c2, c4c4) = (new Cn(2), new Cn(4), FG.Abelian(2, 2), new Cn(5), FG.Abelian(2, 2, 2), FG.Abelian(4, 4));
    
    {
        var (G, r) = (c2, 2);
        var set = GetAll(G, r);
        set = ApplyDr(G, set);
        set = ApplyDr(G, set);
    }

    {
        var (G, r) = (c4, 2);
        var set = GetAll(G, r);
        set = ApplyDr(G, set);
        set = ApplyDr(G, set);
    }
    
    // {
    //     var (G, r) = (c2c2c2, 2);
    //     var set = GetAll(G, r);
    //     set = ApplyDr(G, set);
    //     set = ApplyDr(G, set);
    // }
}