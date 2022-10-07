using FastGoat;
using FastGoat.Gp;
using FastGoat.UserGroup;
using static FastGoat.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

List<Dictionary<T2, Automorphism<T1>>> GroupAction<T1, T2>(IGroup<T1> bn, IGroup<T2> bg)
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
    var nGens = bn.GetGenerators().ToArray();
    var gGens = bg.GetGenerators().ToArray();
    var n = Group.Generate(bn, nGens);
    var g = Group.Generate(bg, gGens);

    var autN = Group.Aut(bn, nGens);
    var autOrders = autN.GroupBy(e => autN.ElementsOrders[e]).Select(e => (ord: e.Key, auts: e.ToArray())).ToArray();
    var gGensOrders = gGens.Select(e => (g: e, ord: g.ElementsOrders[e])).ToArray();
    var gpMap = gGensOrders.Select(e =>
            autOrders.Where(a => a.ord != 1 && e.ord % a.ord == 0).SelectMany(a => a.auts).Select(a => (e.g, a))
                .ToArray())
        .ToArray();


    var maps = new List<Dictionary<T2, Automorphism<T1>>>() { new Dictionary<T2, Automorphism<T1>>() };
    foreach (var tuples in gpMap)
    {
        foreach (var e in tuples)
        {
            var tmpMaps = new List<Dictionary<T2, Automorphism<T1>>>();
            foreach (var map in maps)
            {
                var pmap = new Dictionary<T2, Automorphism<T1>>(map);
                pmap[e.g] = e.a;
                tmpMaps.Add(pmap);
            }

            maps.Clear();
            maps = tmpMaps.ToList();
        }
    }

    List<Dictionary<T2, Automorphism<T1>>> allHom = new();
    var ng = g.Count();
    foreach (var map in maps)
    {
        var hom = Group.HomomorphismMap(g, autN, map);
        if (hom.Count == ng)
            allHom.Add(hom);
    }

    return allHom;
}

// {
//     var N = 102;
//     var lt = Enumerable.Range(2, N).ToArray();
//     var tuples = lt.SelectMany(e1 => lt.Where(e2 => e1 * e2 <= N).Select(e2 => (e1, e2))).ToArray();
//     List<(int n, int g, Dictionary<ZnInt, Automorphism<ZnInt>> y)> solutions = new();
//
//     foreach (var (n0, g0) in tuples)
//     {
//         var n = new Zn(n0);
//         var g = new Zn(g0);
//         var allHom = GroupAction(n, g);
//         if (allHom.Count == 0)
//         {
//             try
//             {
//                 var sdp = Group.SemiDirectProd(new Cn(n0), new Cn(g0));
//             }
//             catch (Exception e)
//             {
//                 Console.WriteLine($"GOOD!!! {e}");
//             }
//         }
//         else
//         {
//             try
//             {
//                 var sdp = Group.SemiDirectProd(new Cn(n0), new Cn(g0));
//             }
//             catch (Exception e)
//             {
//                 Console.WriteLine($"BAD!!! {e}");
//             }
//         }
//         
//         solutions.AddRange(allHom.Select(a => (n0, g0, a)));
//     }
//
//     Console.WriteLine("All solutions : {0}", solutions.Count);
//     foreach (var s in solutions.OrderBy(s0 => s0.g * s0.n))
//     {
//         var cn = new Cn(s.n);
//         var cg = new Cn(s.g);
//         var sdp = Group.SemiDirectProd(cn, cg);
//         Console.WriteLine("Size:{0,-4}   C{1,-2} : C{2,-6} {3}", s.n * s.g, s.n, s.g, sdp.Count() == s.n * s.g);
//         // foreach (var kp in s.y)
//         // {
//         //     Console.WriteLine($"    g={kp.Key} y=({kp.Value})");
//         // }
//
//         // Console.WriteLine();
//     }
// }

{
    var n = Product.Group(new Zn(2), new Zn(2), new Zn(2));
    var g = Product.Group(new Zn(2), new Zn(2));
    var allHom = GroupAction(n, g);
    Console.WriteLine(allHom.Count);
}

{
    var n = Product.Group(new Zn(2), new Zn(2));
    var g = Product.Group(new Zn(2), new Zn(2), new Zn(2));
    var allHom = GroupAction(n, g);
    Console.WriteLine(allHom.Count);
}

{
    var n = Product.Group(new Zn(3), new Zn(3));
    var g = Product.Group(new Zn(2), new Zn(4));
    var allHom = GroupAction(n, g);
    Console.WriteLine(allHom.Count);
}
