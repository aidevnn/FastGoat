using System.CodeDom;
using System.Collections;
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
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

(int m, int n, int r)[] MetaCyclicSdp(int order)
{
    return IntExt.Dividors(order).Where(d => d > 1)
        .SelectMany(m => FG.MetaCyclicSdpGetR(m, order / m).Select(r => (m, n: order / m, r)))
        .ToArray();
}

void AllGensOfMtCycSdpUpToOrder(int maxOrd, int maxDim = 6)
{
    GlobalStopWatch.Restart();
    var missing = new List<(int, int, int)>();
    var allMtCycSdp = (maxOrd - 5).Range(6).SelectMany(ord => MetaCyclicSdp(ord)).ToArray();

    foreach (var e in allMtCycSdp)
    {
        var mt = FG.MetaCyclicSdpWg(e.m, e.n, e.r);
        var mtGL = FG.MetaCyclicSdpMat(e.m, e.n, e.r, maxDim);
        if (mtGL.Count() != 1)
        {
            mtGL.Name = $"M({e.m}x:{e.n}){e.r}";
            var n = mtGL.Neutral().GL.N;
            var p = mtGL.Neutral().GL.P;
            var Up = FG.UnInt(p);
            var e0 = Up.GetGenerators().First();
            var cnf = Cnf.Nth(p - 1);
            var GLnC = FG.GLnK("C", n, cnf);
            var iso = (p - 1).Range().ToDictionary(k => e0.Pow(k).K, k => cnf.Pow(k).Simplify());
            iso[0] = cnf.Zero;
            var gens = mtGL.GetGenerators().ToDictionary(mat => mat.Table.Select(z => iso[z]).ToKMatrix(n), mat => mat);
            var mtGLnC = Group.Generate(mtGL.Name, GLnC, gens.Keys.ToArray());
            DisplayGroup.HeadOrders(mtGLnC);
            
            var ct = FG.CharacterTable(mtGL);
            foreach (var (mat, k) in mtGLnC.GetGenerators().Select((mat, k) => (mat, k + 1)))
            {
                Console.WriteLine(
                    $"gen{k} of order {ct.Classes.GetClassName(gens[mat])} Tr={IntFactorisation.PrettyPrintCnf(mat.Trace.Simplify())}");
                Console.WriteLine(gens[mat]);
                Console.WriteLine();
                Console.WriteLine(mat);
            }

            if (!mtGL.IsIsomorphicTo(mtGLnC))
                throw new();

            var cl_tr = gens.ToDictionary(kv => ct.Classes.GetRepresentative(kv.Value), kv => kv.Key.Trace);
            var chis = ct.AllCharacters
                .Where(chi => cl_tr.All(kv => (kv.Value - chi[kv.Key]!).Value.Simplify().IsZero()))
                .ToArray();
            if (chis.Length == 0)
            {
                if (n > ct.AllCharacters.Max(chi => chi[mtGL.Neutral()]!.Value.E[0].Num))
                    Console.WriteLine("Diagonal block representation");
                else
                {
                    cl_tr.ToDictionary(kv => ct.Classes.GetClassName(kv.Key), kv => kv.Value)
                        .Println("Classe and Trace");
                    Console.WriteLine();
                    foreach (var chi in ct.AllCharacters)
                    {
                        cl_tr.Select(kv => (ct.Classes.GetClassName(kv.Key), kv.Value, chi[kv.Key],
                                (kv.Value - chi[kv.Key])!.Value.Simplify()))
                            .Println($"Chi:{chi}");
                    }

                    throw new();
                }
            }
            else
                chis.Println("Characters");
            
            ct.DisplayCells(tableOnly: true);

            var hom = Group.AllHomomorphisms(mt, mtGL);
            Console.WriteLine($"Nb Endo({mtGL}):{hom.Count}");
            Console.WriteLine();
            continue;
        }

        missing.Add(e);
    }

    var total = allMtCycSdp.Length;
    missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}",
        $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    GlobalStopWatch.Show("END");
    Console.Beep();
}

{
    AllGensOfMtCycSdpUpToOrder(32);
}