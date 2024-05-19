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
        var mtGL = FG.MetaCyclicSdpMat(e.m, e.n, e.r, maxDim);
        if (mtGL.Count() != 1)
        {
            mtGL.Name = $"M({e.m}x:{e.n}){e.r}";
            DisplayGroup.Head(mtGL);
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
                var tr = CharacterTable<Mat>.PrettyPrintCnf(mat.Trace.Simplify()).c;
                Console.WriteLine($"gen{k} of order {ct.Classes.GetClassName(gens[mat])} Tr={tr}");
                Console.WriteLine(gens[mat]);
                Console.WriteLine();
                Console.WriteLine(mat);
            }

            if (!mtGL.IsIsomorphicTo(mtGLnC))
                throw new();

            var m0 = mtGL.Neutral();
            var cl_tr = gens.ToDictionary(kv => ct.Classes.GetRepresentative(kv.Value), kv => kv.Key.Trace);
            var allChis = ct.AllCharacters.Order().Select((chi, k) => (chi, k: k + 1)).ToArray();
            var oneChis = allChis.Where(ck => (ck.chi[m0]! - n)!.Value.IsZero() &&
                                           cl_tr.All(kv => (kv.Value - ck.chi[kv.Key]!).Value.IsZero()))
                .ToArray();
            
            if (oneChis.Length == 0)
            {
                var twoChis = allChis.Grid2D().Where(ck => ck.t1.k < ck.t2.k)
                    .Select(ck => (chi: ck.t1.chi + ck.t2.chi, k1: ck.t1.k, k2: ck.t2.k))
                    .Where(ck => (ck.chi[m0]! - n)!.Value.IsZero() &&
                                 cl_tr.All(kv => (kv.Value - ck.chi[kv.Key]!).Value.IsZero()))
                    .ToArray();
                
                if (twoChis.Length == 0)
                    missing.Add(e);
                else
                    twoChis.Select(ck => $"Ꭓ.{ck.k1} + Ꭓ.{ck.k2} = {ck.chi}").Println("Two Characters");
            }
            else
                oneChis.Select(ck => $"Ꭓ.{ck.k} = {ck.chi}").Println("Characters");
            
            ct.DisplayCells(tableOnly: true);
            Console.WriteLine();
        }
    }

    var total = allMtCycSdp.Length;
    missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}",
        $"Characters Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    GlobalStopWatch.Show("END");
    Console.Beep();
}

Mat GLDiagOrdGenerators(GL gl, int ord)
{
    var id = gl.Neutral().Table;
    var e0 = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(gl.P, ord);
    return gl.At(id, 0, e0);
}

Mat[] GLDiagPermOrdGenerators(int n, int p, int ord)
{
    var gl = new GL(n, p);
    var diag = GLDiagOrdGenerators(gl, ord);
    var perms = FG.SnGensMat(n).Select(e => gl.Create(e.Table)).ToArray();
    return perms.Append(diag).ToArray();
}

ConcreteGroup<Mat> GLDiagPermOrd(int n, int p, int ord)
{
    var gens = GLDiagPermOrdGenerators(n, p, ord);
    return Group.Generate($"DPGL({n},{p})", gens[0].GL, gens);
}

ConcreteGroup<Mat> GLDiagPerm(int n, int p) => GLDiagPermOrd(n, p, p - 1);

{
    AllGensOfMtCycSdpUpToOrder(48);
}