using System.Globalization;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using System.Numerics;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.Arm;
using System.Security.Cryptography.X509Certificates;
using System.Xml;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

(EPoly<K> W, int l) Primitive3<K>(EPoly<K> U, EPoly<K> V) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var n = U.F.Degree;
    var vecV = V.Poly.ToVMatrix(n);
    for (int l = 1; l < 50; l++)
    {
        var W = U + l * V;
        var M = KMatrix<K>.MergeSameRows(n.Range().Select(i => W.Pow(i).Poly.ToVMatrix(n)).ToArray());
        var vM = KMatrix<K>.MergeSameRows(vecV, M);
        var dimKerM = M.NullSpace().nullity;
        var dimKerVM = vM.NullSpace().nullity;
        if (dimKerM == dimKerVM)
        {
            return (W, l);
        }
    }

    throw new();
}

// {
//     Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
//     var x = FG.QPoly();
//     var (z, i) = FG.NumberFieldQ((x.Pow(2) - 3, "√3"), (x.Pow(2) + 1, "i"));
//     var one = z.One;
//     var gl = FG.GLnK($"Q({z})", 2, z);
//     var j = (1 + z * i) / 2;
//     var j2 = (1 - z * i) / 2;
//     
//     var A = gl[j, 0, 0, j2];
//     var B = gl[0, 1, 1, 0];
//     var C = gl[one / 2, z / 2, -z / 2, one / 2];
//     var D = gl[1, 0, 0, -1];
//     var GAB = Group.Generate("D12-AB", gl, A, B);
//     var GCD = Group.Generate("D12-CD", gl, C, D);
//     DisplayGroup.HeadElements(GAB);
//     DisplayGroup.HeadElements(GCD);
//
//     var d12 = FG.DihedralWg(6);
//     var a = d12["a"];
//     var b = d12["b"];
//     DisplayGroup.HeadElements(d12);
//
//     var pMapAB = Group.PartialMap((a, A), (b, B));
//     var homAB = Group.Hom(d12, Group.HomomorphismMap(d12, GAB, pMapAB));
//     Console.WriteLine(homAB.HomMap.Glue("\n"));
//     DisplayGroup.AreIsomorphics(d12, GAB);
//     Console.WriteLine();
//
//     var pMapCD = Group.PartialMap((a, C), (b, D));
//     var homCD = Group.Hom(d12, Group.HomomorphismMap(d12, GCD, pMapCD));
//     Console.WriteLine(homCD.HomMap.Glue("\n"));
//
//     DisplayGroup.AreIsomorphics(GAB, GCD);
// }

{
    Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;

    var rg10 = 9.Range(2);
    var metacycl = rg10.Grid2D(rg10)
        .SelectMany(e => SolveAll_k_pow_m_equal_one_mod_n(e.t1, e.t2).Select(r => FG.MetaCyclicSdp(e.t1, e.t2, r)))
        .OrderBy(e => e.Count()).ThenBy(e => e.ShortName).ToArray();
    
    var gr = (SemiDirectProduct<ZnInt, ZnInt>)metacycl.First(e=>e.Count() == 20);
    var m = gr.G.Count();
    var n = gr.N.Count();
    var am = new Cnf(m);
    var an = new Cnf(n);
    Console.WriteLine(new { m, n, am, an });
    Group.DisplayConjugacyClasses(gr);
    var cl = Group.AllConjugacyClasses(gr);
    var orbx = Group.AllOrbits(gr, Group.ByConjugate(gr));
    var dg = Group.DerivedChain(gr);
    var quo = dg[0].Over(dg[1]);
    
    dg.Select(g => g.ShortName).Println();
    DisplayGroup.HeadElements(quo);
    
    orbx.ToDictionary(e => e.Key, e => e.Value.Orbx.Glue("; ")).Println();
    var Bis = orbx.ToDictionary(e => e.Key, e => e.Value.Orbx.Aggregate(gr.Neutral(), (sum, g) => gr.Op(sum, g)));
    Bis.Println();
    
    var glk = FG.GLnK($"CF({m * n}", 2, am);
    var mat0 = glk[am, 1, 0, an];
    Int32.Max(m, n).Range().Select(i => am.Pow(i)).Println();
    Int32.Max(m, n).Range().Select(i => an.Pow(i)).Println();

    var det = mat0.Det;
    Console.WriteLine(det);
    Console.WriteLine(det.Conj);
    Console.WriteLine(det.Re);
    Console.WriteLine(det.Im);
    Console.WriteLine(det.Module2);

    Console.WriteLine(mat0);
    Console.WriteLine(mat0.Pow(2));
    Console.WriteLine(mat0.Pow(3));
}
