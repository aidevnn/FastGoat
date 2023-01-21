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

{
    var x = FG.QPoly();
    AlgebraicFactorization.Resolvent(x.Pow(4) + 9);
}