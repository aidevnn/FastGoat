using System.Numerics;
using System.Text;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

KMatrix<K> Babai<K>(KMatrix<K> B, KMatrix<K> W) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
{
    if (W.IsRow)
        return Babai(B.T, W.T).T;
    
    if (!B.IsSquare || !W.IsCol || W.M != B.M)
        throw new($"dim(B) = {B.Dim} dim(W) = {W.Dim}");

    var (Bi, Bsi) = (B.Cols, Ring.GramSchmidt(B).O.Cols);
    var v = Bsi[0].Zero;
    var Wi = W.Clone;
    for (int i = W.M - 1; i >= 0; i--)
    {
        var bsi = Bsi[i];
        var li = (Wi.T * bsi)[0, 0] / (bsi.T * bsi)[0, 0];
        var rli = li.RoundEven;
        v += rli * Bi[i];
        if (i > 0)
            Wi -= (li - rli) * bsi + rli * Bi[i];
    }

    return v;
}

KMatrix<K> BabaiRound<K>(KMatrix<K> B, KMatrix<K> W) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
{
    if (W.IsRow)
        return BabaiRound(B.T, W.T).T;
    
    if (!B.IsSquare || !W.IsCol || W.M != B.M)
        throw new($"dim(B) = {B.Dim} dim(W) = {W.Dim}");

    return Ring.ReducedRowsEchelonForm(KMatrix<K>.MergeSameRows(B, W)).A0.GetCol(B.N)
        .Select((e, i) => e.RoundEven * B.GetCol(i)).Aggregate((bi, bj) => bi + bj);
}

{
    var o = BigReal.BrOne(20);
    var B = new[] { 1, 2, 3, 3, 0, -3, 3, -7, 3 }.Select(e => e * o).ToKMatrix(3);
    Console.WriteLine("B");
    Console.WriteLine(B);
    
    var W = new[] { 10, 6, 5 }.Select(e => e * o).ToKMatrix();
    Console.WriteLine($"W    = {W}");
    
    var V = Babai(B, W);
    var Vr = BabaiRound(B, W);
    Console.WriteLine($"V    = {V}");
    Console.WriteLine($"diff = {(V - W).Select(e => e.Absolute).ToKMatrix()}");
    Console.WriteLine($"Vr   = {Vr}");
    Console.WriteLine($"diff = {(Vr - W).Select(e => e.Absolute).ToKMatrix()}");
    Console.WriteLine();
}

{
    var o = Rational.KOne();
    var B = new[] { 1, 2, 3, 3, 0, -3, 3, -7, 3 }.Select(e => e * o).ToKMatrix(3);
    Console.WriteLine("B");
    Console.WriteLine(B);
    
    var W = new[] { 10, 6, 5 }.Select(e => e * o).ToKMatrix();
    Console.WriteLine($"W    = {W}");
    
    var V = Babai(B, W);
    var Vr = BabaiRound(B, W);
    Console.WriteLine($"V    = {V}");
    Console.WriteLine($"diff = {(V - W).Select(e => e.Absolute).ToKMatrix()}");
    Console.WriteLine($"Vr   = {Vr}");
    Console.WriteLine($"diff = {(Vr - W).Select(e => e.Absolute).ToKMatrix()}");
    Console.WriteLine();
}