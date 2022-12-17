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

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Random rnd = new();

{
    var p = 3;
    var o = 3;
    var rg = p.Pow(o).Range();

    foreach (var e1 in rg)
    {
        var z1 = new ZnBInt(p, e1, o);
        var pa1 = new Padic(p, o, e1);
        var znb1 = z1.K;
        var panb1 = pa1.ToBigInteger();

        Console.WriteLine($"{e1} => {z1.PadicNumericString()} {pa1}");

        if (panb1 != znb1)
        {
            Console.WriteLine("#####################");
            throw new Exception();
        }

        if (z1.Opp().K != pa1.Opp().ToBigInteger())
        {
            Console.WriteLine("######### OPP ##########");
            throw new Exception();
        }

        if (z1.K % p != 0 && z1.Inv().K != pa1.Inv().ToBigInteger())
        {
            Console.WriteLine("######### INV ##########");
            throw new Exception();
        }
    }

    foreach (var (e1, e2) in rg.Grid2D(rg))
    {
        var z1 = new ZnBInt(p, e1, o);
        var z2 = new ZnBInt(p, e2, o);
        var pa1 = new Padic(p, o, e1);
        var pa2 = new Padic(p, o, e2);

        var z3 = z1.Add(z2);
        var z4 = z1.Mul(z2);

        var pa3 = pa1.Add(pa2);
        var pa4 = pa1.Mul(pa2);

        if (z3.K != pa3.ToBigInteger())
        {
            Console.WriteLine("######### ADD ##########");
            throw new Exception();
        }

        if (z4.K != pa4.ToBigInteger())
        {
            Console.WriteLine("######### Mul ##########");
            throw new Exception();
        }
    }

    Console.WriteLine("######### PASS ##########");
    Console.WriteLine(new ZnBInt.Infos(p, o));
}
