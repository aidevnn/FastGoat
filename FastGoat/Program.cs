using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
RngSeed(258);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

{
    var (n, t0) = (8, 17);
    var (pm, sk, t, q, pk, rlk) = RLWE.KeyGenBGV(n, t0, t0 * 19 * 23 * 29);
    
    var m1 = RLWE.GenUnif(n, t);
    var c1 = RLWE.EncryptBGV(m1, pm, t, q, pk);
    c1.Show("c1");
    Console.WriteLine();

    var c2 = c1.ModSwitch(17 * 19, t);
    c2.Show("c2");
    Console.WriteLine();

    var c3 = c1.ModSwitch(17 * 19 * 23, t);
    c3.Show("c3");
    Console.WriteLine();

    var c4 = c3.ModSwitch(17 * 19, t);
    c4.Show("c4");
    Console.WriteLine();

    var c5 = c3.ModSwitch(17 * 23, t);
    c5.Show("c5");
    Console.WriteLine();

    var c6 = c5.ModSwitch(17 * 19 * 23, t);
    c6.Show("c6");
    Console.WriteLine();

    Console.WriteLine($"m1:{m1}");
    Console.WriteLine($"d1:{RLWE.DecryptBGV(c1, pm, sk, t)}");
    Console.WriteLine($"d2:{RLWE.DecryptBGV(c2, pm, sk, t)}");
    Console.WriteLine($"d3:{RLWE.DecryptBGV(c3, pm, sk, t)}");
    Console.WriteLine($"d4:{RLWE.DecryptBGV(c4, pm, sk, t)}");
    Console.WriteLine($"d5:{RLWE.DecryptBGV(c5, pm, sk, t)}");
    Console.WriteLine($"d6:{RLWE.DecryptBGV(c6, pm, sk, t)}");
}