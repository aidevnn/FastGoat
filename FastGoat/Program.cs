using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Lattice;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    var lweRegev = new Regev(10);
    Console.WriteLine(lweRegev.Params);
    var bits = DistributionExt.BernouilliSample(100, 0.5).ToArray();
    var encrypted = bits.Select(e => lweRegev.EncryptBit(e)).ToArray();
    var decrypted = encrypted.Select(e => lweRegev.DecryptBit(e)).ToArray();
    Console.WriteLine($"bits     :[{bits.Glue(", ")}]");
    Console.WriteLine($"decrypted:[{decrypted.Glue(", ")}]");
    Console.WriteLine($"Are equal:{bits.SequenceEqual(decrypted)}");
    lweRegev.Err.GroupBy(e => e.Signed).Select(e => new { e = e.Key, Nb = e.Count() })
        .OrderBy(e => e.e).Println("Err sample");
}