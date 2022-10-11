using FastGoat;
using FastGoat.Examples;
using FastGoat.Gp;
using FastGoat.UserGroup;
using static FastGoat.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    for (int n = 3; n < 38; ++n)
    {
        var carm = Carmichael(n);
        var phi = Phi(n);
        Console.WriteLine("N:{0,-2} Phi:{1,-2} L:{2,-2} [{3}]", n,phi,  carm.Min(), carm.Glue(" "));
    }

}