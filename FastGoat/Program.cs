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

{
    for (int n = 2; n < 100; n++)
    {
        var subgr = PowModSubGroups(n);
        Console.WriteLine($"N:{n,3} [{subgr.Glue(" ")}]");
        var un = new Un(n);
        Console.WriteLine("{0,5} [{1}]", un, un.LongestCycles.Select(a => a.Key[un.Cn[1]].K).Glue(" "));
        Console.WriteLine("{0,5} [{1}]", " ", un.PseudoGenerators.Select(e => e[un.Cn[1]].K).Glue(" "));
        Console.WriteLine();
    }
}