using System.Collections;
using FastGoat;
using FastGoat.Examples;
using FastGoat.Gp;
using FastGoat.UserGroup;
using static FastGoat.IntExt;
using static FastGoat.EnumerableExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

IEnumerable<IEnumerable<int>> YieldAllPermutationsOld(int n)
{
    if (n == 0)
        yield return Enumerable.Empty<int>();
    else
        foreach (var perm in YieldAllPermutations(n - 1))
            for (int i = 0; i < n; i++)
                yield return perm.Take(i).Append(n).Concat(perm.Skip(i));
}

{
    for (int i = 0; i < 5; i++)
    {
        
        {
            GlobalStopWatch.Restart();
            Console.WriteLine(YieldAllPermutations(9).Count());
            GlobalStopWatch.Show("Yield New");;
        }

        {
            GlobalStopWatch.Restart();
            Console.WriteLine(YieldAllPermutationsOld(9).Count());
            GlobalStopWatch.Show("Yield Old");
        }

    }
}