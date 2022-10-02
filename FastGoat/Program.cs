using FastGoat;
using FastGoat.Gp;
using FastGoat.UserGroup;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

for (int k = 1; k < 5; ++k)
{
    GlobalStopWatch.Restart();
    var s7 = new Sn(7); // Trivial s7
    var a = s7[(1, 2, 3, 4, 5, 6, 7)];
    var S7 = Group.GenerateElements(s7, s7[(1, 2)], a);
    var allC3 = S7.Where(p => (p ^ 3) == s7.Neutral()).ToArray();
    var b = allC3.First(p => (a ^ 2) == p * a * (p ^ -1));
    GlobalStopWatch.Stop();

    Console.WriteLine("|S7|={0}, |{{b in S7 with b^3 = 1}}| = {1}",S7.Count(), allC3.Count());
    Console.WriteLine("First Solution a^7 = b^3 = 1 and a^2 = b * a * b^-1 : a = {0} and b = {1}", a, b);
    Console.WriteLine();

    var h = Group.Generate("H", a);
    var g21 = Group.Generate(a, b);
    DisplayGroup.Head(g21);
    DisplayGroup.Head(g21.Over(h));
    GlobalStopWatch.Show("Group21");
}