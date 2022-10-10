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
    var s5 = Group.Create(new Sn(5));
    var card = 12; // 10, 20,  
    var allC12 = s5.SelectMany(a => s5.Select(b => (a, b))).Where(e => Group.Generate(s5, e.a, e.b).Count() == card);
    foreach (var (a,b) in allC12)
    {
        Console.WriteLine($"a={a}; b={b}");
        DisplayGroup.HeadElements(Group.Generate("M12", a, b));
    }
}