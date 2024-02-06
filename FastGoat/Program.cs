using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using System.Text.RegularExpressions;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void FactorWord(string word)
{
    var reg = new Regex(@"(?'rel'(?:\w|[a-zA-Z]\-\d+)+?)(\k'rel')+");
    var word2 = StringExt.ReducedWordForm1(word);
    Console.WriteLine($"Word:{word} => {word2}");
    var word3 = $"{word2}";
    foreach (Match match in reg.Matches(word2))
    {
        Console.WriteLine($"New group:{match}");
        var str = match.Groups["rel"]!.Value;
        var pat = match.Groups[1]!.Value;
        var exp = str.Length / pat.Length;
        
        foreach (GroupRegX g in match.Groups)
            Console.WriteLine($"group:{g} x {match.Length / g.Value.Length}");

        Console.WriteLine();
    }
}

string ExpandWord(string word)
{
    var regX = new Regex(@"(?:\()([a-zA-Z0-9\-]*)(?:\))((\-{1}\d{1,})|(\d{0,}))");
    var curr = "";
    var succ = $"{word}";
    do
    {
        curr = $"{succ}";
        foreach (Match m in regX.Matches(curr))
        {
            var str = m.Groups[0]!.Value;
            var pat = m.Groups[1]!.Value;
            var exp = int.Parse(m.Groups[2]!.Value);
            if (exp > 0)
            {
                var ext = Enumerable.Repeat(pat, exp).Glue();
                succ = curr.Replace(str, ext);
                break;
            }
            else
            {
                var patInv = StringExt.Revert(pat).Glue();
                var ext = Enumerable.Repeat(patInv, -exp).Glue();
                succ = curr.Replace(str, ext);
                break;
            }
        }
    } while (!string.Equals(succ, curr));

    return StringExt.ReducedWordForm2(succ);
}

{
    var str1 = "ababababababab";
    var str2 = "abab-1abab-1abab-1abab-1";
    var str3 = "abababab-1ab-1ab-1abababab-1ab-1ab-1";
    var str4 = "a2b2a2b2a2b2a2b2acdacdacdabababab";
    var str5 = "dabcabcd";
    var str6 = "aaabca2bc";
    foreach (var str in new[] { str1, str2, str3, str4, str5, str6 })
    {
        FactorWord(str);
    }
    
    // var str5 = "(ab)-2(cde)3(ab)-2";
    // var str6 = "ab2c(efg)-2efg";
    // var str7 = "ab(a(ac)2b)-3";
    //
    // foreach (var str in new[] { str5, str6, str7 })
    // {
    //     Console.WriteLine($"{str} => {ExpandWord(str)}");
    // }
}
