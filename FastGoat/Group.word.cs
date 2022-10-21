using System.Text.RegularExpressions;
using FastGoat.UserGroup;

namespace FastGoat;

public static partial class Group
{
    static Regex regX = new Regex(@"([a-zA-Z])((\-{1}\d{1,})|(\d{0,}))");

    static (char c, int pow) ToLowerForm(char c, int p)
    {
        var c0 = char.IsLower(c) ? c : char.ToLower(c);
        var p0 = c == c0 ? p : -p;
        return (c0, p0);
    }

    static char Revert(char c) => char.IsLower(c) ? char.ToUpper(c) : char.ToLower(c);
    public static IEnumerable<char> Revert(this IEnumerable<char> letters) => letters.Reverse().Select(Revert);

    public static IEnumerable<char> Add(this IEnumerable<char> w0, IEnumerable<char> w1)
    {
        var s = new Stack<char>(w0);
        foreach (var c in w1)
        {
            if (s.Count == 0)
            {
                s.Push(c);
            }
            else
            {
                var c0 = s.Peek();
                if ((char.IsLower(c) && char.IsUpper(c0) && c0 == char.ToUpper(c)) ||
                    (char.IsUpper(c) && char.IsLower(c0) && c0 == char.ToLower(c)))
                {
                    s.Pop();
                }
                else
                {
                    s.Push(c);
                }
            }
        }

        return s.Reverse();
    }

    static IEnumerable<(char c, int pow)> Reduce(IEnumerable<(char c, int pow)> letters)
    {
        Stack<(char c, int pow)> stack = new();
        foreach (var l in letters)
        {
            if (stack.Count == 0)
            {
                stack.Push(l);
                continue;
            }
            else if (l.c == stack.Peek().c)
            {
                var l0 = stack.Pop();
                var p = l.pow + l0.pow;
                if (p != 0)
                    stack.Push((l.c, p));
            }
            else
                stack.Push(l);
        }

        return stack.Reverse();
    }

    static IEnumerable<(char c, int pow)> ParseReducedWord(string word)
    {
        List<(char c, int pow)> letters = new();
        foreach (Match m in regX.Matches(word))
        {
            var powStr = m.Groups[2].Value;
            var c = char.Parse(m.Groups[1].Value);
            var p = string.IsNullOrEmpty(powStr) ? 1 : int.Parse(powStr);
            letters.Add(ToLowerForm(c, p));
        }

        return Reduce(letters);
    }

    static string ExtendLetters((char c, int pow) e)
    {
        var c0 = e.pow < 0 ? char.ToUpper(e.c) : e.c;
        var p0 = e.pow < 0 ? -e.pow : e.pow;
        return String.Join("", Enumerable.Repeat(c0, p0));
    }

    static string OneLetter((char c, int pow) e)
    {
        return e.pow == 1 ? $"{e.c}" : $"{e.c}{e.pow}";
    }

    static string ExtendLetters(IEnumerable<(char c, int pow)> letters) => letters.Select(ExtendLetters).Glue();
    public static string ReducedWordForm1(string word) => ParseReducedWord(word).Select(OneLetter).Glue();
    public static string ReducedWordForm2(string word) => ExtendLetters(ParseReducedWord(word));

    public static string ExpandRelator(string r)
    {
        var nbEq = r.Count(c => c == '=');
        if (nbEq == 0)
        {
            return ReducedWordForm2(r);
        }
        else if(nbEq == 1)
        {
            var sp = r.Split('=');
            var r1 = ReducedWordForm2(sp[0]);
            var r2 = ReducedWordForm2(sp[1]);
            var rf = r1.Add(r2.Revert());
            return rf.Glue();
        }
        
        throw new GroupException(GroupExceptionType.GroupDef);
    }

    public static ConcreteGroup<Word> Words(string relators)
    {
        return new ConcreteGroup<Word>(new WordGroup(relators));
    }
    public static ConcreteGroup<Word> Words(string name, string relators)
    {
        return new ConcreteGroup<Word>(name, new WordGroup(relators));
    }
}