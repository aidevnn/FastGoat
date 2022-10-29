using System.Collections;
using FastGoat.Commons;

namespace FastGoat.ToddCoxeter;

public class Header : IEnumerable<Generator>
{
    List<Generator> head { get; }

    public Header(IEnumerable<IEnumerable<Generator>> gens)
    {
        List<int> seps = new() { 0 };
        var arrGens = gens.ToArray();
        foreach (var g in arrGens)
            seps.Add(seps.Last() + g.Count());

        Separators = seps.ToArray();
        head = arrGens.SelectMany(g => g).ToList();
        Generators = head.Union(head.Select(g => g.Invert())).ToHashSet();
        Count = head.Count;
    }

    public Header(Header header)
    {
        head = header.head.ToList();
        Separators = header.Separators.ToArray();
        Generators = head.ToHashSet();
        Count = head.Count;
    }

    public int[] Separators { get; }
    public int Count { get; }
    public HashSet<Generator> Generators { get; }
    public Generator this[int k] => head[k];
    public IEnumerator<Generator> GetEnumerator() => head.GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => head.GetEnumerator();

    public string Display(int digits)
    {
        var fmt = $"{{0,{digits + 1}}}";
        return " " + head.Glue("", fmt);
    }

    public void DisplayLineUp(int digits)
    {
        var s1 = Enumerable.Repeat('─', (digits + 1) * (head.Count + 1)).Append(' ').Glue().ToArray();
        foreach (var k in Separators)
            s1[(digits + 1) * k] = s1[(digits + 1) * (k + 1)] = '┴';

        s1[0] = '└';
        s1[^1] = '┘';
        Console.WriteLine(s1.Glue());
    }

    public void DisplayLineDown(int digits)
    {
        var s1 = Enumerable.Repeat('─', (digits + 1) * (head.Count + 1)).Append(' ').Glue().ToArray();
        foreach (var k in Separators)
            s1[(digits + 1) * k] = s1[(digits + 1) * (k + 1)] = '┬';

        s1[0] = '┌';
        s1[^1] = '┐';
        Console.WriteLine(s1.Glue());
    }

    public void DisplayHead(int digits)
    {
        Console.WriteLine(Display(digits));
        DisplayLineDown(digits);
    }

    public void ReDisplayHead(int digits)
    {
        DisplayLineUp(digits);
        DisplayHead(digits);
    }

    public override string ToString() => "  " + head.Glue(" ");
}