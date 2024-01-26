using System.Linq;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.Tools;
using Xunit;

namespace Tests;

public class ToddCoxeterUnitTest
{
    [Fact]
    public void Test1Cyclic()
    {
        var ta5 = Graph.Run("a5");
        var gens = ta5.Generators().ToArray();
        var words = ta5.Words().Select(e => e.Glue()).ToArray();
        Assert.Single(gens);
        Assert.Equal('a', gens[0]);
        Assert.Contains(new[] { "", "a", "aa", "aaa", "A" }, words.Contains);
        Assert.Equal("aa", ta5.Rewrite("aaaAaaaAaaa"));
    }

    [Fact]
    public void Test2Klein()
    {
        var klein = Graph.Run("a2, b2, ab = ba");
        var gens = klein.Generators().ToArray();
        var words = klein.Words().Select(e => e.Glue()).ToArray();
        Assert.Contains(new[] { 'a', 'b' }, gens.Contains);
        Assert.Contains(new[] { "", "a", "b", "ab" }, words.Contains);
        Assert.Equal("ab", klein.Rewrite("ababAbbaAb").Glue());
    }

    [Fact]
    public void Test3H21()
    {
        var h21 = Graph.Run("a7, b3, a2=bab-1");
        var gens = h21.Generators().ToArray();
        var words = h21.Words().Select(e => e.Glue()).ToArray();
        Assert.Contains(new[] { 'a', 'b' }, gens.Contains);
        Assert.Equal(21, words.Length);
    }

    [Fact]
    public void Test4Abelian()
    {
        var ab = Graph.Run("a4, b3, c2, ab=ba, ac=ca, bc=cb");
        var gens = ab.Generators().ToArray();
        var words = ab.Words().Select(e => e.Glue()).ToArray();
        Assert.Contains(new[] { 'a', 'b', 'c' }, gens.Contains);
        Assert.Equal(24, words.Length);
    }

    [Fact]
    public void Test5F20()
    {
        var f20 = Graph.Run("a5, b4, abababab, a2ba-1b-1");
        var gens = f20.Generators().ToArray();
        var words = f20.Words().Select(e => e.Glue()).ToArray();
        Assert.Contains(new[] { 'a', 'b' }, gens.Contains);
        Assert.Equal(20, words.Length);
    }

    [Fact]
    public void Test6WordGroup()
    {
        var wg = new WordGroup("a12, b2=a6, bab-1=a-1");
        var w = wg["bababa"];
        Assert.Equal(24, wg.Count());
        Assert.Equal(GroupType.NonAbelianGroup, wg.GroupType);
        Assert.Equal("AB", w.Get());
    }
}