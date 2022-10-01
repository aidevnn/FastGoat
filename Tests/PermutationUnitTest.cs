using System.Linq;
using FastGoat;
using Xunit;

namespace Tests;

public class PermutationUnitTest
{
    [Fact]
    public void Test1Tuple2Array()
    {
        Tuple2Array c1 = 1;
        Tuple2Array c2 = (2, 1);
        Tuple2Array c3 = (2, 3, 5);
        Tuple2Array c4 = (2, 3, 4, 5);
        Tuple2Array c5 = (8, 6, 2, 3, 5);
        Tuple2Array c6 = (2, 9, 3, 1, 5, 4);
        Tuple2Array c7 = (2, 8, 6, 11, 3, 5, 0);
        Tuple2Array c8 = (2, 8, 6, 11, 3, 10, 5, 8, 0);
        Tuple2Array c9 = (2, (1, 7), 5);

        Assert.Single(c1.Table);
        Assert.Equal(2, c2.Table.Length);
        Assert.Equal(3, c3.Table.Length);
        Assert.Equal(4, c4.Table.Length);
        Assert.Equal(5, c5.Table.Length);
        Assert.Equal(6, c6.Table.Length);
        Assert.Equal(7, c7.Table.Length);
        Assert.Empty(c8.Table);
        Assert.Empty(c9.Table);
    }

    [Fact]
    public void Test2Permutation2Cycles()
    {
        var p1 = IntExt.PermutationToCycles(2, new[] { 0, 1 });
        Assert.Equal(2, p1.Length);
        Assert.Equal(0, p1[0][0]);
        Assert.Equal(1, p1[1][0]);

        var p2 = IntExt.PermutationToCycles(2, new[] { 1, 0 });
        Assert.Single(p2);
        Assert.Equal(0, p2[0][0]);
        Assert.Equal(1, p2[0][1]);

        var p3 = IntExt.PermutationToCycles(5, new[] { 2, 0, 1, 3, 4 });
        Assert.True(p3[0].SequenceEqual(new[] { 0, 2, 1 }));
        Assert.True(p3[1].SequenceEqual(new[] { 3 }));
        Assert.True(p3[2].SequenceEqual(new[] { 4 }));

        var p4 = IntExt.PermutationToCycles(5, new[] { 2, 3, 4, 1, 0 });
        Assert.True(p4[0].SequenceEqual(new[] { 0, 2, 4 }));
        Assert.True(p4[1].SequenceEqual(new[] { 1, 3 }));

        var p5 = IntExt.PermutationToCycles(5, new[] { 2, 3, 4, 3, 0 });
        Assert.Empty(p5);

        var p6 = IntExt.PermutationToCycles(0, new[] { 2, 3, 4, 3, 0 });
        Assert.Empty(p6);
    }

    [Fact]
    public void Test3PermutationsOperations()
    {
        var p1 = new[] { 1, 0, 3, 2 };
        var p2 = new[] { 2, 0, 1, 3 };
        var p3 = new[] { 1, 2, 0, 3 };
        var p4 = new[] { 0, 2, 3, 1 };
        var c = new[] { 0, 0, 0, 0 };

        IntExt.InvertPermutation(p1, c);
        Assert.True(p1.SequenceEqual(c));

        IntExt.InvertPermutation(p2, c);
        Assert.True(p3.SequenceEqual(c));

        IntExt.ComposePermutation(p1, p2, c);
        Assert.True(p4.SequenceEqual(c));
        
        c = new[] { 0, 1, 2, 3 };
        IntExt.ApplyCycle(c, new[] { 2, 1, 0 });
        Assert.True(p2.SequenceEqual(c));
        IntExt.ApplyCycle(c, new[] { 0, 1 });
        IntExt.ApplyCycle(c, new[] { 3, 2 });
        Assert.True(p4.SequenceEqual(c));
    }

    [Fact]
    public void Test4AllPermutations()
    {
        Assert.Equal(6,IntExt.GetPermutations(3).Length);
        Assert.Equal(24,IntExt.GetPermutations(4).Length);
        Assert.Equal(120,IntExt.GetPermutations(5).Length);
        Assert.Equal(720,IntExt.GetPermutations(6).Length);
        Assert.Equal(5040,IntExt.GetPermutations(7).Length);
        Assert.Equal(40320,IntExt.GetPermutations(8).Length);
    }
}